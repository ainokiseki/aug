#pragma once
#include"tri.h"
//#include<vector>
//#include<thread>
class MatrixReplacement;
using Eigen::SparseMatrix;

namespace Eigen {
    namespace internal {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<>
        struct traits<MatrixReplacement> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
        {};
    }
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };

    Index rows() const { return var_num; }
    Index cols() const { return var_num; }

    template<typename Rhs>
    Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
        return Eigen::Product<MatrixReplacement, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
    }

    // Custom API:
    MatrixReplacement(std::vector<int64_t>* aj1, int n, double a, Eigen::VectorXd* aj2, double sigma) :aj1(aj1), n(n), a(a), sigma(sigma), aj2(aj2) {
        var_num = n * (n - 1) / 2;

    }

    std::vector<int64_t>* aj1;
    int var_num;
    double a;
    int n;
    double sigma;
    Eigen::VectorXd* aj2;
};


// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
    namespace internal {

        template<typename Rhs>
        struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
            : generic_product_impl_base<MatrixReplacement, Rhs, generic_product_impl<MatrixReplacement, Rhs> >
        {
            typedef typename Product<MatrixReplacement, Rhs>::Scalar Scalar;

            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha == Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);

                // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
                // but let's do something fancier (and less efficient):
                Eigen::VectorXd ap = Eigen::VectorXd::Zero(rhs.rows());
                if (lhs.aj1->size() > 120000) {

                    //int thnum = 4;
                    int thnum = std::fmin(lhs.aj1->size()/120000+1,4);
                    Eigen::VectorXd* apres = new Eigen::VectorXd[thnum];
                    for (int i = 0; i < thnum; i++)
                        apres[i] = Eigen::VectorXd::Zero(lhs.var_num);
                    Eigen::VectorXd rhscopy(lhs.var_num);
                    for (int i = 0; i < lhs.var_num; i++)
                        rhscopy[i] = rhs[i];
                    std::vector<std::thread>tv;
                    for (int i = 0; i < thnum; i++)
                        tv.push_back(std::thread(cacap, apres + i, i * lhs.aj1->size() / thnum, (i + 1) * lhs.aj1->size() / thnum, lhs.aj1, &rhscopy,lhs.n));

                    for (int i = 0; i < thnum; i++)
                        tv[i].join();

                    for (int i = 0; i < thnum; i++)
                        ap += apres[i];
                    delete[] apres;
                }
                else {

                    int n = lhs.n;
                    //Eigen::VectorXd ap = Eigen::VectorXd::Zero(rhs.rows());
                    //ij,ik,jk
                    for (int i = 0; i < lhs.aj1->size(); i++) {
                        int64_t tem = (*(lhs.aj1))[i];
                        int ci, cj, ck;
                        ck = tem % n;
                        cj = tem / n;
                        ci = cj / n;
                        cj = cj % n;
                        int idx[3];
                        idx[0] = ij(ci, cj, n);
                        idx[1] = ij(ci, ck, n);
                        idx[2] = ij(cj, ck, n);


                        double num =rhs[idx[0]] - rhs[idx[1]]- rhs[idx[2]];

                        ap[idx[0]] += num;

                        ap[idx[1]] -= num;

                        ap[idx[2]] -= num;

                    }
                }
                dst = ap * lhs.sigma + lhs.sigma * lhs.aj2->cwiseProduct(rhs) + lhs.a * rhs;

                //int thnum = 4;
                //Eigen::VectorXd* ap = new  Eigen::VectorXd[thnum];
                //vector<thread>v; 
                //double* drhs = new double[lhs.var_num];
                //for (int i = 0; i < lhs.var_num; i++)
                //    drhs[i] = rhs[i];
                //multi_calc_ap mca = multi_calc_ap( lhs.aj1,drhs);
                //int l = lhs.aj1->size();
                //for (int i = 0; i < thnum; i++) {
                //    ap[i] = Eigen::VectorXd::Zero(lhs.var_num);
                //    v.push_back(thread(mca, l * i / thnum, l * (i + 1) / thnum, ap + i));
                //}
                //for (int i = 0; i < thnum; i++)
                //    v[i].join();
                //for (int i = 1; i < thnum; i++)
                //    ap[0] += ap[i];

                //dst = ap[0]*lhs.sigma + lhs.sigma * (lhs.vaj2.array() * rhs.array()).matrix() +lhs.a*rhs;
                //delete[] ap;
                //delete[]drhs;


            }
        };

    }
}

