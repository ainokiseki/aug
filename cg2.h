#pragma once
#include"tri.h"
#include<vector>
class MatrixReplacement2;
using Eigen::SparseMatrix;

namespace Eigen {
	namespace internal {
		// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
		template<>
		struct traits<MatrixReplacement2> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
		{};
	}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement2 : public Eigen::EigenBase<MatrixReplacement2> {
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
	Eigen::Product<MatrixReplacement2, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement2, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}

	// Custom API:
	MatrixReplacement2(std::vector<int64_t>* aj1, int n, double a, Eigen::VectorXd* x, Eigen::VectorXd* w, double sigma) :aj1(aj1), n(n), a(a), sigma(sigma), x(x), w(w) {

		var_num = n * (n - 1) / 2;
	}

	std::vector<int64_t>* aj1;
	int var_num;
	double a;
	int n;
	double sigma;
	Eigen::VectorXd* x;
	Eigen::VectorXd* w;
};


namespace Eigen {
	namespace internal {

		template<typename Rhs>
		struct generic_product_impl<MatrixReplacement2, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
			: generic_product_impl_base<MatrixReplacement2, Rhs, generic_product_impl<MatrixReplacement2, Rhs> >
		{
			typedef typename Product<MatrixReplacement2, Rhs>::Scalar Scalar;

			template<typename Dest>
			static void scaleAndAddTo(Dest& dst, const MatrixReplacement2& lhs, const Rhs& rhs, const Scalar& alpha)
			{
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				assert(alpha == Scalar(1) && "scaling is not implemented");
				EIGEN_ONLY_USED_FOR_DEBUG(alpha);

				// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
				// but let's do something fancier (and less efficient):

				Eigen::VectorXd ap = Eigen::VectorXd::Zero(rhs.rows());
				int n = lhs.n;
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


					double num = rhs[idx[0]] - rhs[idx[1]] - rhs[idx[2]];

					ap[idx[0]] += num;

					ap[idx[1]] -= num;

					ap[idx[2]] -= num;

				}
				Eigen::VectorXd ix = lhs.w->cwiseProduct(*(lhs.w)).cwiseProduct(rhs);
				//dst = ap * lhs.sigma + lhs.sigma * lhs.aj2->cwiseProduct(rhs) + lhs.a * rhs;
				if (lhs.x->norm() <= 1 / lhs.sigma) {
					dst = ap * lhs.sigma + lhs.sigma * ix + lhs.a * rhs;
				}
				else {
					double xnorm = lhs.x->norm();

					Eigen::VectorXd xx = lhs.x->dot(lhs.w->cwiseProduct(rhs)) * lhs.w->cwiseProduct(*(lhs.x)) / (xnorm * xnorm);
					dst = ap * lhs.sigma + lhs.a * rhs + 1 / xnorm * (ix - xx);
				}
			}
		};

	}
}