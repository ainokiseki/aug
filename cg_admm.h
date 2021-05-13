#pragma once
#include"tri.h"
#include<vector>
#include<Eigen/Sparse>
class MatrixReplacement_admm;
using Eigen::SparseMatrix;

namespace Eigen {
	namespace internal {
		// MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
		template<>
		struct traits<MatrixReplacement_admm> : public Eigen::internal::traits<Eigen::SparseMatrix<double> >
		{};
	}
}

// Example of a matrix-free wrapper from a user type to Eigen's compatible type
// For the sake of simplicity, this example simply wrap a Eigen::SparseMatrix.
class MatrixReplacement_admm : public Eigen::EigenBase<MatrixReplacement_admm> {
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
	Eigen::Product<MatrixReplacement_admm, Rhs, Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
		return Eigen::Product<MatrixReplacement_admm, Rhs, Eigen::AliasFreeProduct>(*this, x.derived());
	}

	// Custom API:
	MatrixReplacement_admm(int n,Eigen::VectorXd* w) :n(n),w(w) {
		var_num = n * (n - 1) / 2;

	}


	int var_num;
	int n;
	Eigen::VectorXd* w;
};


namespace Eigen {
	namespace internal {

		template<typename Rhs>
		struct generic_product_impl<MatrixReplacement_admm, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
			: generic_product_impl_base<MatrixReplacement_admm, Rhs, generic_product_impl<MatrixReplacement_admm, Rhs> >
		{
			typedef typename Product<MatrixReplacement_admm, Rhs>::Scalar Scalar;

			template<typename Dest>
			static void scaleAndAddTo(Dest& dst, const MatrixReplacement_admm& lhs, const Rhs& rhs, const Scalar& alpha)
			{
				// This method should implement "dst += alpha * lhs * rhs" inplace,
				// however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
				assert(alpha == Scalar(1) && "scaling is not implemented");
				EIGEN_ONLY_USED_FOR_DEBUG(alpha);

				// Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
				// but let's do something fancier (and less efficient):

				Eigen::VectorXd ap = Eigen::VectorXd::Zero(rhs.rows());
				auto ij = [](int64_t x, int64_t y, int64_t n) {
					int64_t res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
					return res;
				};
				for (int i = 0; i <lhs.n; i++) {
					for (int j = i + 1; j < lhs.n; j++) {
						for (int k = j + 1; k < lhs.n; k++) {
							
							int idx[3];
							idx[0] = ij(i, j, lhs.n);
							idx[1] = ij(i, k, lhs.n);
							idx[2] = ij(j, k, lhs.n);

							for (int time = 0; time < 3; time++) {
								double num = -(rhs[idx[0]] + rhs[idx[1]] + rhs[idx[2]]) + 2 * rhs[idx[time]];

								ap[idx[0]] -= num;

								ap[idx[1]] -= num;

								ap[idx[2]] -= num;

								ap[idx[time]] += 2 * num;
							}
						}
					}

				}
				dst = ap + lhs.w->cwiseProduct(*(lhs.w)).cwiseProduct(rhs);

			}
		};

	}
}