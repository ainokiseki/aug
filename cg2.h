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
	MatrixReplacement2(std::vector<tri*>* aj1, int var_num, double a, Eigen::VectorXd* x,Eigen::VectorXd *w, double sigma) :aj1(aj1), var_num(var_num), a(a), sigma(sigma), x(x), w(w) {


	}

	std::vector<tri*>* aj1;
	int var_num;
	double a;

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

			}
		};

	}
}