#include "newaxbu.h"
//constexpr int scale_a = 3;


void new_axbu(almfunc& f, concurrency::concurrent_vector<cons_full>::iterator beg, concurrency::concurrent_vector<cons_full>::iterator end, int n, Eigen::VectorXd& x, Eigen::VectorXd& b, double sigma) {
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (auto i = beg; i != end; i++) {
		int* idx = i->idx;
		double axbu = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]])/ scale_a + i->u / sigma;
		f(axbu, idx);
	}
	f.integrate();
}
