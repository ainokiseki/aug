#include "axbu_eigen.h"


axbu_eigen::axbu_eigen(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<long long, double>* u, double sigma) :
	n(n), x(x), b(b), u(u), sigma(sigma)
{
}

axbu_eigen::~axbu_eigen()
{
}
void axbu_eigen::operator() (){

	int var_num = n * (n - 1) / 2;
	long long index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				double axbu = 0;
				for (int time = 0; time < 3; time++) {
					axbu -= (*x)[idx[time]] + (*b)[idx[time]];
				}
				for (int time = 0; time < 3; time++,index++) {
					double axbu2 = axbu;
					if (u->find(index) != u->end()) {
						axbu2 += (*u)[index] / sigma;
					}
					axbu2 += 2* ((*x)[idx[time]] + (*b)[idx[time]]);
					(*f)(axbu2,idx,time,index);
				}

			}
		}
	}
	delete f;

	return;

}
void axbu_eigen::set_function(func_eigen* func) {
	f = func;
}