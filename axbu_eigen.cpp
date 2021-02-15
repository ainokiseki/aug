#include "axbu_eigen.h"
#include<iostream>
inline long long ij(long long x, long long y, long long n) {
	long long res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;

	return res;
}
axbu_eigen::axbu_eigen(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<long long, double>* u, double sigma) :
	n(n), x(x), b(b), u(u), sigma(sigma),f(NULL)
{
}

axbu_eigen::~axbu_eigen()
{

}
void axbu_eigen::operator() (int init,int end){
	if (end < 0)end = n;
	int var_num = n * (n - 1) / 2;
	long long index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int i = init; i < end; i++) {
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