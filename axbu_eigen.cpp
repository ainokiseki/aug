#include "axbu_eigen.h"
#include<iostream>

axbu_eigen::axbu_eigen(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<int64_t, double>* u, double sigma, std::vector<int64_t>* cons, std::vector<int>*colidx) :
	n(n), x(x), b(b), u(u), sigma(sigma),f(NULL),cons(cons),colidx(colidx)
{
}

axbu_eigen::~axbu_eigen()
{

}
void axbu_eigen::operator() (int init,int end){
	if (end < 0)end = cons->size();

	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int istart = init; istart < end; istart++) {
		int64_t longindex = (*cons)[istart];
		int i, j, k;
		k = longindex % n;
		j = longindex / n;
		i = j / n;
		j = j % n;
		//std::cout << i << "\t" << j << "\t" << k <<std::endl;


				int idx[3];
				idx[0] = (*colidx)[ij(i, j, n)];
				idx[1] = (*colidx)[ij(i, k, n)];
				idx[2] = (*colidx)[ij(j, k, n)];


				double axbu = 0;
				for (int time = 0; time < 3; time++) {
					axbu -= (*x)[idx[time]] + (*b)[idx[time]];
				}
				for (int time = 0; time < 3; time++) {
					int64_t index;
					switch (time)
					{
					case 0:
						index = longindex;
						break;
					case 1:
						index = static_cast<int64_t>(i) * n * n + k * n + j;
						break;
					case 2:
						index = static_cast<int64_t>(j) * n * n + k * n + i;
					}

					double axbu2 =0;
					if (u->find(index) != u->end()) {
						axbu2 += (*u)[index] / sigma;
					}
					for (int tt = 0; tt < 3; tt++) {
						if(tt==time)
						axbu2 +=  ((*x)[idx[tt]] + (*b)[idx[tt]]);
						else
							axbu2 -=((*x)[idx[tt]] + (*b)[idx[tt]]);
					}
					//axbu2 += 2 * ((*x)[idx[time]] + (*b)[idx[time]]);
					(*f)(axbu2, idx, time, index);
				}



	}

	delete f;

	return;

}
void axbu_eigen::set_function(func_eigen* func) {

	f = func;
}

au_eigen::au_eigen(int n, std::map<int64_t, double>* u, Eigen::VectorXd v, Eigen::VectorXd w, std::vector<int64_t>* cons, std::vector<int>* colidx) :n(n), u(u), v(v),w(w),cons(cons),colidx(colidx)
{
	int var_num = n * (n - 1) / 2;
}

au_eigen::~au_eigen()
{
}

void au_eigen::operator()(Eigen::VectorXd*au, int init, int end)
{
	if (end < 0)end = cons->size();

	int var_num = n * (n - 1) / 2;
	//int64_t index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int istart = init; istart < end; istart++) {
		int64_t longindex = (*cons)[istart];
		int i, j, k;
		k = longindex % n;
		j = longindex / n;
		i = j / n;
		j = j % n;


				int idx[3];
				idx[0] = (*colidx)[ij(i, j, n)];
				idx[1] = (*colidx)[ij(i, k, n)];
				idx[2] = (*colidx)[ij(j, k, n)];

				for (int time = 0; time < 3; time++) {
					int64_t index ;
					switch (time)
					{
					case 0:
						index = longindex;
						break;
					case 1:
						index = static_cast<int64_t>(i) * n * n + k * n + j;
						break;
					case 2:
						index = static_cast<int64_t>(j) * n * n + k * n + i;
					}

					if (u->find(index) == u->end())
						continue;
					double tem = (*u)[index];
					(*au)[idx[0]] -= tem;
					(*au)[idx[1]] -= tem;
					(*au)[idx[2]] -= tem;
					(*au)[idx[time]] += 2 * tem;

				}


	}
	
}

axbu_pro::axbu_pro(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<int64_t, double>* u, double sigma)
	:n(n),x(x),b(b),u(u),sigma(sigma)
{
	res = Eigen::VectorXd::Zero(n * (n - 1) / 2);
}

axbu_pro::~axbu_pro()
{
}

Eigen::VectorXd axbu_pro::operator()(int init, int end)
{

	if (end < 0)end = n;
	int var_num = n * (n - 1) / 2;
	//int64_t index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int i = init; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				double axbu = 0;
				double bb = 0;
				for (int time = 0; time < 3; time++) {
					axbu -= (*x)[idx[time]] + (*b)[idx[time]];
					bb -= (*b)[idx[time]];
				}

				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;
				for (int time = 0; time < 3; time++) {
					int64_t index = longindex[time];
					double axbu2 = axbu;
					bb += 2 * (*b)[idx[time]];
					if (u->find(index) != u->end()) {
						axbu2 += (*u)[index] / sigma;
					}
					axbu2 += 2 * ((*x)[idx[time]] + (*b)[idx[time]]);
					axbu2 = (axbu2 <= 0) ? axbu2 : 0.0;
					if (u->find(index) != u->end()) {
						axbu2 -= (*u)[index] / sigma;
					}
					axbu2 += bb;
					res[idx[0]] -= axbu2;
					res[idx[1]] -= axbu2;
					res[idx[2]] -= axbu2;
					res[idx[time]] += 2 * axbu2;
				}

			}
		}

	}
	return res;
}
