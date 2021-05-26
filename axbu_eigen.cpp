#include "axbu_eigen.h"
#include<iostream>

axbu_eigen::axbu_eigen(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<int64_t, double>* u, double sigma) :
	n(n), x(x), b(b), u(u), sigma(sigma),f(NULL)
{
}

axbu_eigen::~axbu_eigen()
{

}
void axbu_eigen::operator() (int init,int end){
	if (end < 0)end = n;
	int var_num = n * (n - 1) / 2;
	
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int i = init; i < end; i++) {
		for (int j = i+1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;
				double axbu = 0;
				for (int time = 0; time < 3; time++) {
					axbu -= (*x)[idx[time]] + (*b)[idx[time]];
				}
				for (int time = 0; time < 3; time++) {
					int64_t index = longindex[time];
					double axbu2 = axbu;
					if (u->find(index) != u->end()) {
						axbu2 += (*u)[index] / sigma;
					}
					axbu2 += 2 * ((*x)[idx[time]] + (*b)[idx[time]]);
					(*f)(axbu2, idx, time, index);
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

au_eigen::au_eigen(int n, std::map<int64_t, double>* u, Eigen::VectorXd v, Eigen::VectorXd w) :n(n), u(u), v(v),w(w)
{
	int var_num = n * (n - 1) / 2;
	au = Eigen::VectorXd::Zero(var_num);
}

au_eigen::~au_eigen()
{
}

double au_eigen::operator()(int init, int end)
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
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;
				for (int time = 0; time < 3; time++) {
					int64_t index = longindex[time];
					if (u->find(index) == u->end())
						continue;
					double tem = (*u)[index];
					(au)[idx[0]] -= tem;
					(au)[idx[1]] -= tem;
					(au)[idx[2]] -= tem;
					(au)[idx[time]] += 2 * tem;

				}
			}
		}

	}
	return (au + w.cwiseProduct(v)).norm();
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

axbu_aj1::axbu_aj1(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<int64_t, double>* u, double sigma, Eigen::SparseMatrix<double>* mat,Eigen::VectorXd*rhs):
	n(n), x(x), b(b), u(u), sigma(sigma), mat(mat),rhs(rhs)
{
}

axbu_aj1::~axbu_aj1()
{
}

void axbu_aj1::operator()(int init, int end)
{
	if (end < 0)end = n;
	int var_num = n * (n - 1) / 2;
	Eigen::VectorXd diag=Eigen::VectorXd::Zero(var_num);
	std::vector<Eigen::Triplet<double> >triv;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int i = init; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;
				double axbu = 0;

				axbu = (*x)[idx[0]] + (*b)[idx[0]] - (*x)[idx[1]] - (*b)[idx[1]] - (*x)[idx[2]] - (*b)[idx[2]];

				if (u->find(longindex[0]) != u->end()) {
					axbu += (*u)[longindex[0]] / sigma;
				}
				if (axbu > 0) {
					(*rhs)(idx[0]) += axbu;
					(*rhs)(idx[1]) -= axbu;
					(*rhs)(idx[2]) -= axbu;

					diag[idx[0]] += 1;
					diag[idx[1]] += 1;
					diag[idx[2]] += 1;

					triv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
					triv.push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
				}

				axbu = (*x)[idx[1]] + (*b)[idx[1]] - (*x)[idx[0]] - (*b)[idx[0]] - (*x)[idx[2]] - (*b)[idx[2]];

				if (u->find(longindex[1]) != u->end()) {
					axbu += (*u)[longindex[1]] / sigma;
				}
				if (axbu > 0) {
					(*rhs)(idx[1]) += axbu;
					(*rhs)(idx[0]) -= axbu;
					(*rhs)(idx[2]) -= axbu;

					diag[idx[0]] += 1;
					diag[idx[1]] += 1;
					diag[idx[2]] += 1;

					triv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[0], idx[2], 1.0));
					triv.push_back(Eigen::Triplet<double>(idx[2], idx[0], 1.0));
					triv.push_back(Eigen::Triplet<double>(idx[1], idx[2], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[2], idx[1], -1.0));
				}
				axbu = (*x)[idx[2]] + (*b)[idx[2]] - (*x)[idx[1]] - (*b)[idx[1]] - (*x)[idx[0]] - (*b)[idx[0]];

				if (u->find(longindex[2]) != u->end()) {
					axbu += (*u)[longindex[2]] / sigma;
				}
				if (axbu > 0) {
					(*rhs)(idx[2]) += axbu;
					(*rhs)(idx[1]) -= axbu;
					(*rhs)(idx[0]) -= axbu;

					diag[idx[0]] += 1;
					diag[idx[1]] += 1;
					diag[idx[2]] += 1;

					triv.push_back(Eigen::Triplet<double>(idx[0], idx[1], 1.0));
					triv.push_back(Eigen::Triplet<double>(idx[1], idx[0], 1.0));
					triv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[1], idx[2], -1.0));
					triv.push_back(Eigen::Triplet<double>(idx[2], idx[1], -1.0));
				}
			}
		}

	}
	mat->setFromTriplets(triv.begin(), triv.end());
	(*mat) += diag.asDiagonal();


	return;
}
