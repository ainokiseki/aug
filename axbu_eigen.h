#pragma once
#include<Eigen/Dense>
#include<map>
#include"func_eigen.h"

class axbu_eigen
{
public:
	axbu_eigen(int n,Eigen::VectorXd*x, Eigen::VectorXd* b, std::map<long long,double>*u ,double sigma);
	~axbu_eigen();
	void operator()(int init=0,int end=-1);
	void set_function(func_eigen* );


private:
	int n;
	Eigen::VectorXd* x;
	Eigen::VectorXd* b;
	std::map<long long, double>* u;
	double sigma;
	func_eigen* f;
};
class au_eigen
{
public:
	au_eigen();
	~au_eigen();
	void operator()(int init=0,int end=-1);
private:

};

au_eigen::au_eigen()
{
}

au_eigen::~au_eigen()
{
}

void au_eigen::operator()(int init,int end)
{
	if (end < 0)end = n;
	int var_num = n * (n - 1) / 2;
	long long index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int i = init; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				double axbu = 0;
				for (int time = 0; time < 3; time++) {
					if (u->find(index) == u->end())
						continue;


				}


			}
		}

	}
}


