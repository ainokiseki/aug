#pragma once
#include<Eigen/Dense>
#include<map>
#include"func_eigen.h"
class axbu_eigen
{
public:
	axbu_eigen(int n,Eigen::VectorXd*x, Eigen::VectorXd* b, std::map<long long,double>*u ,double sigma);
	~axbu_eigen();
	void operator()();
	void set_function(func_eigen* );


private:
	int n;
	Eigen::VectorXd* x;
	Eigen::VectorXd* b;
	std::map<long long, double>* u;
	double sigma;
	func_eigen* f;
};



