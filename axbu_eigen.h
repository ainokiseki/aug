#pragma once
#include<Eigen/Dense>
#include<map>
#include"func_eigen.h"
//#include"cg_inf.h"
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
	au_eigen(int, std::map<long long, double>*, Eigen::VectorXd ,Eigen::VectorXd);
	~au_eigen();
	double operator()(int init=0,int end=-1);
private:
	std::map<long long, double>* u;
	int n;
	Eigen::VectorXd v;
	Eigen::VectorXd au;
	Eigen::VectorXd w;
};


class axbu_pro
{
public:
	axbu_pro(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<long long, double>* u, double sigma);
	~axbu_pro();
	Eigen::VectorXd operator()(int init = 0, int end = -1);
private:
	int n;
	Eigen::VectorXd* x;
	Eigen::VectorXd* b;
	std::map<long long, double>* u;
	double sigma;
	Eigen::VectorXd res;
};