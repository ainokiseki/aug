#pragma once
#include<Eigen/Dense>
#include<map>
#include<unordered_map>
#include"func_eigen.h"
//#include"cg_inf.h"
class axbu_eigen
{
public:
	axbu_eigen(int n,Eigen::VectorXd*x, Eigen::VectorXd* b, std::map<int64_t,double>*u ,double sigma, std::vector<int64_t>* cons=NULL, std::unordered_map<int, int>*colidx=NULL);
	~axbu_eigen();
	void operator()(int init=0, int end=-1);


	void set_function(func_eigen* );



	int n;
	Eigen::VectorXd* x;
	Eigen::VectorXd* b;
	std::map<int64_t, double>* u;
	double sigma;
	std::unordered_map<int, int>*colidx;
	func_eigen* f;
	std::vector<int64_t>* cons;
};
class au_eigen
{
public:
	au_eigen(int, std::map<int64_t, double>*, Eigen::VectorXd ,Eigen::VectorXd, std::vector<int64_t>* cons = NULL, std::unordered_map<int, int>* colidx = NULL);
	~au_eigen();
	void operator()(Eigen::VectorXd* au, int init, int end);
	//double operator()(int init=0,int end=-1);

	std::map<int64_t, double>* u;
	int n;
	Eigen::VectorXd v;
	Eigen::VectorXd w;
	std::vector<int64_t>* cons;
	std::unordered_map<int, int>* colidx;
};


class axbu_pro
{
public:
	axbu_pro(int n, Eigen::VectorXd* x, Eigen::VectorXd* b, std::map<int64_t, double>* u, double sigma);
	~axbu_pro();
	Eigen::VectorXd operator()(int init = 0, int end = -1);
private:
	int n;
	Eigen::VectorXd* x;
	Eigen::VectorXd* b;
	std::map<int64_t, double>* u;
	double sigma;
	Eigen::VectorXd res;
};