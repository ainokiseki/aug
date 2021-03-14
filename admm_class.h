#pragma once
#include "metric.h"
#include"axbu_eigen.h"
#include"cg_admm.h"
using namespace std;

class admm_basic:public metric
{
public:
	admm_basic(string graph_name);
	~admm_basic();
	virtual inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma) = 0;
	virtual double func_p(Eigen::VectorXd x) = 0;
	double get_kkt(double sigma);
	void run();

protected:


	map<long long, double>u;
	axbu_eigen* axbu;
};

class admm_1 :public admm_basic
{
public:
	admm_1(string graph_name);
	~admm_1();
	inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);
	double func_p(Eigen::VectorXd x);
};
