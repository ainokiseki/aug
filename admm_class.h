#pragma once
#include "metric.h"
#include"axbu_eigen.h"
#include"cg_admm.h"
#include<thread>
#include<mutex>
#include<Eigen/Sparse>
#include<Eigen/SparseCholesky>
using namespace std;

class admm_basic :public metric
{
public:
	admm_basic(string graph_name);
	~admm_basic();

	double func_p(Eigen::VectorXd x);

	Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma);

	double get_kkt(double sigma);
	double run(double tau = 1.618, double sigma = 1, double scale_b = 1, double scale_c = 1);

	std::pair<int, int> run_store(double tau, double sigma, double scale_b, double scale_c);


	void run_dense();

	Eigen::SparseMatrix<double>ata;

	map<int64_t, double>u;
	axbu_eigen* axbu;
};

class admm_base :public metric
{
public:
	admm_base(string graph_name);
	~admm_base();

	virtual double func_p(Eigen::VectorXd x);

	virtual Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);

	double get_kkt(double sigma);
	double run(double tau = 1.618, double sigma = 1, double scale_b = 1, double scale_c = 1);


	void run_dense();

	Eigen::SparseMatrix<double>ata;

	map<int64_t, double>u;
	axbu_eigen* axbu;
};
class admm_2 :public admm_base
{
public:
	admm_2(string graph_name);
	~admm_2();

	virtual double func_p(Eigen::VectorXd x);

	virtual Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);

};