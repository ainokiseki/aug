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
	std::pair<int, int> run(double tau = 1.618, double sigma = 1, double scale_b = 1, double scale_c = 9);


	void run_dense();

	Eigen::SparseMatrix<double>ata;

	map<int64_t, double>u;
	axbu_eigen* axbu;
};