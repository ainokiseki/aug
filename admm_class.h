#pragma once
#include "metric.h"
#include"axbu_eigen.h"
#include"cg_admm.h"
#include<thread>
#include<mutex>
#include<Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include<queue>
#include<string>
using namespace std;

class admm:public metric
{
public:
	admm(string graph_name, int type = 0);
	~admm();

	int itertime(int x);

	void m_all_in_all(int n, vector<vector<pair<int64_t, double>>>& u, vector<vector<pair<int64_t, double>>>& u_new, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, Eigen::VectorXd& old_x, double& axby, double& y_norm, double sigma, std::vector<int64_t>& cons, vector<int>& colidx, double tau, double scale_a);

	void m_atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd& b, vector<vector<pair<int64_t, double>>>& u, Eigen::VectorXd& res, double sigma, std::vector<int64_t>& cons, vector<int>& colidx, double scale_a);

	void constrain_generate(int istart, int iend, mutex& mut, vector<int>& colidx, int& colidx_num, std::vector<Eigen::Triplet<double>>& acons, std::vector<int64_t>& cons, std::set<int64_t>& viscons);

	//void constrain_generate(mutex& mut, vector<int>& colidx, int& colidx_num, std::vector<Eigen::Triplet<double>>& acons, std::vector<int64_t>& cons, std::set<int64_t>& viscons);

	//void constrain_generate(mutex& mut, vector<int>& colidx, int& colidx_num, std::vector<int64_t>& acons, std::vector<int64_t>& cons, std::set<int64_t>& viscons);


	//void m_atyu(int pnum, int n, Eigen::VectorXd& x, Eigen::VectorXd& b, vector<queue<pair<int64_t, double>>>& u, Eigen::VectorXd& res, double sigma, std::vector<int64_t>& cons, std::unordered_map<int, int>& colidx, double scale_a);

	//void m_all_in_all(int n, vector<vector<pair<int64_t, double>>>& u, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, Eigen::VectorXd& old_x, double& axby, double& y_norm, double sigma, std::vector<int64_t>& cons, std::unordered_map<int, int>& colidx, double tau, double scale_a);

	//void m_all_in_all(int n, vector<vector<pair<int64_t, double>>>& u, vector<vector<pair<int64_t, double>>>& u_new, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, Eigen::VectorXd& old_x, double& axby, double& y_norm, double sigma, std::vector<int64_t>& cons, std::unordered_map<int, int>& colidx, double tau, double scale_a);

	//void m_atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd& b, vector<vector<pair<int64_t, double>>>& u, Eigen::VectorXd& res, double sigma, std::vector<int64_t>& cons, std::unordered_map<int, int>& colidx, double scale_a);

};






class admm_basic :public admm
{
public:
	admm_basic(string graph_name,int type=0);
	~admm_basic();

	double func_p(Eigen::VectorXd x, Eigen::VectorXd w);

	//double func_p(Eigen::VectorXd x);

	Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma, Eigen::VectorXd w);

	



	double get_kkt(double sigma);
	double run(double tau = 1.618, double sigma = 1, double scale_b = 1, double scale_c = 1);

	std::pair<int, int> run_store(double tau, double sigma, double scale_b, double scale_c);


	void run_dense();

	Eigen::SparseMatrix<double>ata;
	Eigen::VectorXd sub_x;
	Eigen::VectorXd sub_b;
	Eigen::VectorXd sub_w;
	Eigen::VectorXd sub_v;

	axbu_eigen* axbu;
};


class admm_square :public admm
{
public:
	admm_square(string graph_name, int type = 0);
	~admm_square();

	double func_p(Eigen::VectorXd x);

	//double func_p(Eigen::VectorXd x, Eigen::VectorXd w);

	//double func_p(Eigen::VectorXd x);






	double run(double tau = 1.618, double sigma = 1, double scale_b = 1, double scale_c = 1);





	Eigen::SparseMatrix<double>ata;
	Eigen::VectorXd sub_x;
	Eigen::VectorXd sub_b;
	Eigen::VectorXd sub_w;


	axbu_eigen* axbu;
};
class admm_inf :public admm
{
public:
	admm_inf(string graph_name, int type = 0);
	~admm_inf();


	double func_p(Eigen::VectorXd x);

	//double func_p(Eigen::VectorXd x);

	Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma, Eigen::VectorXd w);





	double get_kkt(double sigma);
	Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);
	double run(double tau = 1.618, double sigma = 1, double scale_b = 1, double scale_c = 1);

	std::pair<int, int> run_store(double tau, double sigma, double scale_b, double scale_c);


	void run_dense();

	Eigen::SparseMatrix<double>ata;
	Eigen::VectorXd sub_x;
	Eigen::VectorXd sub_b;
	Eigen::VectorXd sub_w;
	Eigen::VectorXd sub_v;

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