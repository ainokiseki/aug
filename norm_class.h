


#include"metric.h"
#include <Eigen/PardisoSupport>

#include<thread>
#include<math.h>
#include<numeric>
#include<string>
#include<sstream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>

#include"cg.h"
#include"cg2.h"
#include"cg_inf.h"
#include"cg3.h"
#include"cgdiag.h"
#include"axbu_eigen.h"
//#include<omp.h>


class norm_basic:public metric
{
public:
	norm_basic(string graph_name);
	~norm_basic();
	virtual inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma)=0;
	virtual Eigen::VectorXd get_direction(double sigma,Eigen::VectorXd rhs,vector<int64_t>*aj1)=0;
	virtual double func_p(Eigen::VectorXd x)=0;
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);;
	void run();




	Eigen::VectorXd rhs;
	map<int64_t, double>u;
	axbu_eigen* axbu;
};
class norm_ori :public metric
{
public:
	norm_ori(string graph_name);
	norm_ori(string graph_name, string data);
	norm_ori(int n,Eigen::VectorXd b, Eigen::VectorXd w);

	~norm_ori();

	inline Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma) ;
	inline Eigen::VectorXd prox_p_star (Eigen::VectorXd x, double sigma);
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);

	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1, Eigen::VectorXd*diag);

	Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag);

	Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag, int ld);

	void run(int looplim=100,double sig=1,double kkt_limit=1e-6);

	inline double func_p(Eigen::VectorXd x);
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma, int numer);
	Eigen::VectorXd get_direction_hand(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1, Eigen::VectorXd* diag);

	double time_record = 0;
	double cst = 0;
	int iter_all = 0;
	double sigma = 1;
	int mthread = 8;
	ofstream prt_res;
	Eigen::VectorXd rhs,x0;
	map<int64_t, double>u;
	axbu_eigen* axbu;
};

class norm_1:public norm_basic
{
public:
	norm_1(string graph_name);
	~norm_1();
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);
	double func_p(Eigen::VectorXd x);
	inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);
private:

};

class norm_2 :public norm_basic
{
public:
	norm_2(string graph_name);
	~norm_2();
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);
	double func_p(Eigen::VectorXd x);
	inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);
	void init();
private:

};

class norm_inf :public norm_basic
{
public:
	norm_inf(string graph_name);
	~norm_inf();
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);
	double func_p(Eigen::VectorXd x);
	inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);
private:

};


class norm_2_if
{
public:
	norm_2_if(string graph_name);
	~norm_2_if();
	virtual void init();
	virtual Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);
	virtual double func_p(Eigen::VectorXd x);
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);;
	void run();

protected:
	string graph_name;
	int n, var_num;
	int64_t cons_num;
	Eigen::VectorXd w;
	Eigen::VectorXd b;
	Eigen::VectorXd x;
	Eigen::VectorXd rhs;
	map<int64_t, double>u;
	axbu_eigen* axbu;
};
