


#include"metric.h"

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
#include"axbu_eigen.h"
//#include<omp.h>
#include<mat.h>

class norm_basic:public metric
{
public:
	norm_basic(string graph_name);
	~norm_basic();
	virtual inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma)=0;
	virtual Eigen::VectorXd get_direction(double sigma,Eigen::VectorXd rhs,vector<tri*>*aj1)=0;
	virtual double func_p(Eigen::VectorXd x)=0;
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);;
	void run();




	Eigen::VectorXd rhs;
	map<long long, double>u;
	axbu_eigen* axbu;
};
class norm_ori :public metric
{
public:
	norm_ori(string graph_name);
	~norm_ori();

	inline Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma) ;
	inline Eigen::VectorXd prox_p_star (Eigen::VectorXd x, double sigma);
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<tri*>* aj1);
	inline double func_p(Eigen::VectorXd x);
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);
	void run();

	double cst = 0;


	ofstream prt_res;
	Eigen::VectorXd rhs,x0;
	map<long long, double>u;
	axbu_eigen* axbu;
};

class norm_1:public norm_basic
{
public:
	norm_1(string graph_name);
	~norm_1();
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<tri*>* aj1);
	double func_p(Eigen::VectorXd x);
	inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma);
private:

};

class norm_2 :public norm_basic
{
public:
	norm_2(string graph_name);
	~norm_2();
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<tri*>* aj1);
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
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<tri*>* aj1);
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
	virtual Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<tri*>* aj1);
	virtual double func_p(Eigen::VectorXd x);
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);;
	void run();

protected:
	string graph_name;
	int n, var_num;
	long long cons_num;
	Eigen::VectorXd w;
	Eigen::VectorXd b;
	Eigen::VectorXd x;
	Eigen::VectorXd rhs;
	map<long long, double>u;
	axbu_eigen* axbu;
};
