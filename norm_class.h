


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
#include"cgdiag.h"
#include"axbu_eigen.h"
//#include<omp.h>

class norm_ori :public metric
{
public:
	norm_ori(string graph_name);
	norm_ori(string graph_name, string data);
	~norm_ori();

	inline Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma) ;
	inline Eigen::VectorXd prox_p_star (Eigen::VectorXd x, double sigma);
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);

	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1, Eigen::VectorXd*diag);

	Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag);

	inline double func_p(Eigen::VectorXd x);
	double get_kkt(double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma);
	double get_step_length(Eigen::VectorXd& direction, double sigma, int numer);
	//Eigen::VectorXd get_direction_hand(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1, Eigen::VectorXd* diag);
	void run();
	double time_record = 0;
	double cst = 0;
	int iter_all = 0;
	double sigma = 1;
	int mthread = 8;
	ofstream prt_res;
	Eigen::VectorXd rhs,x0;
	map<int64_t, double>u;
	axbu_eigen* axbu;
	double beta;
};
