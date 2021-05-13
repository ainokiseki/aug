#pragma once
#include <vector>
#include<Eigen/Dense>
#include<map>
#include"tri.h"
class func_eigen
{
public:
	func_eigen();
	~func_eigen();
	virtual void operator()(double axbu,int*idx,int time, int64_t index)=0;
private:

};

class func_eigen_rhs_aj1:public func_eigen
{
public:
	func_eigen_rhs_aj1(Eigen::VectorXd* rhs,std::vector<int64_t>*aj1);
	~func_eigen_rhs_aj1();
	void operator()(double axbu, int* idx, int time,int64_t index);
private:
	Eigen::VectorXd* rhs;
	std::vector<int64_t>* aj1;
};
class func_eigen_phi:public func_eigen
{
public:
	func_eigen_phi(double,Eigen::VectorXd*,double*,int debug=0);
	~func_eigen_phi();
	void operator()(double axbu, int* idx, int time, int64_t index);
private:
	double alpha;
	Eigen::VectorXd* direction;
	double* norm;
	int debug;
};

class func_eigen_y:public func_eigen
{
public:
	func_eigen_y(std::map<int64_t,double>*,double sigma, double*,double*,double*);
	~func_eigen_y();
	void operator()(double axbu, int* idx, int time, int64_t index);
private:
	std::map<int64_t, double>* u;
	double sigma;
	double* y_norm;
	double* kkt_y;
	double* kkt_u;
};

class func_eigen_rhs:public func_eigen
{
public:
	func_eigen_rhs(double alpha, Eigen::VectorXd* direction,Eigen::VectorXd*rhs);
	~func_eigen_rhs();
	void operator()(double axbu, int* idx, int time, int64_t index);
private:
	double alpha;
	Eigen::VectorXd* direction;
	Eigen::VectorXd* rhs;
};

class func_eigen_u :public func_eigen
{
public:

	func_eigen_u(std::map<int64_t, double>* u, double sigma);
	~func_eigen_u();
	void operator()(double axbu, int* idx, int time, int64_t index);
private:
	std::map<int64_t, double>* u;
	double sigma;
};
