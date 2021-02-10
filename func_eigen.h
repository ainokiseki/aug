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
	virtual void operator()(double axbu,int*idx,int time,long long index)=0;
private:

};

class func_eigen_rhs_aj1:public func_eigen
{
public:
	func_eigen_rhs_aj1(Eigen::VectorXd* rhs,std::vector<tri*>*aj1);
	~func_eigen_rhs_aj1();
	void operator()(double axbu, int* idx, int time,long long index);
private:
	Eigen::VectorXd* rhs;
	std::vector<tri*>* aj1;
};
class func_eigen_phi:public func_eigen
{
public:
	func_eigen_phi(double,Eigen::VectorXd*,double*);
	~func_eigen_phi();
	void operator()(double axbu, int* idx, int time, long long index);
private:
	double alpha;
	Eigen::VectorXd* direction;
	double* norm;
};

class func_eigen_u:public func_eigen
{
public:
	func_eigen_u(std::map<long long,double>*,double sigma, double*,double*,double*);
	~func_eigen_u();
	void operator()(double axbu, int* idx, int time, long long index);
private:
	std::map<long long, double>* u;
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
	void operator()(double axbu, int* idx, int time, long long index);
private:
	double alpha;
	Eigen::VectorXd* direction;
	Eigen::VectorXd* rhs;
};

