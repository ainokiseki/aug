#pragma once
#include <vector>
#include<Eigen/Dense>
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
