#include<iostream>
#include<fstream>
#include<string>
#include"func_eigen.h"
func_eigen::func_eigen()
{
}

func_eigen::~func_eigen()
{
}



func_eigen_rhs_aj1::func_eigen_rhs_aj1(Eigen::VectorXd* rhs, std::vector<tri*>* aj1) :rhs(rhs), aj1(aj1)
{
}

func_eigen_rhs_aj1::~func_eigen_rhs_aj1()
{
}

void func_eigen_rhs_aj1::operator()(double axbu, int* idx, int time, long long index)
{
	if (axbu > 0) {
		for (int i = 0; i < 3; i++) {
			if (i == time)continue;
			(*rhs)(idx[i]) -= axbu;
		}
		(*rhs)(idx[time]) += axbu;
		if(aj1!=NULL)
		aj1->push_back(new tri(index,idx[0], idx[1], idx[2],  time));
	}
}

func_eigen_phi::func_eigen_phi(double alpha, Eigen::VectorXd* direction, double* norm,int debug) :alpha(alpha), direction(direction), norm(norm),debug(debug)
{

}

func_eigen_phi::~func_eigen_phi()
{
}

void func_eigen_phi::operator()(double axbu, int* idx, int time, long long index)
{
	if (alpha > 0 && direction != NULL) {
		axbu -= alpha * (*direction)[idx[0]];
		axbu -= alpha * (*direction)[idx[1]];
		axbu -= alpha * (*direction)[idx[2]];
		axbu += 2 * alpha * (*direction)[idx[time]];
	}

	if (axbu > 0) {
		*norm += axbu * axbu;
	}
	//std::ofstream prt_axbu("axbu"+std::to_string(debug)+".txt", std::ios::out|std::ios::app);
	//if(alpha==0)
	//prt_axbu << index<<"\t"<<axbu << "\n";
}
func_eigen_u::func_eigen_u(std::map<long long, double>* u, double sigma,double*y_norm,double *kkt_y,double*kkt_u) :
	u(u), sigma(sigma),y_norm(y_norm),kkt_y(kkt_y),kkt_u(kkt_u)
{
}

func_eigen_u::~func_eigen_u()
{
}
void func_eigen_u::operator()(double axbu, int* idx, int time, long long index)
{
	double axb= axbu - (*u)[index] / sigma;
	double naxbu = axb;
	if (axbu > 0) {
		naxbu   += axbu;
		(*u)[index] = sigma * axbu;
	}
	else {
		if (u->find(index) != u->end())
			u->erase(index);
	}
	double y = (naxbu < 0 ? naxbu : 0.0);
	*y_norm += y * y;
	*kkt_y += (axb - y) * (axb - y);
	double y_plus_u;
	if (u->find(index)==u->end())
	{
		y_plus_u = (y > 0 ? y : 0);
	}
	else {
		y_plus_u = (y + (*u)[index] > 0 ? y : -(*u)[index]);
	}

	*kkt_u += y_plus_u * y_plus_u;

	
}


func_eigen_rhs::func_eigen_rhs(double alpha, Eigen::VectorXd* direction, Eigen::VectorXd* rhs) :alpha(alpha), direction(direction), rhs(rhs)
{

}

func_eigen_rhs::~func_eigen_rhs()
{
}

void func_eigen_rhs::operator()(double axbu, int* idx, int time, long long index)
{
	if (alpha > 0 && direction != NULL) {
		axbu -= alpha * (*direction)[idx[0]];
		axbu -= alpha * (*direction)[idx[1]];
		axbu -= alpha * (*direction)[idx[2]];
		axbu += 2 * alpha * (*direction)[idx[time]];
	}

	if (axbu > 0) {
		for (int i = 0; i < 3; i++) {
			if (i == time)continue;
			(*rhs)(idx[i]) -= axbu;
		}
		(*rhs)(idx[time]) += axbu;
	}
}
