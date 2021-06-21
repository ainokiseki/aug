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



func_eigen_rhs_aj1::func_eigen_rhs_aj1(Eigen::VectorXd* rhs, std::vector<int64_t>* aj1) :rhs(rhs), aj1(aj1)
{
}

func_eigen_rhs_aj1::~func_eigen_rhs_aj1()
{
}

void func_eigen_rhs_aj1::operator()(double axbu, int* idx, int time, int64_t index)
{
	if (axbu > 0) {
		
			for (int i = 0; i < 3; i++) {
				if (i == time)continue;
				(*rhs)(idx[i]) -= axbu;
			}
			(*rhs)(idx[time]) += axbu;
		
		if(aj1!=NULL)
		aj1->push_back(index);

	}
}

func_eigen_rhs_aj1_diag::func_eigen_rhs_aj1_diag(Eigen::VectorXd* rhs, std::vector<int64_t>* aj1,Eigen::VectorXd*diag, std::vector<Eigen::Triplet<double> >* triv) :rhs(rhs), aj1(aj1),diag(diag),triv(triv)
{
}

func_eigen_rhs_aj1_diag::~func_eigen_rhs_aj1_diag()
{
}

void func_eigen_rhs_aj1_diag::operator()(double axbu, int* idx, int time, int64_t index)
{
	if (axbu > 0) {

		for (int i = 0; i < 3; i++) {
			if (i == time)continue;
			(*rhs)(idx[i]) -= axbu;
		}
		(*rhs)(idx[time]) += axbu;
		(*diag)[idx[0]] += 1;
		(*diag)[idx[1]] += 1;
		(*diag)[idx[2]] += 1;

		if (time == 0) {


			triv->push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
			triv->push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
			triv->push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
		}
		else if (time == 1) {
			triv->push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
			triv->push_back(Eigen::Triplet<double>(idx[0], idx[2], 1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[2], idx[0], 1.0));
			triv->push_back(Eigen::Triplet<double>(idx[1], idx[2], -1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[2], idx[1], -1.0));
		}
		else {
			triv->push_back(Eigen::Triplet<double>(idx[0], idx[1], 1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[1], idx[0], 1.0));
			triv->push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
			triv->push_back(Eigen::Triplet<double>(idx[1], idx[2], -1.0));
			//triv->push_back(Eigen::Triplet<double>(idx[2], idx[1], -1.0));
		}
	}
}
func_eigen_phi::func_eigen_phi(double alpha, Eigen::VectorXd* direction, double* norm,int debug) :alpha(alpha), direction(direction), norm(norm),debug(debug)
{

}

func_eigen_phi::~func_eigen_phi()
{
}

void func_eigen_phi::operator()(double axbu, int* idx, int time, int64_t index)
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
func_eigen_y::func_eigen_y(std::map<int64_t, double>* u, double sigma,double*y_norm,double *kkt_y,double*kkt_u) :
	u(u), sigma(sigma),y_norm(y_norm),kkt_y(kkt_y),kkt_u(kkt_u)
{
}

func_eigen_y::~func_eigen_y()
{
}
void func_eigen_y::operator()(double axbu, int* idx, int time, int64_t index)
{

	double axb = axbu;
	if (u->find(index)!=u->end())
		axb-= (*u)[index] / sigma;

	double y = (axb < 0 ? axb : 0.0);
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

void func_eigen_rhs::operator()(double axbu, int* idx, int time, int64_t index)
{
	if (alpha > 0 && direction != NULL) {
		axbu -= alpha * (*direction)[idx[0]];
		axbu -= alpha * (*direction)[idx[1]];
		axbu -= alpha * (*direction)[idx[2]];
		axbu += 2 * alpha * (*direction)[idx[time]];
	}

	if (axbu >= 0) {
		for (int i = 0; i < 3; i++) {
			if (i == time)continue;
			(*rhs)(idx[i]) -= axbu;
		}
		(*rhs)(idx[time]) += axbu;
	}
}
func_eigen_u::func_eigen_u(std::map<int64_t, double>* u, double sigma) :
	u(u), sigma(sigma)
{
}

func_eigen_u::~func_eigen_u()
{
}
void func_eigen_u::operator()(double axbu, int* idx, int time, int64_t index)
{
	
	double y = (axbu < 0 ? axbu : 0.0);
	


	double nu = sigma * (axbu - y);

	if (nu != 0)
		(*u)[index] = nu;
	else if (u->find(index)!=u->end()){
		u->erase(index);
	}
}

func_eigen_superphi::func_eigen_superphi(double* alpha, Eigen::VectorXd* direction, double* norm, int length) :alpha(alpha), direction(direction), norm(norm), length(length)
{

}

func_eigen_superphi::~func_eigen_superphi()
{
}


void func_eigen_superphi::operator()(double axbu, int* idx, int time, int64_t index)
{
	for (int i = 0; i < length; i++) {
		double naxbu = axbu;
		naxbu -= alpha[i] * (*direction)[idx[0]];
		naxbu -= alpha[i] * (*direction)[idx[1]];
		naxbu -= alpha[i] * (*direction)[idx[2]];
		naxbu += 2 * alpha[i] * (*direction)[idx[time]];
		if (naxbu > 0) {
			norm[i] += naxbu * naxbu;
		}
	}

	//std::ofstream prt_axbu("axbu"+std::to_string(debug)+".txt", std::ios::out|std::ios::app);
	//if(alpha==0)
	//prt_axbu << index<<"\t"<<axbu << "\n";
}