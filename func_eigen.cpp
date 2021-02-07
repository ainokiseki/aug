
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
		aj1->push_back(new tri(index,idx[0], idx[1], idx[2],  time));
	}
}

func_eigen_phi::func_eigen_phi(double alpha, Eigen::VectorXd* direction, double* norm) :alpha(alpha), direction(direction), norm(norm)
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
}
func_eigen_u::func_eigen_u(std::map<long long, double>* u, double sigma) :u(u), sigma(sigma)
{
}

func_eigen_u::~func_eigen_u()
{
}
void func_eigen_u::operator()(double axbu, int* idx, int time, long long index)
{
	if (axbu > 0) {
		(*u)[index] = sigma * axbu;
	}
	else {

		if (u->find(index) != u->end())
			u->erase(index);
	}
}
