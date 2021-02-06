
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
		aj1->push_back(new tri(idx[0], idx[1], idx[2], index, time));
	}
}

