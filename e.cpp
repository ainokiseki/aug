#include"norm_class.h"



int main() {

	norm_inf nm("123");

	Eigen::VectorXd x(2);
	x << 1, -1.5;
	cout<<nm.pi_c(x,1);


    return 0;
}
