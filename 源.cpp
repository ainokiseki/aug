#include<iostream>
#include<Eigen/Dense>
using namespace std;
int main() {
	Eigen::VectorXd a=Eigen::VectorXd::Zero(50);
	

	for (int i = 0; i < 50; i++)a(i)=i ;
	for (int i = 0; i < 50; i++)cout<<a(i)<<endl;

	return 0;


}