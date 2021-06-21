#pragma once

#include<string>
#include<map>
#include <iostream>
#include <fstream>
#include<vector>
#include<set>
#include<algorithm>
#include<Eigen/Dense>
#include <iterator> 
using namespace std;
class metric
{
public:
	metric(std::string graph_name);
	metric(int n,Eigen::VectorXd b, Eigen::VectorXd w);
	~metric();
	void init();

	string graph_name;
	int n, var_num;
	int64_t cons_num;
	Eigen::VectorXd w;
	Eigen::VectorXd b;
	Eigen::VectorXd v;
	Eigen::VectorXd x;
};
