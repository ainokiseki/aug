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
	~metric();
	void init();
protected:
	string graph_name;
	int n, var_num;
	long long cons_num;
	Eigen::VectorXd w;
	Eigen::VectorXd b;
	Eigen::VectorXd v;
	Eigen::VectorXd x;
};
