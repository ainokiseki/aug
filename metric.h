#pragma once

#include<string>
#include<map>
#include<unordered_map>
#include <iostream>
#include <fstream>
#include<vector>
#include<set>
#include<algorithm>
#define EIGEN_USE_MKL_ALL   
#include<Eigen/Core>
#include<Eigen/Dense>
#include <iterator> 
#include<condition_variable>
#include<mutex>
#include<atomic>
#include<thread>
using namespace std;
class metric
{
public:
	metric(std::string graph_name,int type=0);
	~metric();
	void init();
	void matrix_generate(char type);


	string graph_name;
	int n, var_num;
	int64_t cons_num;
	int mthread = 32;
	Eigen::VectorXd w;
	Eigen::VectorXd b;
	Eigen::VectorXd v;
	Eigen::VectorXd x;
};
