#pragma once
#include<Eigen/Dense>
#include<vector>
#include<unordered_map>
#include<mutex>
#include<set>
int ij(int x, int y, int n);

void cacap(Eigen::VectorXd* res, int start, int end, std::vector<int64_t>* aj1, Eigen::VectorXd* rhs,int n);
