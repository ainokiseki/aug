#pragma once
#include<Eigen/Dense>
#include<vector>


inline int64_t ij(int64_t x, int64_t y, int64_t n);

void cacap(Eigen::VectorXd* res, int start, int end, std::vector<int64_t>* aj1, Eigen::VectorXd* rhs,int n);
