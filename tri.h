#pragma once
#include<Eigen/Dense>
#include<vector>
class tri
{
public:
    long long index;
    int symbol;
    int idx[3];
    tri(long long x,int y,int z, int w, int t) {
        index = x;
        idx[0] = y;
        idx[1] = z;
        idx[2] = w;
        symbol = t;
    }
    ~tri() {}

};

void cacap(Eigen::VectorXd* res, int start, int end, std::vector<tri*>* aj1, Eigen::VectorXd* rhs);