#include "tri.h"
void cacap(Eigen::VectorXd* res, int start, int end, std::vector<tri*>* aj1, Eigen::VectorXd* rhs) {
    for (int i = start; i < end; i++) {
        tri* tem = (*aj1)[i];
        double num = -((*rhs)[tem->idx[0]] + (*rhs)[tem->idx[1]] + (*rhs)[tem->idx[2]]) + 2 * (*rhs)[tem->idx[tem->symbol]];

        (*res)[tem->idx[0]] -= num;

        (*res)[tem->idx[1]] -= num;

        (*res)[tem->idx[2]] -= num;

        (*res)[tem->idx[tem->symbol]] += 2 * num;
    }
}