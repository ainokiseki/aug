#include "tri.h"

inline int64_t ij(int64_t x, int64_t y, int64_t n) {
    if (x > y)std::swap(x, y);
	int64_t res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;

	return res;
}


void cacap(Eigen::VectorXd* res, int start, int end, std::vector<int64_t>* aj1, Eigen::VectorXd* rhs,int n) {
    for (int i = start; i < end; i++) {
        int64_t tem = (*aj1)[i];


        int ci, cj, ck;
        ck = tem % n;
        cj = tem / n;
        ci = cj / n;
        cj = cj % n;
        int idx[3];
        idx[0] = ij(ci, cj, n);
        idx[1] = ij(ci, ck, n);
        idx[2] = ij(cj, ck, n);


        double num = (*rhs)[idx[0]] - (*rhs)[idx[1]] - (*rhs)[idx[2]];

        (*res)[idx[0]] += num;

        (*res)[idx[1]] -= num;

        (*res)[idx[2]] -= num;
    }
}