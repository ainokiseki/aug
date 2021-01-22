#pragma once
#include<vector>
#include<map>
#include<iostream>
#include<time.h>
using namespace std;
inline double sign(double x) {
    return (x > 0) - (x < 0);
}


long long ij(long long x, long long y, long long n) {
    long long res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
    if (res < 0)
        cout << "break\n";
    return res;
}
inline long long ij2(long long x, long long y, long long n) {
    return n * (n - 1) / 2 - n - (n - x) * (n - 1 - x) / 2 + y - 1;
}

class tri
{
public:
    long long index;
    int symbol;
    long long idx[3];
    tri(long long x, long long y, long long z, long long w, int t) {
        index = x;
        idx[0] = y;
        idx[1] = z;
        idx[2] = w;
        symbol = t;
    }
    ~tri() {}

};


template <typename T, typename A1, typename A2>
void tri_loop(int n, double* x, bool* b, map<long long, double>& u, double sigma, T func, A1& arg1 = NULL, A2& arg2 = NULL) {
    long long index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j + 1; k < n; k++) {

                int idx[3];

                idx[0] = ij(i, j, n);
                idx[1] = ij(i, k, n);
                idx[2] = ij(j, k, n);
                for (int three = 0; three < 3; three++) {
                    int ia[3];
                    ia[0] = ia[1] = ia[2] = -1;
                    ia[three] = 1;
                    double axbu = 0;
                    axbu += ia[0] * (x[idx[0]] + b[idx[0]]);
                    axbu += ia[1] * (x[idx[1]] + b[idx[1]]);
                    axbu += ia[2] * (x[idx[2]] + b[idx[2]]);

                    if (u.find(index) != u.end()) {
                        axbu += u[index] / sigma;
                    }
                    //calculation finish
                    func(idx, three, index, axbu, sigma, arg1, arg2);




                    index++;
                }
            }

        }

    }

}

void simul_loop(int init,int end,int n, double* x, bool* b, map<long long, double>& u, double sigma , double* rhs ,vector<tri*>&aj1 ) {
    long long index = n*(n-1)*(n-2)/2-(n-init)*(n-init-1)*(n-init-2)/2;
    cout <<init<<"\t" <<clock() << endl;
    for (int i = init; i < end; i++) {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j + 1; k < n; k++) {

                long long idx[3];

                idx[0] = ij(i, j, n);
                idx[1] = ij(i, k, n);
                idx[2] = ij(j, k, n);
                for (int three = 0; three < 3; three++) {
                    int ia[3];
                    ia[0] = ia[1] = ia[2] = -1;
                    ia[three] = 1;
                    double axbu = 0;
                    axbu += ia[0] * (x[idx[0]] + b[idx[0]]);
                    axbu += ia[1] * (x[idx[1]] + b[idx[1]]);
                    axbu += ia[2] * (x[idx[2]] + b[idx[2]]);

                    if (u.find(index) != u.end()) {
                        axbu += u[index] / sigma;
                    }
                    //calculation finish
                    if (axbu > 0) {
                        rhs[idx[0]] += ia[0] * axbu * sigma;
                        rhs[idx[1]] += ia[1] * axbu * sigma;
                        rhs[idx[2]] += ia[2] * axbu * sigma;
                        aj1.push_back(new tri(index, idx[0], idx[1], idx[2], three));
                    }




                    index++;
                }
            }

        }

    }
    cout <<init<< "\tfinished"<<clock()<<"\n";

}
void loop_update_lx(int* idx, int three, long long index, const double& axbu, const double& sigma, double& norm,int& aj1) {
    if (axbu > 0) {
        norm += axbu * axbu;
    }


    return;
}


void calc_ap(double* ap, vector<tri*>& aj1, double* p) {
    for (int i = 0; i < aj1.size(); i++) {
        tri* tem = aj1[i];
        double num = -(p[tem->idx[0]] + p[tem->idx[1]] + p[tem->idx[2]]) + 2 * p[tem->idx[tem->symbol]];
        ap[tem->idx[0]] -= num;
        ap[tem->idx[1]] -= num;
        ap[tem->idx[2]] -= num;
        ap[tem->idx[tem->symbol]] += 2 * num;
    }
    return;
}

double calc_phi(double* x, bool* b, double* w, double alpha, double* direction, map<long long, double>& u, double* v, int n, double sigma) {
    double norm = 0;
    int var_num = n * (n - 1) / 2;
    long long index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j + 1; k < n; k++) {

                long long idx[3];

                idx[0] = ij(i, j, n);
                idx[1] = ij(i, k, n);
                idx[2] = ij(j, k, n);
                for (int three = 0; three < 3; three++, index++) {


                    int ia, ib, ic;
                    ia = ib = ic = -1;
                    switch (three) {
                    case 0:
                        ia = 1;
                        break;
                    case 1:
                        ib = 1;
                        break;
                    case 2:
                        ic = 1;
                        break;
                    }
                    double axbu = 0;
                    axbu += ia * (x[idx[0]] + b[idx[0]] + alpha * direction[idx[0]]);
                    axbu += ib * (x[idx[1]] + b[idx[1]] + alpha * direction[idx[1]]);
                    axbu += ic * (x[idx[2]] + b[idx[2]] + alpha * direction[idx[2]]);

                    if (u.find(index) != u.end()) {
                        axbu += u[index] / sigma;
                    }
                    //calculation finish




                    if (axbu > 0) {
                        norm += axbu * axbu;
                    }

                }
            }



        }
    }
    //norm(x+v/sig-z)
    norm = norm * sigma / 2;
    //cout << "norm1:" << norm << endl;
    double nm2 = 0;
    for (int i = 0; i < var_num; i++) {
        double x_v = x[i] + alpha * direction[i] + v[i] / sigma;
        if (abs(x_v) < w[i] / sigma)
            norm += x_v * x_v * sigma / 2;

        else {
            //nm2 += w[i] * abs(x_v - sign(x_v) * w[i] / sigma);
            norm += w[i] * w[i] / (sigma * 2) + w[i] * abs(x_v - sign(x_v) * w[i] / sigma);
        }
    }
    //cout << "norm2:" << nm2 << endl;
    //cout << "final_norm:" << norm << endl;
    return norm;
}


double calc_phi(double* x, bool* b, double* w, map<long long, double>& u, double* v, int n, double sigma) {
    double norm = 0;
    int var_num = n * (n - 1) / 2;
    long long index = 0;

    tri_loop(n, x, b, u, sigma, loop_update_lx, norm, n);
    //norm(x+v/sig-z)
    norm = norm * sigma / 2;

    for (int i = 0; i < var_num; i++) {
        double x_v = x[i] + v[i] / sigma;
        if (abs(x_v) < w[i] / sigma)
            norm += x_v * x_v * sigma / 2;

        else
            norm += w[i] * w[i] / (sigma * 2) + w[i] * abs(x_v);
    }
    return norm;
}
