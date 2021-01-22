
#include <iostream>
#include <fstream>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include <iterator> 
#include<math.h>
#include<numeric>
#include "prox.h"
using namespace std;
typedef map<int, set<int> > mis;

map<int, set<int> > g;





double gamma = 2;
double epsi = 0.01;

//inline double sign(double x) {
//    return (x > 0) - (x < 0);
//}
//
//
//inline long long ij(long long x, long long y, long long n) {
//    return n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
//}
//inline long long ij2(long long x, long long y, long long n) {
//    return n * (n - 1) / 2 - n - (n - x) * (n - 1 - x) / 2 + y - 1;
//}

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

void loop_update_u(int* idx, int three, long long index, const double& axbu, const double& sigma, map<long long, double>& u, double& y_norm) {
    if (axbu > 1e-8) {
        u[index] = sigma * axbu;
    }
    else {
        y_norm += axbu * axbu;
        if (u.find(index) != u.end())
            u.erase(index);
    }
    return;

}

void loop_update_rhs(int* idx, int three, long long index, const double& axbu, const double& sigma, double*& rhs, vector<tri*>& aj1) {
    int ia[3];
    ia[0] = ia[1] = ia[2] = -1;
    ia[three] = 1;
    if (axbu > 0) {
        rhs[idx[0]] += ia[0] * axbu * sigma;
        rhs[idx[1]] += ia[1] * axbu * sigma;
        rhs[idx[2]] += ia[2] * axbu * sigma;
        aj1.push_back(new tri(index, idx[0], idx[1], idx[2], three));
    }
    return;

}
void calc_rhs(double* x, bool* b, double* w, map<long long, double>& u, double* v, double* rhs, int n, double sigma) {
    int var_num = n * (n - 1) / 2;
    long long index = 0;
    std::fill(rhs, rhs + var_num, 0);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j + 1; k < n; k++) {

                int idx[3];

                idx[0] = ij(i, j, n);
                idx[1] = ij(i, k, n);
                idx[2] = ij(j, k, n);
                for (int three = 0; three < 3; three++) {


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
                    axbu += ia * (x[idx[0]] + b[idx[0]]);
                    axbu += ib * (x[idx[1]] + b[idx[1]]);
                    axbu += ic * (x[idx[2]] + b[idx[2]]);

                    if (u.find(index) != u.end()) {
                        axbu += u[index] / sigma;
                    }
                    //calculation finish




                    if (axbu > 0) {
                        rhs[idx[0]] += ia * axbu * sigma;
                        rhs[idx[1]] += ib * axbu * sigma;
                        rhs[idx[2]] += ic * axbu * sigma;

                    }
                    index++;
                }
            }



        }
    }
    //calculate rhs
    for (int i = 0; i < var_num; i++) {
        if (abs(x[i] * sigma + v[i]) < w[i]) {
            rhs[i] += x[i] * sigma + v[i];

        }
        else {
            rhs[i] += sign(x[i] * sigma + v[i]) * w[i];

        }

    }

}



void calc_z(double* x, double* z, double* w, double* v, double sigma, int n) {
    int var_num = n * (n - 1) / 2;
    for (int i = 0; i < var_num; i++) {
        double x_v = x[i] + v[i] / sigma;
        if (abs(x_v) > w[i] / sigma) {
            z[i] = x_v - sign(x_v) * w[i] / sigma;
        }
        else {
            z[i] = 0;
        }
    }

    return;
}


int main()
{


    ifstream graph("graph14.txt", ios::in);
    ofstream prt_res("prt_res.txt", ios::out);
    int tmp_a, tmp_b;
    //build graph
    int input_ct = 0;
    while (graph >> tmp_a >> tmp_b) {
        g[tmp_a].insert(tmp_b);
        g[tmp_b].insert(tmp_a);
        input_ct++;
        if (input_ct > 300)
            break;
    }
    //point list
    vector<int> target;
    //int n= 50;
    for (mis::iterator it = g.begin(); it != g.end(); it++) {
        target.push_back(it->first);
        //if (target.size() >= n)
        //   break;
    }
    int n = target.size();
    cout << "n:" << n << endl;
    int var_num = n * (n - 1) / 2;
    long long cons_num = n * (n - 1) * (n - 2) / 6;
    double* w = new double[var_num];
    std::fill(w, w + var_num, 0);

    bool* b = new bool[var_num];
    double* z = new double[var_num];

    int tmp_it = 0;
    double delta = 0.05;
    //calculate w
    for (int i = 0; i < n; i++)
    {
        int t_i = target[i];
        for (int j = i + 1; j < n; j++, tmp_it++) {
            //calculate w
            set<int>tmp1;
            set<int>tmp2;

            int t_j = target[j];
            set_intersection(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp1, tmp1.begin()));
            set_union(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp2, tmp2.begin()));
            double tmp = 1.0 * tmp1.size() / tmp2.size();
            tmp = log((1 + tmp - delta) / (1 - tmp + delta));
            if (tmp > 0) {
                tmp += epsi;
            }
            else if (tmp < 0) {
                tmp -= epsi;
            }
            else {
                if (g[t_i].find(t_j) != g[t_i].end())
                    tmp = epsi;
                else
                    tmp = -epsi;
            }
            b[tmp_it] = (tmp < 0);
            w[tmp_it] = abs(tmp);
            //calculate d


        }
    }

    cout << "w:" << accumulate(w, w + var_num, 0.0) << endl;
    cout << "b:" << accumulate(b, b + var_num, 0.0) << endl;
    for (int i = 0; i < var_num; i++)
        prt_res << w[i] << "\t" << b[i] << endl;
    //map for y
    //u
    map<long long, double>u;
    //v
    double* v = new double[var_num];
    std::fill(v, v + var_num, 0);
    double* x = new double[var_num];
    std::fill(x, x + var_num, 0);

    double sigma = 1;
    while (sigma < 5) {
        cout << "start main loop, sigma=" << sigma << "\n";
        //main loop
        for (int iter_1 = 0; iter_1 < 30; iter_1++) {//using ssn to find minimum of phi(x)
            //calculate grad,hessie
            cout << "start find minimum of phi(x),";
            cout << "iter time:" << iter_1 << endl;
            long long index = 0;
            double* rhs = new double[var_num];
            std::fill(rhs, rhs + var_num, 0);
            vector<tri*>aj1;
            double* aj2 = new double[var_num];

            tri_loop(n, x, b, u, sigma, loop_update_rhs, rhs, aj1);

            //for (int i = 0; i < n; i++) {
            //    for (int j = i + 1; j < n; ++j)
            //    {
            //        for (int k = j + 1; k < n; k++) {

            //            int idx[3];

            //            idx[0] = ij(i , j , n);
            //            idx[1] = ij(i , k , n);
            //            idx[2] = ij(j, k, n);
            //            for (int three = 0; three < 3; three++) {


            //                int ia, ib, ic;
            //                ia = ib = ic = -1;
            //                switch (three) {
            //                case 0:
            //                    ia = 1;
            //                    break;
            //                case 1:
            //                    ib = 1;
            //                    break;
            //                case 2:
            //                    ic = 1;
            //                    break;
            //                }
            //                double axbu = 0;
            //                axbu += ia * (x[idx[0]] + b[idx[0]]);
            //                axbu += ib * (x[idx[1]] + b[idx[1]]);
            //                axbu += ic * (x[idx[2]] + b[idx[2]]);

            //                if (u.find(index) != u.end()) {
            //                    axbu += u[index] / sigma;
            //                }
            //                //calculation finish




            //                if (axbu > 0) {
            //                    rhs[idx[0]] += ia * axbu*sigma;
            //                    rhs[idx[1]] += ib * axbu*sigma;
            //                    rhs[idx[2]] += ic * axbu*sigma;
            //                    aj1.push_back(new tri(index, idx[0], idx[1], idx[2], three));
            //                }
            //                index++;
            //            }
            //        }



            //    }
            //}
            //calculate rhs
            for (int i = 0; i < var_num; i++) {
                if (abs(x[i] * sigma + v[i]) < w[i]) {
                    rhs[i] += x[i] * sigma + v[i];
                    aj2[i] = 1;
                }
                else {
                    rhs[i] += sign(x[i] * sigma + v[i]) * w[i];
                    aj2[i] = 0;
                }

            }
            //cout << "rhs:" << accumulate(rhs, rhs + var_num, 0.0) << endl;
            //cg
            double* ap = new double[var_num];
            double* direction = new double[var_num];
            std::fill(direction, direction + var_num, 0);
            double* p = new double[var_num];
            //std::fill(p,p+var_num,0);
            double* r = new double[var_num];
            copy(rhs, rhs + var_num, r);
            double r_old, r_new;
            r_old = 0;
            r_new = inner_product(r, r + var_num, r, 0.0);
            double r_norm = inner_product(rhs, rhs + var_num, rhs, 0.0);
            //cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;
            double a = min(r_norm, 20.0);
            for (int i = 0; i < 50; i++) {
                if (i == 0) {
                    copy(r, r + var_num, p);
                }
                else {
                    double beta = r_new / r_old;
                    for (int j = 0; j < var_num; j++) {
                        p[j] = beta * p[j] + r[j];
                    }
                }
                std::fill(ap, ap + var_num, 0);
                //debug
                //for (int j = 0; j < var_num; j++)
                //    p[j] = j;
                //calculate A*p
                calc_ap(ap, aj1, p);
                for (int j = 0; j < var_num; j++)
                    ap[j] = ap[j] * sigma + sigma * aj2[j] * p[j] + a * p[j];
                //calculate finish!

                //for (int j = 0; j <var_num; j++)
                //    cout << ap[j]<< endl;
                //cout << "sigma:" << sigma << endl;
                //return 0;

                double alpha = r_new / inner_product(p, p + var_num, ap, 0.0);
                //cout << "alpha:" << alpha << endl;


                for (int j = 0; j < var_num; j++) {
                    direction[j] += alpha * p[j];
                    r[j] -= alpha * ap[j];
                }
                r_old = r_new;
                r_new = inner_product(r, r + var_num, r, 0.0);
                //cout << "cg num:" << inner_product(p, p + var_num, ap, 0.0) << endl;
                if (inner_product(p, p + var_num, ap, 0.0) < 1e-30) {
                    cout << "cg finished,times=" << i << endl;
                    //for (int j = 0; j < var_num; j++)
                    //    cout << direction[j] << endl;
                    //return 0;
                    break;
                }

            }

            for (int i = 0; i < var_num; i++)
                direction[i] = -direction[i];

            delete[]ap;

            delete[]p;

            delete[]r;
            //cg end,result saved at direction
            //find step length
            double alpha = 1;
            double c1 = 1e-4;

            double lx0 = calc_phi(x, b, w, u, v, n, sigma);
            double down = lx0 + c1 * alpha * inner_product(rhs, rhs + var_num, direction, 0.0);
            //cout << "down:" << down << "\tlx0:" << lx0 << endl;

            for (int i = 0; i < 10; i++) {
                //cout << "start find step,alpha=" << alpha << endl;
                double norm = calc_phi(x, b, w, alpha, direction, u, v, n, sigma);
                //cout << "LX=" << norm << endl;
                //norm(Ax-b+sigma*u)
                if (norm <= down)
                    break;
                //p(z)
                alpha *= 0.5;
            }


            for (int i = 0; i < var_num; i++)
                x[i] += alpha * direction[i];
            calc_rhs(x, b, w, u, v, rhs, n, sigma);
            if (inner_product(rhs, rhs + var_num, rhs, 0.0) < 1e-10)
                break;
            delete[]direction;
            delete[]aj2;

        }
        //update y,z
        long long index = 0;
        calc_z(x, z, w, v, sigma, n);
        double y_norm, z_norm;
        y_norm = z_norm = 0;
        z_norm = inner_product(z, z + var_num, z, 0.0);

        double kkt_y, kkt_z, kkt_u, kkt_v, kkt_x;
        kkt_y = kkt_u = kkt_v = kkt_x = kkt_z = 0;


        ////update u
        tri_loop(n, x, b, u, sigma, loop_update_u, u, y_norm);
        //for (int i = 0; i < n; i++) {
        //    for (int j = i + 1; j < n; ++j)
        //    {
        //        for (int k = j + 1; k < n; k++) {

        //            int idx[3];

        //            idx[0] = ij(i , j , n);
        //            idx[1] = ij(i , k, n);
        //            idx[2] = ij(j , k, n);
        //            for (int three = 0; three < 3; three++) {
        //                int ia, ib, ic;
        //                ia = ib = ic = -1;
        //                switch (three) {
        //                case 0:
        //                    ia = 1;
        //                    break;
        //                case 1:
        //                    ib = 1;
        //                    break;
        //                case 2:
        //                    ic = 1;
        //                    break;
        //                }
        //                double axbu = 0;
        //                axbu += ia * (x[idx[0]] + b[idx[0]]);
        //                axbu += ib * (x[idx[1]] + b[idx[1]]);
        //                axbu += ic * (x[idx[2]] + b[idx[2]]);
        //                double ax_b = axbu;
        //                if (u.find(index) != u.end()) {
        //                    axbu += u[index] / sigma;
        //                }
        //                //calculation finish

        //                if (axbu > 0) {
        //                    u[index] = sigma*axbu;

        //                }
        //                else {
        //                    y_norm += axbu * axbu;

        //                }
        //                index++;
        //            }
        //        }

        //    }

        //}
        cout << "size(u):" << u.size() << endl;
        //update v,z_norm


        //for (int i = 0; i < var_num; i++) {
        //    if (abs(x[i] * sigma + v[i]) > w[i]) {
        //        test_v[i] = sign(x[i] * sigma + v[i])*w[i];
        //        
        //    }
        //    else {
        //        test_v[i] = x[i]*sigma + v[i];
        //    }
        //}
        for (int i = 0; i < var_num; i++)
            v[i] += sigma * (x[i] - z[i]);





        sigma = sigma * 2;
        //calculate kkt
        cout << "kkt u:" << kkt_u << "\tv:" << kkt_v << "\ty:" << kkt_y << endl;

    }
    double res_f = 0;
    for (int i = 0; i < var_num; i++)
        res_f += abs(x[i] * w[i]);
    cout << "function value:" << res_f;

    delete[]x;
    delete[]w;
    delete[]v;
    delete[]b;
    delete[]z;





    return 0;
}
