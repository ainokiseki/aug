
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
#include<ctime>
#include <thread>
#include<future>
#include"tri_class.h"
#include"multi.h"
#include<mutex>
#include<chrono>
#include<Eigen/Dense>
using namespace std;
typedef map<int, set<int> > mis;

map<int, set<int> > g;




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


void loop_update_u(int*idx,int three,long long index,const double &axbu,const double &sigma, map<long long, double>&u,double&y_norm) {
    if (axbu > 1e-8) {
        u[index] = sigma * axbu;
    }
    else {

        if (u.find(index) != u.end())
            u.erase(index);
    }
    return;

}

void loop_update_kkt(int*idx,int three, long long index, const double& axbu, const double& sigma, map<long long, double>& u, double**& kkt) {
    int ia[3];
    ia[0] = ia[1] = ia[2] = -1;
    ia[three] = 1;
    //kkt 0:kkt_y, 1:y_norm 2:kkt_u
    if (axbu < 0) {
        if (u.find(index) != u.end() && u[index] + axbu > 0) {
            *kkt[0] += axbu * axbu;
        }
        else {
            if(u.find(index) != u.end())
            *kkt[0] += u[index] * u[index];
        }
    } 
        


        
    
    if (axbu < 0) {
        *kkt[1] += axbu * axbu;
    }
    double tem=0;
    tem += max(axbu, 0.0);
    if (u.find(index) != u.end())
        tem -= u[index] / sigma;
    *kkt[2] += tem * tem;
    return;

}

void loop_update_rhs(int* idx, int three, long long index, double axbu, double sigma, double* rhs, vector<tri*> aj1) {
    int ia[3];
    ia[0] = ia[1] = ia[2] = -1;
    ia[three] = 1;
    if (axbu > 0) {
        rhs[idx[0]] += ia[0] * axbu * sigma;
        rhs[idx[1]] += ia[1] * axbu * sigma;
        rhs[idx[2]] += ia[2] * axbu * sigma;
        //aj1.push_back(new tri(index, idx[0], idx[1], idx[2], three));
    }
    return;

}
void loop_update_aj(int* idx, int three, long long index, const double& axbu, const double& sigma, double*& rhs, vector<tri*>& aj1) {

    if (axbu > 0) {

        aj1.push_back(new tri(index, idx[0], idx[1], idx[2], three));
    }
    return;
}
void te(int n) {
    cout << n << endl;
}
void calc_z(double* x, double* z,double* w,double* v,double sigma, int n) {
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
time_t start, t_end;



void ipd(double*res,double*a,double*b,int n) {
    double r = 0.0;
    for (int i = 0; i < n; i++)
        r += a[i] * b[i];
    *res = r;
    return;
}

double multi_inner_product(double* a, double* b, int n) {
    vector<std::thread>v;
    double* res = new double[8];
    for (int i = 0; i < 8; i++)
        v.push_back(std::thread(ipd,res+i, a+i*n/8,b+i*n/8,n/8));
    for (int i=0;i<8;i++)
        v[i].join();
    return accumulate(res, res + 8,0.0);
}


int main()
{
    std::chrono::high_resolution_clock::time_point ck,ck2;


    vector<double>test_norm;
    start = clock();
    ifstream graph("graph14.txt", ios::in);
    ofstream prt_res("prt_res.txt", ios::out);
    prt_res.setf(ios_base::scientific, ios_base::floatfield);
    prt_res.precision(2);
    int tmp_a, tmp_b;
    //build graph
    int input_ct = 0;
    while (graph >> tmp_a >> tmp_b) {
        g[tmp_a].insert(tmp_b);
        g[tmp_b].insert(tmp_a);
        input_ct++;
        //if (input_ct > 300)
        //    break;
    }
    cout << "line:"<<input_ct << endl;
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
        for (int j = i + 1; j < n; j++,tmp_it++) {
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
            b[tmp_it] = (tmp<0);
            w[tmp_it] =  abs(tmp);
            //calculate d
            

        }
    }

    cout << "w:" << accumulate(w, w + var_num, 0.0) << endl;
    cout << "b:" << accumulate(b, b + var_num, 0.0)<<endl;

    //map for y
    //u
    map<long long, double>u;
    //v
    double* v = new double[var_num];
    std::fill(v, v + var_num, 0);
    double* x = new double[var_num];
    std::fill(x, x + var_num, 0);
    double* rhs = new double[var_num];
    double sigma = 1;

    for (int iter_all = 0; iter_all < 1000;iter_all++) {
        cout << "start main loop, sigma="<<sigma<<"\n";
        prt_res << "\nmain loop,time=" << iter_all+1 << endl;
        prt_res << "iter_time\trhs\tcg\talpha\n" ;
        //main loop
        for (int iter_1 = 0; iter_1 < 150;iter_1++) {//using ssn to find minimum of phi(x)
            //calculate grad,hessie
            cout << "start find minimum of phi(x)," ;
            cout << "iter time:" << iter_1 << endl;
            long long index = 0;
            prt_res << iter_1<<"\t";
            std::fill(rhs, rhs + var_num, 0);
            vector<tri*>aj1;
            double* aj2 = new double[var_num];
            double* rhs_t[4];
            //for (auto& i : rhs_t) {
            //    i = new double[var_num];
            //    std::fill(i,i + var_num, 0);
            //}

            //std::thread t[3];
            //for(int i=4;i>1;i--)
            //    t[4-i]= std::thread(simul_loop, n-(int)(sqrt(i)*n/2), n - (int)(sqrt(i-1) * n/2), n, x, b, u, sigma,  rhs_t[4-i]);
            //std::thread t = std::thread(simul_loop, 0, 275, n, x, b, u, sigma, rhs_t[0]);
            
            //tri_rhs(n, x, b, w, &u, sigma, rhs, &aj1)();
            //multi_rhs ttt(n, x, b, w, &u, sigma, rhs, &aj1);
            multi tttt(n, x, b, w, &u, sigma, rhs, &aj1);
            tttt(8);

            //simul_loop(0, n, n, x, b, u, sigma, rhs,aj1);
            
            
            //t.join();
            //for (int i = 0; i < 3; i++)
            //    t[i].join();
            //for (int i = 0; i < var_num; i++) {

            //    for (auto j : rhs_t) {
            //        rhs[i] += j[i];
            //    }
            //}
            //tri_loop(n,  x,b,  u,  sigma,loop_update_rhs, rhs,aj1);


            //tri_loop(n, x, b, u, sigma, loop_update_aj, rhs, aj1);
            //tic();
            cout << "aj1:" << aj1.size() << endl;

            //calculate rhs
            for (int i = 0; i < var_num; i++) {
                if (abs(x[i] * sigma + v[i]) < w[i]) {
                    rhs[i] += x[i] * sigma + v[i];
                    aj2[i] = 1;
                }
                else {
                    rhs[i] += sign(x[i] * sigma + v[i])*w[i];
                    aj2[i] = 0;
                }

            }
            prt_res << inner_product(rhs, rhs + var_num, rhs, 0.0) << "\t";
            if (inner_product(rhs, rhs + var_num, rhs, 0.0) < 1e-12) {
                delete[]aj2;
                for (int i = 0; i < aj1.size(); i++)
                    delete aj1[i];

                break;
            }
            

            Eigen::VectorXd vdirection=Eigen::VectorXd::Zero(var_num);

            Eigen::VectorXd p(var_num);
            //std::fill(p,p+var_num,0);
            Eigen::VectorXd r(var_num);
            for (int i = 0; i < var_num; i++)
                r[i]= rhs[i];
            for(int i=0;i<var_num;i++)
               p[i]=rhs[i];

            double r_old, r_new;
            r_old = 0;
            r_new = r.dot(r);
            double r_norm = r_new;

            //cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;
            double a = min(r_norm, 20.0);
            int maxit = max(50, (int)sqrt(cons_num));

            for (int i = 0; i < maxit; i++) {

                Eigen::VectorXd ap=Eigen::VectorXd::Zero(var_num);
                //debug
                //for (int j = 0; j < var_num; j++)
                //    p[j] = j;
                //calculate A*p

                calc_ap(ap, aj1, p);


                for (int j = 0; j < var_num; j++)
                    ap[j] = ap[j] * sigma + sigma * aj2[j] * p[j] + a * p[j];

                double alpha = r_new / p.dot(ap);
                //$double alpha = r_new / multi_inner_product(p,ap,var_num);
                //cout << "alpha:" << alpha << endl;

                vdirection += alpha * p;
                r -= alpha * ap;

                r_old = r_new;

                r_new = r.dot(r);
                //r_new = multi_inner_product(r,r,var_num);

                //cout << "cg num:" << inner_product(p, p + var_num, ap, 0.0) << endl;
                if (p.norm() < 1e-8) {
                    cout << "cg finished,times=" << i << endl;
                    prt_res << i << "\t";
                    //for (int j = 0; j < var_num; j++)
                    //    cout << direction[j] << endl;
                    //return 0;
                    break;
                }

                double beta = r_new / r_old;
                p = beta * p + r;
                //cout << i<<"\t"<<alpha << "\t" << beta << endl;

            }

            double* direction = new double[var_num];

            for (int i = 0; i < var_num; i++)
                direction[i] = vdirection[i];

            for (int i = 0; i < var_num; i++)
                direction[i] = -direction[i];


            //cg end,result saved at direction
            //find step length
            double alpha = 1;
            double c1 = 1e-4;

            //double lx0 = calc_phi(x, b, w, u, v, n, sigma);
            double lx0 = multi_calc_phi(n, x, b, w, &u, sigma, v)(8);

            double down=lx0+c1*alpha* inner_product(rhs, rhs + var_num,direction, 0.0);
            cout << "down:" << down << "\tlx0:" << lx0 << endl;
            

            for (int i = 0; i < 20;i++) {
                //cout << "start find step,alpha=" << alpha << endl;
                double norm=multi_calc_phi(n, x, b, w, &u, sigma, v,alpha,direction)(8);
                //double norm = calc_phi(x,b,w,alpha,direction,u,v,n,sigma);

                //cout << "LX=" << norm << endl;
                //norm(Ax-b+sigma*u)
                if (norm <= down) {
                    cout << "alpha=" << alpha << endl;
                    break;
                }
                //p(z)
                alpha *= 0.5;
            }
            prt_res << alpha << "\n";
            cout<<"time of phi"<<endl;


            for (int i = 0; i < var_num; i++)
                x[i] += alpha * direction[i];
            delete[]direction;
            delete[]aj2;
            for (int i = 0; i <aj1.size(); i++)
                delete aj1[i];



        }
        //update y,z
        long long index=0;
        calc_z(x, z, w, v, sigma, n);
        double y_norm, z_norm;
        y_norm = z_norm = 0;
        z_norm = inner_product(z, z + var_num, z, 0.0);


        //calculate kkt
        double kkt_y,kkt_z,kkt_u,kkt_v,kkt_x;
        kkt_y = kkt_u = kkt_v = kkt_x = kkt_z = 0;


        ////update u
        cout << "start u\n";
        tri_update_u(n, x, b, w, &u, sigma)();
        //tri_loop(n, x, b, u, sigma, loop_update_u, u,y_norm);
        cout << "start x\n";

        cout << "size(u):" << u.size() << endl;
        //update v,z_norm

        for (int i = 0; i < var_num; i++)
            v[i] += sigma * (x[i] - z[i]);


        cout << "start kkt\n";
        tri_kkt(n, x, b, w, &u, sigma, kkt_y, kkt_u, y_norm)();
        //tri_loop(n, x, b, u, sigma, loop_update_kkt, u, kkt);
        kkt_x = multi_calc_au(n, &u, v)();
        kkt_y = sqrt(kkt_y) / (1 + sqrt(y_norm));
        kkt_u = sqrt(kkt_u) / (1 + sqrt(y_norm));
        for (int i = 0; i < var_num; i++)
            kkt_v += (x[i] - z[i]) * (x[i] - z[i]);
        kkt_v = sqrt(kkt_v) / (1+sqrt(inner_product(x,x + var_num, x, 0.0)));
        for (int i = 0; i < var_num; i++)
            kkt_z += (z[i] - sign(x[i]+v[i])*max(0.0,abs(x[i]+v[i])-w[i])) * (z[i] - sign(x[i] + v[i]) * max(0.0, abs(x[i] + v[i]) - w[i]));
        kkt_z = sqrt(kkt_z) / (1 + sqrt(inner_product(z, z + var_num, x, 0.0)));
        sigma = min(sigma *2,1e7);
        cout << "kkt v:" << kkt_v << "\tz:" << kkt_z << endl;

        cout << "kkt u:" << kkt_u  << "\ty:" << kkt_y <<"\tx:"<<kkt_x<<endl;

        prt_res << "\n" << "kkt:" << max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x)) << endl;
        double t_res = 0.0;
        for (int i = 0; i < var_num; i++)
            t_res += abs(x[i] * w[i]);
        prt_res << "f(Ax-b)+p(x) \t" <<t_res << endl;
        if (max(max(kkt_u, kkt_y),max(max(kkt_v,kkt_z),kkt_x)) < 1e-5)
            break;

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
    delete[]rhs;


    return 0;
}

