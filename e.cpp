
#include <iostream>
#include <fstream>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include <iterator> 
#include<math.h>
#include<numeric>
#include<ctime>
#include <thread>
#include<string>
#include<sstream>


#include<mutex>
#include<chrono>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
#include"cg.h"
#include "tri.h"
#include"axbu_eigen.h"

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
inline long long ij(long long x, long long y, long long n) {
    long long res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
    if (res < 0)
        cout << "break\n";
    return res;
}

void loop_update_u(int* idx, int three, long long index, const double& axbu, const double& sigma, map<long long, double>& u, double& y_norm) {
    if (axbu > 1e-8) {
        u[index] = sigma * axbu;
    }
    else {

        if (u.find(index) != u.end())
            u.erase(index);
    }
    return;

}

void loop_update_kkt(int* idx, int three, long long index, const double& axbu, const double& sigma, map<long long, double>& u, double**& kkt) {
    int ia[3];
    ia[0] = ia[1] = ia[2] = -1;
    ia[three] = 1;
    //kkt 0:kkt_y, 1:y_norm 2:kkt_u
    if (axbu < 0) {
        if (u.find(index) != u.end() && u[index] + axbu > 0) {
            *kkt[0] += axbu * axbu;
        }
        else {
            if (u.find(index) != u.end())
                *kkt[0] += u[index] * u[index];
        }
    }





    if (axbu < 0) {
        *kkt[1] += axbu * axbu;
    }
    double tem = 0;
    tem += max(axbu, 0.0);
    if (u.find(index) != u.end())
        tem -= u[index] / sigma;
    *kkt[2] += tem * tem;
    return;

}




time_t start, t_end;



void ipd(double* res, double* a, double* b, int n) {
    double r = 0.0;
    for (int i = 0; i < n; i++)
        r += a[i] * b[i];
    *res = r;
    return;
}
Eigen::VectorXd prox_p_star(Eigen::VectorXd x,Eigen::VectorXd w,double sigma ) {
    
    return ;

}


int main()
{
    std::chrono::high_resolution_clock::time_point ck, ck2;

    ck = std::chrono::high_resolution_clock::now();
    string prt_name = to_string(ck.time_since_epoch().count() / 10000);
    vector<double>test_norm;
    start = clock();
    ifstream graph("graph14.txt", ios::in);
    ofstream prt_res("prt_res"+ prt_name+".txt", ios::out);
    prt_res.setf(ios_base::scientific, ios_base::floatfield);
    prt_res.precision(6);
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
    cout << "line:" << input_ct << endl;
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
    
    Eigen::VectorXd w = Eigen::VectorXd::Zero(var_num);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(var_num);
    
    Eigen::VectorXd z = Eigen::VectorXd::Zero(var_num);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(var_num);
    Eigen::VectorXd x = Eigen::VectorXd::Zero(var_num);
    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(var_num);
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

    cout << "w:" << w.sum() << endl;
    cout << "b:" << b.sum() << endl;

    //map for y
    //u
    map<long long, double>u;
    //v


    double sigma = 1;

    for (int iter_all = 0; iter_all < 1000&&sigma<8; iter_all++) {
        cout << "start main loop, sigma=" << sigma << "\n";
        prt_res << "\nmain loop,time=" << iter_all + 1 << endl;
        prt_res << "iter_time\trhs\tcg\talpha\n";
        //main loop
        axbu_eigen axbu(n, &x,&b, &u, sigma);
        for (int iter_1 = 0; iter_1 < 150; iter_1++) {//using ssn to find minimum of phi(x)
            //calculate grad,hessie
            
            cout << "start find minimum of phi(x),";
            cout << "iter time:" << iter_1 << endl;
            long long index = 0;
            prt_res << iter_1 << "\t";
            rhs = Eigen::VectorXd::Zero(var_num);
            vector<tri*>aj1;


            axbu.set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
            axbu();
            rhs = sigma * rhs + sigma*prox_p_star(x + v / sigma, w, sigma);

            cout << "aj1:" << aj1.size() << endl;




            if (rhs.norm() < 1e-6) {

                for (int i = 0; i < aj1.size(); i++)
                    delete aj1[i];
                break;
            }
            Eigen::VectorXd aj2 =(((sigma * x + v).cwiseAbs() - w).array() > 0).matrix().cast<double>();
            

            //cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;

            double a = min(pow(rhs.norm(),1.5)*0.5,0.5);

            MatrixReplacement A(&aj1,var_num,a,&aj2,sigma);
            //Eigen::BiCGSTAB<MatrixReplacement,  Eigen::IdentityPreconditioner> cg;
            Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
            cg.compute(A);


            Eigen::VectorXd direction = cg.solve(-rhs);
            
            cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
            prt_res << cg.iterations() << "\t";


            //cg end,result saved at direction
            //find step length
            double alpha = 1;
            double c1 = 1e-4;
            double lx0 = 0;
            axbu.set_function(new func_eigen_phi(alpha,&direction,&lx0));
            axbu();
            double down = lx0 + c1 * alpha*rhs.dot(direction);
            //cout << "down:" << down << "\tlx0:" << lx0 << endl;


            for (int i = 0; i < 20; i++) {
                //cout << "start find step,alpha=" << alpha << endl;
                double lx = 0;
                axbu.set_function(new func_eigen_phi(alpha, &direction, &lx));
                axbu();
                if (lx <= down) {
                    break;
                }
                alpha *= 0.5;
            }
            cout << "alpha=" << alpha << endl;
            prt_res << alpha << "\n";


            x = x + alpha * direction;
            for (int i = 0; i < aj1.size(); i++)
                delete aj1[i];
        }

        axbu.set_function(new func_eigen_u(&u, sigma));
        axbu();
        v = sigma * prox_p_star(x + v / sigma,w,sigma);

        ////update y,z
        //long long index = 0;
        //calc_z(x, z, w, v, sigma, n);
        //double y_norm, z_norm;
        //y_norm = z_norm = 0;
        //z_norm = inner_product(z, z + var_num, z, 0.0);


        ////calculate kkt
        //double kkt_y, kkt_z, kkt_u, kkt_v, kkt_x;
        //kkt_y = kkt_u = kkt_v = kkt_x = kkt_z = 0;


        //////update u
        //cout << "start u\n";
        //tri_update_u(n, x, b, w, &u, sigma)();
        ////tri_loop(n, x, b, u, sigma, loop_update_u, u,y_norm);
        //cout << "start x\n";

        //cout << "size(u):" << u.size() << endl;
        ////update v,z_norm

        //for (int i = 0; i < var_num; i++)
        //    v[i] += sigma * (x[i] - z[i]);


        //cout << "start kkt\n";
        //tri_kkt(n, x, b, w, &u, sigma, kkt_y, kkt_u, y_norm)();
        ////tri_loop(n, x, b, u, sigma, loop_update_kkt, u, kkt);
        //kkt_x = multi_calc_au(n, &u, v)();
        //kkt_y = sqrt(kkt_y) / (1 + sqrt(y_norm));
        //kkt_u = sqrt(kkt_u) / (1 + sqrt(y_norm));
        //for (int i = 0; i < var_num; i++)
        //    kkt_v += (x[i] - z[i]) * (x[i] - z[i]);
        //kkt_v = sqrt(kkt_v) / (1 + sqrt(inner_product(x, x + var_num, x, 0.0)));
        //for (int i = 0; i < var_num; i++)
        //    kkt_z += (z[i] - sign(x[i] + v[i]) * max(0.0, abs(x[i] + v[i]) - w[i])) * (z[i] - sign(x[i] + v[i]) * max(0.0, abs(x[i] + v[i]) - w[i]));
        //kkt_z = sqrt(kkt_z) / (1 + sqrt(inner_product(z, z + var_num, x, 0.0)));
        //sigma = min(sigma * 2, 1e7);
        //cout << "kkt v:" << kkt_v << "\tz:" << kkt_z << endl;

        //cout << "kkt u:" << kkt_u << "\ty:" << kkt_y << "\tx:" << kkt_x << endl;

        //prt_res << "\n" << "kkt:" << max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x)) << endl;
        //double t_res = 0.0;
        //for (int i = 0; i < var_num; i++)
        //    t_res += abs(x[i] * w[i]);
        //prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
        //if (max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x)) < 1e-5)
        //    break;

    }

    cout << "function value:" << w.dot(x.cwiseAbs())<<endl;


    return 0;
}

