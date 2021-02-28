
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
#include<omp.h>

using namespace std;
typedef map<int, set<int> > mis;

map<int, set<int> > g;




double epsi = 0.01;


inline long long ij(long long x, long long y, long long n) {
    long long res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
    if (res < 0)
        cout << "break\n";
    return res;
}






time_t start, t_end;



 inline Eigen::VectorXd prox_p_star(Eigen::VectorXd x,double sigma ) {
    Eigen::VectorXd res;

    res = x.cwiseAbs().cwiseMin(1/ sigma);
    res = res.cwiseProduct(x.cwiseSign());
    res = res.cast<double>();
    return res;

}
inline Eigen::VectorXd prox_p(Eigen::VectorXd x, Eigen::VectorXd w, double sigma) {
    Eigen::VectorXd res, res2;
    res = (x.cwiseAbs() - w / sigma).cwiseMax(0);
    res = res.cwiseProduct(x.cwiseSign());
    res = res.cast<double>();
    return res;



}
void m_calc_rhs_and_aj1(int pnum,Eigen::VectorXd&rhs,vector<tri*>&aj1,int n,axbu_eigen axbu) {
    int var_num = n * (n - 1) / 2;
    Eigen::VectorXd* mrhs=new Eigen::VectorXd[pnum];
    vector<tri*>*maj1=new vector<tri*>[pnum];
    for (int i = 0; i < pnum; i++)
        mrhs[i] = Eigen::VectorXd::Zero(var_num);
    vector<thread>tv;
    for (int i = 0; i < pnum; i++) {
        axbu.set_function(new func_eigen_rhs_aj1(mrhs+i, maj1+i));
        tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
    }
    for (int i = 0; i < pnum; i++)
        tv[i].join();
    for (int i=0;i<pnum;i++)
        rhs += mrhs[i];
    for (int i=0;i<pnum;i++)
        aj1.insert(aj1.end(), maj1[i].begin(), maj1[i].end());
    delete[]mrhs;
    delete[]maj1;
}
void m_calc_rhs(int pnum, Eigen::VectorXd& rhs,int n, axbu_eigen axbu,double alpha,Eigen::VectorXd*direction) {
    int var_num = n * (n - 1) / 2;
    Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];
    for (int i = 0; i < pnum; i++)
        mrhs[i] = Eigen::VectorXd::Zero(var_num);
    vector<thread>tv;
    for (int i = 0; i < pnum; i++) {
        axbu.set_function(new func_eigen_rhs(alpha,direction,mrhs + i));
        tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
    }
    for (int i = 0; i < pnum; i++)
        tv[i].join();
    for (int i = 0; i < pnum; i++)
        rhs += mrhs[i];
    delete[]mrhs;
}
void m_calc_phi(int pnum, double& phi, int n, axbu_eigen axbu, double alpha, Eigen::VectorXd* direction) {
    int var_num = n * (n - 1) / 2;
    double*mphi=new double[pnum];
    for (int i = 0; i < pnum; i++)
        mphi[i] = 0;
    vector<thread>tv;
    for (int i = 0; i < pnum; i++) {
        axbu.set_function(new func_eigen_phi(alpha, direction,mphi+i,i));
        int a, b;
        a = n - (int)(cbrt((i + 1) * 1.0 / pnum) * n);
        b = n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n);
        tv.push_back(thread(axbu, a,b));

    }
    for (int i = 0; i < pnum; i++)
        tv[i].join();
    for (int i = 0; i < pnum; i++)
        phi += mphi[i];
    
    delete[]mphi;
}
int main()
{
    omp_set_num_threads(8);
    std::chrono::high_resolution_clock::time_point ck, ck2;

    ck = std::chrono::high_resolution_clock::now();
    string prt_name = to_string(ck.time_since_epoch().count() / 10000);
    vector<double>test_norm;
    start = clock();
    ifstream graph("graph500.txt", ios::in);
    ofstream prt_res("prt_res"+ prt_name+".txt", ios::out);

    prt_res.setf(ios_base::scientific, ios_base::floatfield);
    prt_res.precision(6);
    ofstream prt_x("prt_x" + prt_name + ".txt", ios::out);
    prt_res.setf(ios_base::scientific, ios_base::floatfield);
    prt_res.precision(10);

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

    for (int iter_all = 0; iter_all < 1000; iter_all++) {
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
            

            rhs = Eigen::VectorXd::Zero(var_num);
            vector<tri*>aj1;
            //old version

            axbu.set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
            axbu();
            //old version end
            //m_calc_rhs_and_aj1(8, rhs, aj1, n, axbu);
            //Eigen::VectorXd mrhs[8];
            //vector<tri*>maj1[8];
            //for (int i = 0; i < 8; i++)
            //    mrhs[i] = Eigen::VectorXd::Zero(var_num);

            //vector<thread>tv;
            //for (int i = 0; i < 8; i++) {
            //    axbu.set_function(new func_eigen_rhs_aj1(&mrhs[i], &maj1[i]));
            //    tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / 8) * n), n - (int)(cbrt(i * 1.0 / 8) * 1.0 * n)));
            //}
            //for (int i = 0; i < 8; i++)
            //    tv[i].join();
            //for (auto i : mrhs)
            //    rhs += i;
            //for (auto i : maj1)
            //    aj1.insert(aj1.end(), i.begin(), i.end());
            //
            rhs = sigma * rhs;



            rhs+=sigma*w.cwiseProduct( prox_p_star(x + v / sigma,  sigma));


            
            cout << "aj1:" << aj1.size() << endl;

            cout << "rhs:" << rhs.norm() << endl;

            if (rhs.norm() < 1e-6) {

                for (int i = 0; i < aj1.size(); i++)
                    delete aj1[i];
                break;
            }
            Eigen::VectorXd aj2 =((w.cwiseProduct(x)+v/sigma).cwiseAbs().array() < 1/sigma).matrix().cast<double>();
            aj2 = aj2.cwiseProduct(w.cwiseProduct(w));
            
            //cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;

            double a = min(pow(rhs.norm(),1.5)*0.5,0.5);

            MatrixReplacement A(&aj1,var_num,a,&aj2,sigma);
            //Eigen::BiCGSTAB<MatrixReplacement,  Eigen::IdentityPreconditioner> cg;
            Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
            cg.compute(A);


            Eigen::VectorXd direction = cg.solve(-rhs);
            
            //cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
            cout << "CG:" << cg.iterations() << std::endl;





            //cg end,result saved at direction
            //find step length
            double alpha = 1;
            double c1 = 1e-4;
            double lx0 = 0;
            //old version


            //axbu.set_function(new func_eigen_phi(0,&direction,&lx0));
            //axbu();

            m_calc_phi(8, lx0, n, axbu, 0, &direction);
            auto square = [](double x) {return x * x; };
            lx0 = sigma / 2 * lx0 + sigma/2*square((w.cwiseProduct(x)+v/sigma).cwiseMin(1/sigma).norm());
            lx0 += w.dot( prox_p(x + v / sigma, w, sigma).cwiseAbs());
            double down = lx0 + c1 * alpha*rhs.dot(direction);
            //cout << "down:" << down << "\tlx0:" << lx0 << endl;

            cout <<"lx0:"<<lx0<<"\tdown"<< down << endl;
            Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(var_num);
            axbu.set_function(new func_eigen_rhs(alpha, &direction, &nrhs));
            axbu();
            //m_calc_rhs(8, nrhs, n, axbu, alpha, &direction);
            nrhs = sigma * nrhs;
            nrhs += sigma * prox_p_star(x+alpha*direction + v / sigma, sigma);
            if (nrhs.dot(direction) < 0)
                goto AL;

            for (int i = 0; i < 20; i++) {
                //cout << "start find step,alpha=" << alpha << endl;
                double lx = 0;
                axbu.set_function(new func_eigen_phi(alpha, &direction, &lx));
                axbu();
                //m_calc_phi(8, lx, n, axbu, alpha, &direction);
                lx = sigma / 2 * lx;
                lx+=sigma / 2 * prox_p_star(x + alpha * direction + v / sigma, w, sigma).norm() * prox_p_star(x + alpha * direction + v / sigma, w, sigma).norm();
                lx += w.dot(prox_p(x + v / sigma + alpha * direction, w, sigma).cwiseAbs());

                if (lx <= down) {
                    break;
                }
                alpha *= 0.5;
            }
        AL:
            cout << "alpha=" << alpha << endl;
            prt_res << iter_1<<" "<<lx0<<" "<<rhs.norm()<<"\t"<<cg.iterations()<<"\t"<<alpha << "\n";
            x = x + alpha * direction;
            for (int i = 0; i < aj1.size(); i++)
                delete aj1[i];
        }
        double y_norm = 0;
        double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
        kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
        axbu.set_function(new func_eigen_u(&u, sigma,&y_norm,&kkt_y,&kkt_u));
        axbu();
        cout << "u:" << u.size() << endl;
        z = prox_p(x + v / sigma, w, sigma);
        v = sigma * prox_p_star(x + v / sigma,w,sigma);
        sigma *= 2;
        y_norm = sqrt(y_norm);
        kkt_u = sqrt(kkt_u) / (1 + y_norm);
        kkt_y = sqrt(kkt_y) / (1 + y_norm);
        //cout << "kkt_u:" << kkt_u  << "\tkkt_y:" <<kkt_y << endl;
        kkt_x = au_eigen(n, &u, v)();
        kkt_x = kkt_x / (1 + v.norm());

        //cout << "z_norm" << z.norm() << endl;
        kkt_v = (x - z).norm() / (1 + x.norm());
        kkt_z = (z - prox_p(x + v,w,1)).norm() / (1 + z.norm());
        //for (int i = 0; i < var_num; i++)
        //    t_res += abs(x[i] * w[i]);
        //prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
        double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
        //cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
        cout << "kkt:" << kkt_final << endl;
        //cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;
        if (kkt_final < 1e-5)
            break;

    }

    cout << "function value:" << w.dot(x.cwiseAbs())<<endl;


    return 0;
}

