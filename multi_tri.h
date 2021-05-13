#pragma once
#include"tri_class.h"
#include<mutex>
class multi_rhs :public tri_rhs
{
public:
    multi_rhs(int, double*, bool*, double*, map<int64_t, double>*, double, double*, vector<int64_t>*, mutex*, mutex*);
    ~multi_rhs();

    virtual void func(int64_t* idx, int three, int64_t index, double axbu) {

        if (axbu > 0) {

            mr->lock();
            rhs[idx[0]] += (three == 0 ? 1 : -1) * axbu * sigma;
            rhs[idx[1]] += (three == 1 ? 1 : -1) * axbu * sigma;
            rhs[idx[2]] += (three == 2 ? 1 : -1) * axbu * sigma;
            mr->unlock();
            ma->lock();
            aj1->push_back(new tri(index, idx[0], idx[1], idx[2], three));
            ma->unlock();
        }
    }
    mutex* mr;
    mutex* ma;


private:

};

multi_rhs::multi_rhs(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* rhs, vector<int64_t>* aj1, mutex* mr, mutex* ma) :
    tri_rhs(n, x, b, w, u, sigma, rhs, aj1), mr(mr), ma(ma)
{

}

multi_rhs::~multi_rhs()
{

}
class multi_calc_ap {
public:
    multi_calc_ap(vector<int64_t>* aj1,double* p)
    :aj1(aj1),p(p)
    {
    
    }
    ~multi_calc_ap() {}
    void operator()(int start, int end, Eigen::VectorXd* ap) {
        for (int i =start; i <end ; i++) {
            int64_t tem = (*aj1)[i];
            double num = -(p[tem->idx[0]] + p[tem->idx[1]] + p[tem->idx[2]]) + 2 * p[tem->idx[tem->symbol]];

            (*ap)[tem->idx[0]] -= num;
            (*ap)[tem->idx[1]] -= num;
            (*ap)[tem->idx[2]] -= num;
            (*ap)[tem->idx[tem->symbol]] += 2 * num;
        }
    }

    vector<int64_t>* aj1;
    double* p;
};


class multi_calc_au
{
public:
    multi_calc_au(int,  map<int64_t, double>*,double*);
    ~multi_calc_au();
    int n;
    map<int64_t, double>* u;
    double* au;
    double* v;
    double operator()(int init = 0, int end = -1) {
        if (end == -1)end = n;
        int var_num = n * (n - 1) / 2;
        int64_t index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
        //cout << init << "\t" << end << "\tstart\n";
        for (int i = init; i < end; i++) {
            for (int j = i + 1; j < n; ++j)
            {
                for (int k = j + 1; k < n; k++) {

                    int64_t idx[3];

                    idx[0] = ij(i, j, n);
                    idx[1] = ij(i, k, n);
                    idx[2] = ij(j, k, n);
                    for (int three = 0; three < 3; three++, index++) {
                        if (u->find(index) == u->end())
                            continue;

                        int ia[3];
                        ia[0] = ia[1] = ia[2] = -1;
                        ia[three] = 1;
                        
                        au[idx[0]] += (*u)[index] * ia[0];
                        au[idx[1]] += (*u)[index] * ia[1];
                        au[idx[2]] += (*u)[index] * ia[2];
                    }
                }

            }

        }
        double res = 0.0;
        for (int i = 0; i < var_num; i++)
            res += (au[i] + v[i]) * (au[i] + v[i]);
        return sqrt(res)/(sqrt(inner_product(v,v+var_num,v,0.0))+1);
        //cout << init <<"\t"<<end<<"\tfinished" << clock() << "\n";




    }





};

multi_calc_au::multi_calc_au(int n, map<int64_t, double>* u,double*v) :n(n), u(u),v(v)
{
    au = new double[n*(n-1)/2];
    fill(au, au + n * (n - 1) / 2, 0.0);
}

multi_calc_au::~multi_calc_au()
{
    delete[]au;
}

class multi_only_rhs :public multi_rhs
{
public:
    multi_only_rhs(int, double*, bool*, double*, map<int64_t, double>*, double, double*, vector<int64_t>*, mutex*, mutex*);
    ~multi_only_rhs();

    virtual void func(int64_t* idx, int three, int64_t index, double axbu) {

        if (axbu > 0) {

            mr->lock();
            rhs[idx[0]] += (three == 0 ? 1 : -1) * axbu * sigma;
            rhs[idx[1]] += (three == 1 ? 1 : -1) * axbu * sigma;
            rhs[idx[2]] += (three == 2 ? 1 : -1) * axbu * sigma;
            mr->unlock();
        }
    }
private:

};

multi_only_rhs::multi_only_rhs(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* rhs, vector<int64_t>* aj1, mutex* mr,mutex*ma) :
    multi_rhs(n, x, b, w, u, sigma, rhs, aj1,mr,ma)
{

}

multi_only_rhs::~multi_only_rhs()
{
    
}