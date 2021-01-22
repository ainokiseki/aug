#pragma once
#include"tri_class.h"
#include<mutex>
class multi_rhs :public tri_rhs
{
public:
    multi_rhs(int, double*, bool*, double*, map<long long, double>*, double, double*, vector<tri*>*,mutex*,mutex*);
    ~multi_rhs();

    virtual void func(long long* idx, int three, long long index, double axbu) {

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
    mutex*ma;


private:

};

multi_rhs::multi_rhs(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma, double* rhs, vector<tri*>* aj1,mutex*mr, mutex*ma) :
    tri_rhs(n, x, b, w, u, sigma,rhs,aj1),mr(mr),ma(ma)
{

}

multi_rhs::~multi_rhs()
{

}