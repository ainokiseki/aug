#pragma once
#include "prox.h"

using namespace std;

class triple_loop
{
public:
	triple_loop(int,double*,bool*,double*, map<long long, double>* , double);
	~triple_loop();
    int n;
    double* x;
    double* w;
    bool* b;
    double sigma;
    map<long long, double>* u;
    virtual void func(long long*,int,long long,double);
	void operator()(int init=0,int end=-1) {
        if (end == -1)end = n;
       
            long long index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
            //cout << init << "\t" << end << "\tstart\n";
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

                            if (u->find(index) != u->end()) {
                                axbu += (*u)[index] / sigma;
                            }
                            //calculation finish

                            func(idx,three,index,axbu);

                            index++;
                        }
                    }

                }

            }
            //cout << init <<"\t"<<end<<"\tfinished" << clock() << "\n";

        
	
	
	}
    




};

triple_loop::triple_loop(int n,double*x,bool*b,double*w,map<long long,double>*u,double sigma):n(n),x(x),b(b),w(w),u(u),sigma(sigma)
{
}

triple_loop::~triple_loop()
{
}

void triple_loop::func(long long*idx, int three, long long index, double axbu)
{

}
class tri_rhs:public triple_loop
{
public:
    tri_rhs(int , double*, bool* , double* , map<long long, double>* , double, double*, vector<tri*>* );
    ~tri_rhs();
    double* rhs;
    vector<tri*>* aj1;
    void init_var(double* rhs, vector<tri*>* aj1) {
        rhs = rhs;
        aj1 = aj1;
    }
    virtual void func(long long* idx, int three, long long index, double axbu) {

        if (axbu > 0) {
            rhs[idx[0]] += (three==0?1:-1)* axbu * sigma;
            rhs[idx[1]] += (three == 1 ? 1 : -1) * axbu * sigma;
            rhs[idx[2]] += (three == 2 ? 1 : -1) * axbu * sigma;
            aj1->push_back(new tri(index, idx[0], idx[1], idx[2], three));
        }
    }
private:

};

tri_rhs::tri_rhs(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma,double*rhs,vector<tri*>*aj1):
    triple_loop(n,x,b,w,u,sigma),rhs(rhs),aj1(aj1)
{

}

tri_rhs::~tri_rhs()
{

}
class tri_kkt:public triple_loop
{
public:
    tri_kkt(int , double* , bool*, double*, map<long long, double>* , double,double&,double&,double&);
    ~tri_kkt();
    double& kkt_y;
    double& kkt_u;
    double& y_norm;
    virtual void func(long long* idx, int three, long long index, double axbu) {

        //kkt 0:kkt_y, 1:y_norm 2:kkt_u
        if (axbu < 0) {
            if (u->find(index) != u->end() && (*u)[index] + axbu > 0) {
                kkt_y += axbu * axbu;
            }
            else {
                if (u->find(index) != u->end())
                    kkt_y += (*u)[index] * (*u)[index];
            }
        }





        if (axbu < 0) {
            y_norm += axbu * axbu;
        }
        double tem = 0;
        tem += max(axbu, 0.0);
        if (u->find(index) != u->end())
            tem -= (*u)[index] / sigma;
        kkt_u += tem * tem;
        return;
    }
};

tri_kkt::tri_kkt(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma,double&kkt_y,double&kkt_u,double&y_norm) :
    triple_loop(n, x, b, w, u, sigma), kkt_y(kkt_y),kkt_u(kkt_u),y_norm(y_norm)
{

}

tri_kkt::~tri_kkt()
{
}


class tri_update_u:public triple_loop
{
public:
    tri_update_u(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma);
    ~tri_update_u();
    virtual void func(long long* idx, int three, long long index, double axbu) {

        if (axbu > 1e-8) {
            (*u)[index] = sigma * axbu;
        }
        else {

            if (u->find(index) != u->end())
                u->erase(index);
        }
        return;
    
    }
private:

};

tri_update_u::tri_update_u(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma) :
    triple_loop(n, x, b, w, u, sigma)
{

}
tri_update_u::~tri_update_u()
{
}


class tri_calc_phi:public triple_loop
{
public:
    tri_calc_phi(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma,double*v,double*norm,double alpha,double*direction);
    ~tri_calc_phi();
    virtual void func(long long* idx, int three, long long index, double axbu) {
        if(alpha>0&&direction!=NULL){
            axbu-=alpha * direction[idx[0]];
            axbu-=alpha * direction[idx[1]];
            axbu-=alpha * direction[idx[2]];
            axbu+=2*alpha * direction[idx[three]];
        }
        if (axbu > 0) {
            *norm += axbu * axbu;
        }
    }

    double*v;
    double alpha;
    double *direction;
    double *norm;
private:

};

tri_calc_phi::tri_calc_phi(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma,double*v,double*norm,double alpha=0,double*direction=NULL) :
    triple_loop(n, x, b, w, u, sigma),v(v),alpha(alpha),direction(direction),norm(norm)
{

}
tri_calc_phi::~tri_calc_phi()
{
}

