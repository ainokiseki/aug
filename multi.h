#pragma once
#include"multi_tri.h"
#include<iostream>
#include<thread>
#include<mutex>
#include<vector>
#include"prox.h"
#include"math.h"
class multi
{
public:
	multi(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma, double* rhs, vector<tri*>* aj1);
	~multi();
	multi_rhs* main_tri;

	mutex* mr; 
	mutex* ma; 
	
	void operator()(int p_num) {
		int n = main_tri->n;
		vector<std::thread>threads;

		cout << "n=" << n << endl;
		for (int i = 0; i < p_num; i++) {
	
			cout << "thread\t"<<n - (int)(cbrt((i + 1)*1.0/p_num) * n ) << "\t" << n - (int)(cbrt(i *1.0 / p_num) * n) << endl;
			threads.push_back(thread(*main_tri, n - (int)(cbrt((i + 1)*1.0 / p_num) * n), n - (int)(cbrt(i *1.0/ p_num) *1.0* n)));
		}
		for (int i = 0; i < p_num; i++) {
			threads[i].join();
		}
		delete mr;
		delete ma;

	}
private:

};

multi::multi(int n, double* x, bool* b, double* w, map<long long, double>* u, double sigma, double* rhs, vector<tri*>* aj1)
{
	mr = new mutex;
	ma = new mutex;
	main_tri = new multi_rhs(n, x, b, w, u, sigma, rhs, aj1,mr,ma);
}

multi::~multi()
{
}