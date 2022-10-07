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
	multi(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* rhs, vector<int64_t>* aj1);
	~multi();
	multi_rhs* main_tri;

	mutex* mr;
	mutex* ma;

	void operator()(int p_num) {
		int n = main_tri->n;
		vector<std::thread>threads;

		//cout << "n=" << n << endl;
		for (int i = 0; i < p_num; i++) {

			//cout << "thread\t"<<n - (int)(cbrt((i + 1)*1.0/p_num) * n ) << "\t" << n - (int)(cbrt(i *1.0 / p_num) * n) << endl;
			threads.push_back(thread(*main_tri, n - (int)(cbrt((i + 1) * 1.0 / p_num) * n), n - (int)(cbrt(i * 1.0 / p_num) * 1.0 * n)));
		}
		for (int i = 0; i < p_num; i++) {
			threads[i].join();
		}


	}
private:

};

multi::multi(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* rhs, vector<int64_t>* aj1)
{
	mr = new mutex;
	ma = new mutex;
	main_tri = new multi_rhs(n, x, b, w, u, sigma, rhs, aj1, mr,ma);
}

multi::~multi()
{
	delete mr;
	delete ma;
	delete main_tri;
}

class multi_calc_phi
{
public:
	multi_calc_phi(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* v, double alpha=0, double* direction=NULL) :
		sigma(sigma), v(v), n(n), x(x), b(b), w(w), u(u), alpha(alpha), direction(direction)
	{
	}
	~multi_calc_phi() {

	}
	int n;
	double* x;
	double* w;
	bool* b;
	double sigma;
	map<int64_t, double>* u;
	double* v;
	double alpha;
	double* direction;
	double operator()(int p_num) {

		vector<std::thread>threads;
		vector<tri_calc_phi*>v2;
		//cout << "n=" << n << endl;
		double* norm = new double[10];
		for (int i = 0; i < p_num; i++)
			norm[i] = 0.0;
		for (int i = 0; i < p_num; i++)
			v2.push_back(new tri_calc_phi(n, x, b, w, u, sigma, v, norm + i, alpha, direction));
		for (int i = 0; i < p_num; i++) {
			//cout << "thread\t"<<n - (int)(cbrt((i + 1)*1.0/p_num) * n ) << "\t" << n - (int)(cbrt(i *1.0 / p_num) * n) << endl;
			threads.push_back(
				thread(
					*v2[i], n - (int)(cbrt((i + 1) * 1.0 / p_num) * n),
					n - (int)(cbrt(i * 1.0 / p_num) * 1.0 * n)
				)
			)
				;
		}
		for (int i = 0; i < p_num; i++) {
			threads[i].join();
		}
		for (int i = 0; i < p_num; ++i)
		{
			delete v2[i];
		}


		double r1 = 0.0;
		for (int i = 0; i < p_num; i++)
			r1 += norm[i];
		r1 = r1 * sigma / 2;
		//cout << "norm1:" << norm << endl;
		int var_num = n * (n - 1) / 2;
		for (int i = 0; i < var_num; i++) {

			double x_v = x[i] + v[i] / sigma;
			if (alpha > 0 && direction != NULL)
				x_v += +alpha * direction[i];
			if (abs(x_v) < w[i] / sigma)
				r1 += x_v * x_v * sigma / 2;

			else {
				//nm2 += w[i] * abs(x_v - sign(x_v) * w[i] / sigma);
				r1 += w[i] * w[i] / (sigma * 2) + w[i] * abs(x_v - sign(x_v) * w[i] / sigma);
			}
		}
		//cout << "norm2:" << nm2 << endl;
		//cout << "final_norm:" << norm << endl;
		delete[] norm;
		return r1;
	}


};

class multi_once
{
public:
	multi_once(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* rhs, vector<int64_t>* aj1);
	~multi_once();
	multi_rhs* main_tri;

	mutex* mr;
	mutex* ma;

	void operator()(int p_num) {
		int n = main_tri->n;
		vector<std::thread>threads;

		//cout << "n=" << n << endl;
		for (int i = 0; i < p_num; i++) {

			//cout << "thread\t"<<n - (int)(cbrt((i + 1)*1.0/p_num) * n ) << "\t" << n - (int)(cbrt(i *1.0 / p_num) * n) << endl;
			threads.push_back(thread(*main_tri, n - (int)(cbrt((i + 1) * 1.0 / p_num) * n), n - (int)(cbrt(i * 1.0 / p_num) * 1.0 * n)));
		}
		for (int i = 0; i < p_num; i++) {
			threads[i].join();
		}


	}
private:

};

multi_once::multi_once(int n, double* x, bool* b, double* w, map<int64_t, double>* u, double sigma, double* rhs, vector<int64_t>* aj1)
{
	mr = new mutex;
	ma = new mutex;
	main_tri = new multi_only_rhs(n, x, b, w, u, sigma, rhs, aj1, mr, ma);
}

multi_once::~multi_once()
{
	delete mr;
	delete ma;
	delete main_tri;
}
