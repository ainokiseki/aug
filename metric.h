#pragma once

#include<string>
#include<map>
#include<unordered_map>
#include<unordered_set>
#include <iostream>
#include <fstream>
#include<vector>
#include<set>
#include<algorithm>

#define EIGEN_USE_BLAS   
#define EIGEN_USE_MKL_VML

#include <Eigen/PardisoSupport>
#include <Eigen/IterativeLinearSolvers>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <iterator> 
#include<condition_variable>
#include<mutex>
#include<atomic>
#include<thread>
#include<mat.h>
using namespace std;

inline int ij(int x, int y, int64_t n) {
	if (x > y)return static_cast<int>(x - 1 - (static_cast<int64_t>(y) * y + 3 * y) / 2 + static_cast<int64_t>(n) * y);
	return static_cast<int>(y - 1 - (static_cast<int64_t>(x) * x + 3 * x) / 2 + static_cast<int64_t>(n) * x);
	//return n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
}
inline int ij(int x, int y, int n) {
	if (x > y)return x - 1 - (y * y + 3 * y) / 2 + n * y;
	return  y - 1 - (x * x + 3 * x) / 2 + n * x;
	//return n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
}
inline pair<int, int> rij(int ijvalue, int n) {

	int i=floor((2 * n - 1 - sqrt((1 - 2 * n) * (1 - 2 * n) - 8 * (ijvalue + 1) + 1)) / 2.0);
	int j = ijvalue - n * i + 1 + (i * i + 3 * i) / 2;
	return pair<int, int>(i, j);
}



class metric
{
public:
	metric(string graph_name, ifstream& fin);
	metric(std::string graph_name,int type=0) :sv_x(sparse_w, 0), sv_v(sparse_w, 0), sv_w(sparse_w, 0), sv_b(sparse_b) {
		ifstream graph(graph_name, ios::in);
		map<int, set<int> > g;
		int tmp_a, tmp_b;
		int input_ct = 0;
		while (graph >> tmp_a >> tmp_b) {
			g[tmp_a].insert(tmp_b);
			g[tmp_b].insert(tmp_a);
			input_ct++;
		}

		cout << "line:" << input_ct << endl;
		//point list
		vector<int> target;
		//int n= 50;
		for (auto it = g.begin(); it != g.end(); it++) {
			target.push_back(it->first);
		}
		n = target.size();
		cout << "n:" << n << endl;
		var_num = n * (n - 1) / 2;
		cons_num = n * (n - 1) * (n - 2) / 6;
		w = Eigen::VectorXd::Zero(var_num);
		b = Eigen::VectorXd::Zero(var_num);


		const double delta = 0.05;
		const double epsi = 0.01;
		//calculate w
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			int t_i = target[i];
			
			for (int j = i + 1; j < n; j++) {
				//calculate w
				int tmp_it = j - 1 - (i * i + 3 * i) / 2 + n * i;
				set<int>tmp1;
				set<int>tmp2;
				int t_j = target[j];
				set_intersection(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp1, tmp1.begin()));
				set_union(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp2, tmp2.begin()));

				double tmp = 1.0 * tmp1.size() / tmp2.size();
				if (tmp1.size() == 0) {
					b[tmp_it] = 1;
					w[tmp_it] = fabs(log((1 - delta) / (1 + delta)));
					//w[tmp_it] = 1;
					continue;
				}


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



		if (type == 1) {

			int tmp_it = 0;
			for (int i = 0; i < n; i++)
			{
				int t_i = target[i];
				for (int j = i + 1; j < n; j++, tmp_it++) {
					int t_j = target[j];
					if (g[t_i].find(t_j) != g[t_i].end()) {
						b[tmp_it] = 0;
						w[tmp_it] = 1;
					}
					else {
						b[tmp_it] = 1;
						w[tmp_it] = 1;
					}

					//calculate d
				}
			}
		}


		cout << "w:" << w.sum() << endl;
		cout << "b:" << b.sum() << endl;
		v = Eigen::VectorXd::Zero(var_num);
		x = Eigen::VectorXd::Zero(var_num);
		cout << "x,v" << x[0] << "\t" << v[0] << endl;
	}

	metric(int xnum,std::string graph_name, int type = 0) :sv_x(sparse_x, 0), sv_v(sparse_v, 0), sv_w(sparse_w, 0), sv_b(sparse_b) {
		ifstream graph(graph_name, ios::in);
		map<int, set<int> > g;
		int tmp_a, tmp_b;
		int input_ct = 0;
		while (graph >> tmp_a >> tmp_b) {
			g[tmp_a].insert(tmp_b);
			g[tmp_b].insert(tmp_a);
			input_ct++;
		}

		cout << "line:" << input_ct << endl;
		//point list
		vector<int> target;
		//int n= 50;
		for (auto it = g.begin(); it != g.end(); it++) {
			target.push_back(it->first);
		}
		n = target.size();
		cout << "n:" << n << endl;
		var_num =static_cast<int>( static_cast<int64_t>(n) * (n - 1) / 2);


		const double delta = 0.05;
		const double epsi = 0.01;
		//calculate w
		sv_w.default_value=fabs(log((1 - delta) / (1 + delta)));
		mutex mu,mu2;
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			int t_i = target[i];

			for (int j = i + 1; j < n; j++) {
				//calculate w
				int64_t tmp_it = int64_t(ij(i,j,n));
				set<int>tmp1;
				set<int>tmp2;
				int t_j = target[j];
				set_intersection(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp1, tmp1.begin()));
				set_union(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp2, tmp2.begin()));

				double tmp = 1.0 * tmp1.size() / tmp2.size();
				if (tmp1.size() == 0) {
					//w[tmp_it] = 1;
					continue;
				}


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
				if (tmp >= 0) {
					mu2.lock();
					sparse_b.insert(tmp_it);
					mu2.unlock();
				}
				mu.lock();
				sparse_w[tmp_it] = abs(tmp);
				w_value_map.push_back(pair<int, int>(i,j));
				mu.unlock();
				//calculate d
			}
		}
		double wrs=0;
		for (int i = 0; i < var_num; i++)wrs += sv_w[i];

		cout << "w:"<<wrs<<"\nb:" << sparse_b.size() <<"\nwdef:"<<sv_w.default_value << endl;

	}

	metric(MATFile* pmatFile):sv_x(sparse_w, 0), sv_v(sparse_w, 0), sv_w(sparse_w, 0),sv_b(sparse_b)
	{
		mxArray* pMxArray = matGetVariable(pmatFile, "A");

		double* initA = (double*)mxGetData(pMxArray);
		n = mxGetM(pMxArray);
		var_num = n * (n - 1) / 2;
		cons_num = n * (n - 1) * (n - 2) / 6;
		w = Eigen::VectorXd::Ones(var_num);
		b = Eigen::VectorXd::Zero(var_num);
		int ct = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				b[ct] = initA[i * n + j];
				ct++;
			}
		}

		v = Eigen::VectorXd::Zero(var_num);
		x = Eigen::VectorXd::Zero(var_num);
	}
	void savex(string s) {
		auto pmatfile = matOpen((s+string(".mat")).c_str(), "w");
		auto pmxarray = mxCreateDoubleMatrix(x.size(), 1, mxREAL);
		memcpy((void*)(mxGetPr(pmxarray)), (void*)x.data(), sizeof(double)*x.size());
		matPutVariable(pmatfile, "x", pmxarray);
		matClose(pmatfile);
		cout << "write x in " << s << ".mat" << endl;
	}
	~metric()=default;
	void init() {};
	void matrix_generate(char type);

	double kkt;
	string graph_name;
	int n, var_num;
	int64_t cons_num;
	int mthread = std::thread::hardware_concurrency();
	Eigen::VectorXd w;
	Eigen::VectorXd b;
	Eigen::VectorXd v;
	Eigen::VectorXd x;
	vector<pair<int,int>>w_value_map;
	vector<int>itervec;
	class hash_b {
	public:
		hash_b(unordered_set<int64_t> &data):data(data){}
		int default_value = 1;
		unordered_set<int64_t>&data;
		int operator[](int x) {
			auto y = data.find(x);
			if (y == data.end())return default_value;
			return 0;
		}

	};

	class hash_x {
	public:
		hash_x(unordered_map<int64_t, double>& data,double default_value) :data(data),default_value(default_value) {}
		unordered_map<int64_t, double>& data;
		double default_value;
		double operator[](int x) {
			auto y = data.find(x);
			if (y == data.end())return default_value;
			return y->second;
		}
	};

	unordered_set<int64_t> sparse_b;

	unordered_map<int64_t, double>sparse_x, sparse_v, sparse_w;
	hash_b sv_b;
	hash_x sv_x, sv_v, sv_w;
	bool cgflg = false;
	std::map<string, double>runhist;

	void SparseSave(ofstream& fout);


};
double calc_bnorm(int mthread, Eigen::VectorXd& b, int n);
