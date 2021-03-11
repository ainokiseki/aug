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

class norm_basic
{
public:
	norm_basic();
	~norm_basic();
	inline Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma) {};
	void init() {
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

		int tmp_it = 0;
		const double delta = 0.05;
		const double epsi = 0.01;
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

		return;
	}
	Eigen::VectorXd get_direction(double sigma,Eigen::VectorXd rhs) {};
	double func_p(Eigen::VectorXd x){};
	double get_kkt(double sigma) {
		double y_norm = 0;
		double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
		kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
		axbu->set_function(new func_eigen_u(&u, sigma, &y_norm, &kkt_y, &kkt_u));
		axbu->operator();
		cout << "u:" << u.size() << endl;
		Eigen::VectorXd z = w.cwiseProduct(x) + v / sigma - pi_c(w.cwiseProduct(x) + v / sigma, sigma);
		v = sigma * pi_c(w.cwiseProduct(x) + v / sigma, sigma);

		y_norm = sqrt(y_norm);
		kkt_u = sqrt(kkt_u) / (1 + y_norm);
		kkt_y = sqrt(kkt_y) / (1 + y_norm);
		cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
		kkt_x = au_eigen(n, &u, v, w)();
		kkt_x = kkt_x / (1 + v.norm());

		//cout << "z_norm" << z.norm() << endl;
		kkt_v = (w.cwiseProduct(x) - z).norm() / (1 + z.norm());
		kkt_z = (v - pi_c(x + v,1)).norm() / (1 + z.norm());
		//for (int i = 0; i < var_num; i++)
		//    t_res += abs(x[i] * w[i]);
		//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
		return max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));

		//cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;




	}
	double get_step_length(Eigen::VectorXd& direction, double sigma) {

		//cg end,result saved at direction
		//find step length
		double alpha = 1;
		const double c1 = 1e-4;
		double lx0 = 0;
		axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
		axbu->operator();
		auto square = [](double x) {return x * x; };


		lx0 = sigma / 2 * lx0 + sigma / 2 * square(pi_c(w.cwiseProduct(x) + v / sigma, sigma).norm());

		lx0 += func_p(w.cwiseProduct(x) + v / sigma - pi_c(w.cwiseProduct(x) + v / sigma, sigma));

		double down = lx0 + c1 * alpha * rhs.dot(direction);
		//cout << "down:" << down << "\tlx0:" << lx0 << endl;

		cout << "lx0:" << lx0 << "\tdown" << down << endl;

		Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(var_num);
		axbu->set_function(new func_eigen_rhs(alpha, &direction, &nrhs));
		axbu->operator();
		//m_calc_rhs(8, nrhs, n, axbu, alpha, &direction);
		nrhs = sigma * nrhs;
		nrhs += sigma * w.cwiseProduct(pi_c(w.cwiseProduct(x + alpha * direction) + v / sigma, sigma));
		if (nrhs.dot(direction) < 0)
			return 1.0;

		for (int i = 0; i < 20; i++) {
			//cout << "start find step,alpha=" << alpha << endl;
			double lx = 0;
			axbu->set_function(new func_eigen_phi(alpha, &direction, &lx));
			axbu->operator();
			//m_calc_phi(8, lx, n, axbu, alpha, &direction);
			lx = sigma / 2 * lx;
			lx += sigma / 2 * square(pi_c(w.cwiseProduct(x + alpha * direction) + v / sigma, sigma).norm());
			lx += func_p(w.cwiseProduct(x + alpha * direction) + v / sigma - pi_c(w.cwiseProduct(x + alpha * direction) + v / sigma, sigma));
			if (lx <= down) {
				break;
			}
			alpha *= 0.5;
		}
		return alpha;
	
	
	};
	void run() {
		init();
		double sigma = 1;
		for (int iter_all = 0; iter_all < 1000; iter_all++, sigma *= 2) {
			axbu = new axbu_eigen(n, &x, &b, &u, sigma);
			for (int iter_1 = 0; iter_1 < 50; iter_1++) {
				rhs = Eigen::VectorXd::Zero(var_num);
				vector<tri*>aj1;
				//old version

				axbu->set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
				axbu->operator();
				rhs = sigma * rhs+sigma * w.cwiseProduct(pi_c(w.cwiseProduct(x) + v / sigma, sigma));
				if (rhs.norm() < 1e-6) {
					for (int i = 0; i < aj1.size(); i++)
						delete aj1[i];
					break;
				}
				Eigen::VectorXd direction = get_direction(sigma,rhs);
				double alpha = get_step_length(direction, sigma);
				cout << "alpha=" << alpha << endl;
				x = x + alpha * direction;
				for (int i = 0; i < aj1.size(); i++)
					delete aj1[i];
			}
			if (get_kkt(sigma) < 1e-5)
				break;



		}




	}

private:
	string graph_name;
	int n,var_num;
	long long cons_num;
	Eigen::VectorXd w;
	Eigen::VectorXd b;

	Eigen::VectorXd v;
	Eigen::VectorXd x;
	Eigen::VectorXd rhs;
	map<long long, double>u;
	axbu_eigen* axbu;
};

norm_basic::norm_basic()
{
}

norm_basic::~norm_basic()
{
}




class norm_1
{
public:
	norm_1();
	~norm_1();
	
private:

};

norm_1::norm_1()
{
}

norm_1::~norm_1()
{
}