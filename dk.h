#pragma once
#include "metric.h"

class dk : public metric
{
public:
	dk(string graph_name);
	~dk();
	void triloop();
	void biloop();
	void run();
	void init();
	double gam;
	std::map<int64_t, double>correction;
	std::map<int64_t, double>next_correction;
	Eigen::MatrixXd y, m,inv_w,mat_b,p;

};

dk::dk(string graph_name) :metric(graph_name)
{
	gam = 10;
}

dk::~dk()
{
}

void dk::triloop()
{
	
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				double fm = inv_w(i, j) + inv_w(i, k) + inv_w(j, k);
				double mu;
				int64_t ijk;
				ijk = i * n * n + j * n + k;
				if (correction.find(ijk) != correction.end()) {
					y(i, j) += correction[ijk] * inv_w(i, j) / fm;
					y(i, k) -= correction[ijk] * inv_w(i, k) / fm;
					y(j, k) -= correction[ijk] * inv_w(j, k) / fm;
				}
				mu = y(i, j) + mat_b(i, j) - (y(i, k) + mat_b(i, k)) - (y(j, k) + mat_b(j, k));
				if (mu > 0) {
					y(i, j) -= mu * inv_w(i, j) / fm;
					y(i, k) += mu * inv_w(i, k) / fm;
					y(j, k) += mu * inv_w(j, k) / fm;
					next_correction.insert(std::pair<int64_t, double>(ijk,mu));
				}
				ijk = i * n * n + k * n + j;
				if (correction.find(ijk) != correction.end()) {
					y(i, k) += correction[ijk] * inv_w(i, k) / fm;
					y(i, j) -= correction[ijk] * inv_w(i, j) / fm;
					y(j, k) -= correction[ijk] * inv_w(j, k) / fm;
				}
				mu = y(i, k) + mat_b(i, k) - (y(i, j) + mat_b(i, j)) - (y(j, k) + mat_b(j, k));
				if (mu > 0) {
					y(i, k) -= mu * inv_w(i, k) / fm;
					y(i, j) += mu * inv_w(i, j) / fm;
					y(j, k) += mu * inv_w(j, k) / fm;
					next_correction.insert(std::pair<int64_t, double>(ijk, mu));
				}
				ijk = j * n * n + i * n + k;
				if (correction.find(ijk) != correction.end()) {
					y(j, k) += correction[ijk] * inv_w(j, k) / fm;
					y(i, k) -= correction[ijk] * inv_w(i, k) / fm;
					y(i, j) -= correction[ijk] * inv_w(i, j) / fm;
				}
				mu = y(j, k) + mat_b(j, k) - (y(i, k) + mat_b(i, k)) - (y(i,j) + mat_b(i,j));
				if (mu > 0) {
					y(j,k) -= mu * inv_w(j, k) / fm;
					y(i, k) += mu * inv_w(i, k) / fm;
					y(i, j) += mu * inv_w(i, j) / fm;
					next_correction.insert(std::pair<int64_t, double>(ijk, mu));
				}
			}
		}
	}

}

void dk::biloop()
{
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			double mu;
			y(i, j) -= p(i, j);
			m(i, j) -= p(i, j);
			mu = -y(i, j) - m(i, j);
			if (mu > 0) {
				y(i, j) += mu/2;
				m(i, j) += mu/2;
				p(i, j) = mu/2;
			}
			else {
				p(i, j) = 0;
			}
			y(i, j) += p(j, i);
			m(i, j) -= p(j, i);
			mu = y(i, j) - m(i, j);
			if (mu > 0) {
				y(i, j) -= mu / 2;
				m(i, j) += mu / 2;
				p(j, i) = mu / 2;
			}
			else {
				p(j, i) = 0;
			}



		}
	}
}

void dk::run()
{
	init();

	while (true) {
		triloop();
		biloop();


	}



}

void dk::init()
{
	metric::init();
	inv_w = Eigen::MatrixXd::Zero(n, n);
	mat_b = Eigen::MatrixXd::Zero(n, n);
	int index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++,index++) {
			inv_w(i, j) = 1 / w[index];
		}
		
	}
	index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++, index++) {
			mat_b(i, j) = b[index];
		}

	}
	y = Eigen::MatrixXd::Zero(n,n);
	m = -gam * Eigen::MatrixXd::Ones(n, n);
	p = Eigen::MatrixXd::Zero(n, n);
	
}
