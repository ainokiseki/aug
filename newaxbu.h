#pragma once
#include<vector>
#include"cons_generate.h"

#include<mutex>

class almfunc
{
public:
	almfunc()=default;
	~almfunc()=default;
	virtual void operator()(double axbu, int* idx) = 0;
	virtual void integrate() {};
	std::mutex* mu;
private:

};
class almfunc_rhs :public almfunc {
public:
	almfunc_rhs(double alpha, Eigen::VectorXd* direction, Eigen::VectorXd* rhs, std::mutex* muu) :alpha(alpha), direction(direction), rhs(rhs) {
		mu = muu;
		mrhs = Eigen::VectorXd::Zero(rhs->size());
	
	};
	~almfunc_rhs()=default;
	double alpha;
	Eigen::VectorXd* direction;
	Eigen::VectorXd* rhs;
	std::mutex* mu;
	Eigen::VectorXd mrhs;
	void operator() (double axbu, int* idx){
		if (alpha > 0 && direction != NULL) {
			
			axbu += alpha * (*direction)[idx[0]]/ scale_a;
			axbu -= alpha * (*direction)[idx[1]]/ scale_a;
			axbu -= alpha * (*direction)[idx[2]]/ scale_a;
		}

		if (axbu >= 0) {

			(mrhs)(idx[1]) -= axbu/ scale_a;
			(mrhs)(idx[2]) -= axbu/ scale_a;
			
			(mrhs)(idx[0]) += axbu/ scale_a;
		}
	}
	void integrate() {
		mu->lock();
		*rhs += mrhs;
		mu->unlock();
	}

};

class almfunc_rhs_aj1_diag :public almfunc
{
public:

	almfunc_rhs_aj1_diag(Eigen::VectorXd* rhs,  Eigen::VectorXd* diag, std::vector<Eigen::Triplet<double> >* triv, std::mutex* muu) : rhs(rhs), aj1(aj1), diag(diag), triv(triv) {
		mu = muu;
		mrhs = Eigen::VectorXd::Zero(rhs->size());
		mdiag = Eigen::VectorXd::Zero(rhs->size());
	}

	~almfunc_rhs_aj1_diag()=default;
	void operator()(double axbu, int* idx) {
		if (axbu > 0) {

			mrhs(idx[1]) -= axbu/ scale_a;
			mdiag[idx[1]] += 1;
			mrhs(idx[2]) -= axbu/ scale_a;
			mdiag[idx[2]] += 1;
			

			mrhs(idx[0]) += axbu/ scale_a;
			mdiag[idx[0]] += 1;

				mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
				mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
				mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
				mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
				mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
				mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
			
			

		}
	}
	void integrate() {
		mu->lock();
		*rhs += mrhs;
		*diag += mdiag;
		triv->insert(triv->end(), mtriv.begin(), mtriv.end());
		mu->unlock();
	}
	Eigen::VectorXd* rhs;
	std::vector<int64_t>* aj1;
	Eigen::VectorXd* diag;
	std::vector<Eigen::Triplet<double> > mtriv;
	std::vector<Eigen::Triplet<double> > *triv;
	Eigen::VectorXd mrhs;

	Eigen::VectorXd mdiag;
};

class almfunc_superphi :public almfunc
{
public:
	almfunc_superphi(double* alpha, Eigen::VectorXd* direction, double* norm, int length) :alpha(alpha), direction(direction), norm(norm), length(length) {
	
		v.resize(length);
		fill(v.begin(), v.end(), 0.0);
	}
	~almfunc_superphi()=default;

	void operator()(double axbu, int* idx) {
		for (int i = 0; i < length; i++) {
			double naxbu = axbu;
			naxbu += alpha[i] * (*direction)[idx[0]]/ scale_a;
			naxbu -= alpha[i] * (*direction)[idx[1]]/ scale_a;
			naxbu -= alpha[i] * (*direction)[idx[2]]/ scale_a;
			if (naxbu > 0) {
				v[i] += naxbu * naxbu;
			}
		}
	}
	void integrate() {
		for (int i = 0; i < length; i++)norm[i] += v[i];
	}
private:
	double* alpha;
	Eigen::VectorXd* direction;
	double* norm;
	int length;
	std::vector<double>v;
};

class almfunc_u :public almfunc
{
public:

	almfunc_u(double sigma) :sigma(sigma) {}
	~almfunc_u()=default;
	void operator()(double axbu, int* idx) {
		double y = (axbu < 0 ? axbu : 0.0);


		float* u = &(((cons_full*)(void*)idx)->u);
		 *u = sigma * (axbu - y);

		

	}
private:

	double sigma;
};

class almfunc_y :public almfunc
{
public:
	almfunc_y(double sigma, double* y_norm, double* kkt_y, double* kkt_u) :sigma(sigma), y_norm(y_norm), kkt_y(kkt_y), kkt_u(kkt_u) {}

	~almfunc_y()=default;
	void operator()(double axbu, int* idx) {
		double axb = axbu;
	//	cout << "whp\n";
		double u = ((cons_full*)(void*)idx)->u;
		//axb -= u/ sigma;
		double y = (axb < 0 ? axb : 0.0);
		//double y = axbu < 0 ? axbu : 0;
		*y_norm += y * y;
	//	cout << "wtf\n";
		//*kkt_u += (axb - y) * (axb - y);
		
		//*kkt_u += u*u;
		double y_plus_u;
		y_plus_u = (y + u > 0 ? y : -u);
		*kkt_y += y_plus_u * y_plus_u;
	}

private:
	double sigma;
	double* y_norm;
	double* kkt_y;
	double* kkt_u;
};

class almfunc_uy :public almfunc
{
public:
	almfunc_uy(double sigma, double* y_norm, double* kkt_y, double* kkt_u) :sigma(sigma), y_norm(y_norm), kkt_y(kkt_y), kkt_u(kkt_u) {}

	~almfunc_uy() = default;
	void operator()(double axbu, int* idx) {

		double y = (axbu < 0 ? axbu : 0.0);


		float* u = &(((cons_full*)(void*)idx)->u);
		double new_u = sigma * (axbu - y);

		//double y = axbu < 0 ? axbu : 0;
		*y_norm += y * y;
		//	cout << "wtf\n";
			//*kkt_u += (axb - y) * (axb - y);

		*kkt_u += (axbu-*u/sigma-y)* (axbu - *u / sigma - y);
		double y_plus_u = (y + new_u > 0 ? y : -new_u);
		*kkt_y += y_plus_u * y_plus_u;

		*u = new_u;

	}

private:
	double sigma;
	double* y_norm;
	double* kkt_y;
	double* kkt_u;
};


void new_axbu(almfunc& f, concurrency::concurrent_vector<cons_full>::iterator beg, concurrency::concurrent_vector<cons_full>::iterator end, int n, Eigen::VectorXd& x, Eigen::VectorXd& b, double sigma);
void new_atu(concurrency::concurrent_vector<cons_full>::iterator beg, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& res);