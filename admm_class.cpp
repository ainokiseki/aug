#include"admm_class.h"

void calc_au(int n, std::map<int64_t, double>& u, Eigen::VectorXd& res, int start = 0, int end = -1) {

	if (end == -1)end = n;
	for (int i = start; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;
				if (u.find(longindex[0]) != u.end()) {
					res[idx[0]] += u[longindex[0]];
					res[idx[1]] -= u[longindex[0]];
					res[idx[2]] -= u[longindex[0]];
				}
				if (u.find(longindex[1]) != u.end()) {
					res[idx[0]] -= u[longindex[1]];
					res[idx[1]] += u[longindex[1]];
					res[idx[2]] -= u[longindex[1]];
				}
				if (u.find(longindex[2]) != u.end()) {
					res[idx[0]] -= u[longindex[2]];
					res[idx[1]] -= u[longindex[2]];
					res[idx[2]] += u[longindex[2]];
				}
			}
		}
	}
}

void calc_u(int n, Eigen::VectorXd& x,Eigen::VectorXd& b, std::map<int64_t, double>& u,double sigma, int start = 0, int end = -1) {

	if (end == -1)end = n;
	for (int i = start; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);

				double axbu;
				double axbu2;
				bool flag;
				double tem;
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;

				axbu2 = axbu = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
				flag = (u.find(longindex[0]) != u.end());

				if (flag)
					axbu2 += u[longindex[0]] / sigma;

				if (axbu2 > 0) {
					u[longindex[0]] = sigma*axbu2;
				}
				else if(flag){
					u.erase(longindex[0]);
				}

				//second scene
				axbu2 = axbu = x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]];
				flag = (u.find(longindex[1]) != u.end());

				if (flag)
					axbu2 += u[longindex[1]] / sigma;

				if (axbu2 > 0) {
					u[longindex[1]] = sigma * axbu2;
				}
				else if (flag) {
					u.erase(longindex[1]);
				}
				//third scene
				axbu2 = axbu = x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]];
				flag = (u.find(longindex[2]) != u.end());

				if (flag)
					axbu2 += u[longindex[2]] / sigma;

				if (axbu2 > 0) {
					u[longindex[2]] = sigma * axbu2;
				}
				else if (flag) {
					u.erase(longindex[2]);
				}
			}
		}
	}
}

void atax(int n, Eigen::VectorXd& x, Eigen::VectorXd& res, int start = 0, int end = -1) {

	if (end == -1)end = n;
	for (int i = start; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				res[idx[0]] += 3 * x[idx[0]] - x[idx[1]] - x[idx[2]];
				res[idx[1]] += 3 * x[idx[1]] - x[idx[0]] - x[idx[2]];
				res[idx[2]] += 3 * x[idx[2]] - x[idx[1]] - x[idx[0]];
			}
		}
	}
}

void atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd&b,std::map<int64_t,double>*u,Eigen::VectorXd& res,double sigma, int start = 0, int end = -1) {

	if (end == -1)end = n;
	for (int i = start; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				
				double axbu;

				bool flag;
				double tem;
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;

				axbu = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
				flag = (u->find(longindex[0]) != u->end());
				
				if (flag)
					axbu +=(*u)[longindex[0]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu - (b[idx[0]] - b[idx[1]] - b[idx[2]]);
				if (flag) {
					tem -= (*u)[longindex[0]] / sigma;
				}

				res[idx[0]] += tem;
				res[idx[1]] -= tem;
				res[idx[2]] -= tem;
				
				//second scene
				axbu = x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]];
				flag = (u->find(longindex[1]) != u->end());

				if (flag)
					axbu += (*u)[longindex[1]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu - (b[idx[1]] - b[idx[0]] - b[idx[2]]);
				if (flag) {
					tem -= (*u)[longindex[1]] / sigma;
				}

				res[idx[1]] += tem;
				res[idx[0]] -= tem;
				res[idx[2]] -= tem;
				//third scene
				 axbu = x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]];
				flag = (u->find(longindex[2]) != u->end());

				if (flag)
					axbu += (*u)[longindex[2]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu - (b[idx[2]] - b[idx[1]] - b[idx[0]]);
				if (flag) {
					tem -= (*u)[longindex[2]] / sigma;
				}

				res[idx[2]] += tem;
				res[idx[1]] -= tem;
				res[idx[0]] -= tem;



			}
		}
	}
}
void atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd& b, std::map<int64_t, double>* u, Eigen::VectorXd& res, double sigma,int test, int start = 0, int end = -1) {

	if (end == -1)end = n;
	for (int i = start; i < end; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);

				double axbu;

				bool flag;
				double tem;
				int64_t longindex[3];
				longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
				longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
				longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;

				axbu = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
				flag = (u->find(longindex[0]) != u->end());

				if (flag)
					axbu += (*u)[longindex[0]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu ;
				if (flag) {
					tem -= (*u)[longindex[0]] / sigma;
				}

				res[idx[0]] += tem;
				res[idx[1]] -= tem;
				res[idx[2]] -= tem;

				//second scene
				axbu = x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]];
				flag = (u->find(longindex[1]) != u->end());

				if (flag)
					axbu += (*u)[longindex[1]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu ;
				if (flag) {
					tem -= (*u)[longindex[1]] / sigma;
				}

				res[idx[1]] += tem;
				res[idx[0]] -= tem;
				res[idx[2]] -= tem;
				//third scene
				axbu = x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]];
				flag = (u->find(longindex[2]) != u->end());

				if (flag)
					axbu += (*u)[longindex[2]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu;
				if (flag) {
					tem -= (*u)[longindex[2]] / sigma;
				}

				res[idx[2]] += tem;
				res[idx[1]] -= tem;
				res[idx[0]] -= tem;



			}
		}
	}
}
admm_basic::admm_basic(string graph_name):metric(graph_name)
{


	
}

admm_basic::~admm_basic()
{

}
inline double admm_basic::func_p(Eigen::VectorXd x) {
	return w.cwiseProduct(x.cwiseAbs()).sum();
}
inline Eigen::VectorXd admm_basic::prox_p(Eigen::VectorXd x, double sigma) {
	return (x.cwiseAbs() - w / sigma).cwiseMax(0.0).cwiseProduct(x.cwiseSign());
}
double admm_basic::get_kkt(double sigma) {

	return 0.0;

	//cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;
}

void admm_basic::run() {
	Eigen::VectorXd rhs(var_num);
	double sigma = 1;
	for(int k=0;k<200;k++)
	{
		rhs = Eigen::VectorXd::Zero(var_num);
		atyu(n, x, b, &u, rhs, sigma);
		rhs += prox_p(x + v / sigma, sigma) - v / sigma;
		//cg


		double tol = 1e-10;
		Eigen::VectorXd r = rhs;
		Eigen::VectorXd p(var_num);
		Eigen::VectorXd res = Eigen::VectorXd::Zero(var_num);
		double row, row_old;
		row = r.squaredNorm();
		double delta = tol * rhs.norm();
		for (int i = 0; sqrt(row) > delta; i++) {
			if (i == 0) {
				p = r;
			}
			else {
				p = r + row / row_old * p;
			}

			Eigen::VectorXd ap = Eigen::VectorXd::Zero(var_num);
			atax(n, p, ap);
			ap += p;

			double mu = row / p.dot(ap);
			res += mu * p;
			r -= mu * ap;
			row_old = row;
			row = r.squaredNorm();
		}
		//cg end

 		//cout << "res:" << (rhs - test).norm() << endl;



		calc_u(n, x, b, u, sigma);
		x = res;
		Eigen::VectorXd z = prox_p(x + v / sigma, sigma);
		v += sigma * (x - z);
		cout << func_p(x) << endl;
		Eigen::VectorXd au = Eigen::VectorXd::Zero(var_num);
		calc_au(n, u, au);
		double rd = (au + v).norm() / (1 + v.norm());
		double rp = (x - z).norm() / (1 + z.norm());
		if (rd > 10 * rp) {
			sigma = std::fmax(sigma/2,1e-5);
		}
		else if (rp > 10 * rd) {
			sigma *= 2;
		}
		//cout << "sig:"<<sigma << endl;
		//cout << rd << "\t" << rp << endl;
			
	}




}

void admm_basic::run_dense() {
	Eigen::SparseMatrix<double>A(n * (n - 1) * (n - 2)/2, var_num);

	{
		int64_t index = 0;
		std::vector<Eigen::Triplet<double> >tmp_v;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				for (int k = j + 1; k < n; k++) {
					int idx[3];
					idx[0] = ij(i, j, n);
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);

					tmp_v.push_back(Eigen::Triplet<double>(index, idx[0], 1));
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[1], -1));
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[2], -1));
					index++;
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[1], 1));
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[0], -1));
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[2], -1));
					index++;
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[2], 1));
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[1], -1));
					tmp_v.push_back(Eigen::Triplet<double>(index, idx[0], -1));
					index++;
				}
			}
		}
		A.setFromTriplets(tmp_v.begin(), tmp_v.end());
	}
	Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ll;
	{
		Eigen::VectorXd tool = Eigen::VectorXd::Ones(var_num);
		Eigen::SparseMatrix<double>ata = A.transpose() * A;
		ata+=tool.asDiagonal();
		ll.analyzePattern(ata);
		ll.factorize(ata);
	}

	Eigen::VectorXd rhs(var_num);
	double sigma = 1;
	int consnum = n * (n - 1) * (n - 2) / 2;
	Eigen::VectorXd y = Eigen::VectorXd::Zero(consnum);
	Eigen::VectorXd z = Eigen::VectorXd::Zero(var_num);
	Eigen::VectorXd vu = Eigen::VectorXd::Zero(consnum);
	Eigen::VectorXd vb = -A * b;
	for (int k = 0; k < 2000; k++)
	{
		y = A * (x+b) + vu / sigma;
		y = y.cwiseMin(0.0);


		z = prox_p(x + v / sigma, sigma);

		rhs = A.transpose() * (y  +vb-vu / sigma)+z-v/sigma;

		Eigen::VectorXd rhs2 = Eigen::VectorXd::Zero(var_num);
		atyu(n, x, b, &u, rhs2, sigma);
		cout << "pt1:" << (A.transpose() * (y + vb)).norm() << "\t" << rhs2.norm() << endl;
		rhs2 += prox_p(x + v / sigma, sigma) - v / sigma;
		//cg
		cout << rhs.norm() << "\t" << rhs2.norm() << endl;

		x = ll.solve(rhs);
		y = A * x - vb + vu / sigma;
		y = y.cwiseMin(0.0);

		vu += sigma * (A * x - y - vb);
		z = prox_p(x + v / sigma, sigma);
		v += sigma * (x - z);
		cout <<k<<"\t"<< func_p(x) << endl;
		
		double rd = (A.transpose()*vu + v).norm() / (1 + v.norm());
		double rp = (x - z).norm() / (1 + z.norm());
		if (rd > 10 * rp) {
			sigma = std::fmax(sigma / 2, 1e-5);
		}
		else if (rp > 10 * rd) {
			sigma *= 2;
		}
		//cout << "sig:" << sigma << endl;
		//cout << rd << "\t" << rp << endl;

	}




}