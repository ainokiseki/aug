#include"admm_class.h"

int itertime(int x) {
	if (x > 500) {
		if (x > 3000)
			if (x > 8000)
				return 300;
			else
				return 200;
		else
			return 100;
	}
	else {
		if (x > 150)
			return 50;
		if (x > 50)
			return 25;
		if (x > 20)
			return 12;
		if (x > 10)
			return 6;
		return 3;
	}

}
void all_in_all(mutex* mu, int n, std::map<int64_t, double>& u, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, double& axby, double& y_norm, double sigma, double tau = 1.0, double scale_a = 1.0, int start = 0, int end = -1) {
	auto square = [](double x) {return x * x; };
	std::vector<std::pair<int64_t,double> >ins;
	std::vector<std::map<int64_t, double>::iterator>del;
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
				{
					double u_value = 0;
					auto t = u.find(longindex[0]);
					if (t != u.end()) {
						u_value = t->second;
					}

					double axb, axbu;


					//calc y,ax
					axbu = axb = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;


					axbu += u_value / sigma;
					y_norm += square(std::fmax(0.0, axbu));
					if (axbu < 0)
						axbu = 0;
					axby += square(axbu - u_value / sigma);
					double u_new = axbu * tau * sigma + (1 - tau) * u_value;

					if (u_new == 0) {
						if (t != u.end()) {
							del.push_back(t);
						}
					}
					else if (t!=u.end()) {
						//mu->lock();
						
						u[longindex[0]] = u_new;
						//mu->unlock();
					}
					else {
						ins.push_back(std::pair<int64_t, double>(longindex[0], u_new));
					}
					
					atu[idx[0]] += u_new / scale_a;
					atu[idx[1]] -= u_new / scale_a;
					atu[idx[2]] -= u_new / scale_a;

				}

				{
					double u_value = 0;
					auto t = u.find(longindex[1]);
					if (t != u.end()) {
						u_value = t->second;
					}

					double axb, axbu;


					//calc y,ax
					axbu = axb = (x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]]) / scale_a;


					axbu += u_value / sigma;
					y_norm += square(std::fmax(0.0, axbu));
					if (axbu < 0)
						axbu = 0;
					axby += square(axbu - u_value / sigma);
					double u_new = axbu * tau * sigma + (1 - tau) * u_value;

					if (u_new == 0) {
						if (t != u.end()) {
							del.push_back(t);
						}
					}
					else if (t != u.end()) {
						//mu->lock();

						u[longindex[1]] = u_new;
						//mu->unlock();
					}
					else {
						ins.push_back(std::pair<int64_t, double>(longindex[1], u_new));
					}
					atu[idx[1]] += u_new / scale_a;
					atu[idx[0]] -= u_new / scale_a;
					atu[idx[2]] -= u_new / scale_a;

				}

				{
					double u_value = 0;
					auto t = u.find(longindex[2]);
					if (t != u.end()) {
						u_value = t->second;
					}

					double axb, axbu;


					//calc y,ax
					axbu = axb = (x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]]) / scale_a;


					axbu += u_value / sigma;
					y_norm += square(std::fmax(0.0, axbu));
					if (axbu < 0)
						axbu = 0;
					axby += square(axbu - u_value / sigma);
					double u_new = axbu * tau * sigma + (1 - tau) * u_value;
					if (u_new == 0) {
						if (t != u.end()) {
							del.push_back(t);
						}
					}
					else if (t != u.end()) {
						//mu->lock();

						u[longindex[2]] = u_new;
						//mu->unlock();
					}
					else {
						ins.push_back(std::pair<int64_t, double>(longindex[2], u_new));
					}
					atu[idx[2]] += u_new / scale_a;
					atu[idx[1]] -= u_new / scale_a;
					atu[idx[0]] -= u_new / scale_a;

				}

			}
		}
	}
	mu->lock();
	for (auto i : ins) {
		u.insert(i);
	}
	for (auto i : del) {
		u.erase(i);
	}
	mu->unlock();
}
void m_all_in_all(int pnum,int n, std::map<int64_t, double>& u, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, double& axby, double& y_norm, double sigma, double tau = 1.0, double scale_a = 1.0) {
	std::vector<std::thread>v;
	std::vector<Eigen::VectorXd>vatu;
	std::vector<double>vaxby(pnum),vy_norm(pnum);
	std::mutex mu;
	for (int i = 0; i < pnum; i++) {
		vatu.push_back(Eigen::VectorXd::Zero(n * (n - 1) / 2));
	}
	for (int i = 0; i < pnum; i++) {

		v.push_back(std::thread(all_in_all,&mu, n, ref(u), ref(x), ref(b), ref(vatu[i]), ref(vaxby[i]), ref(vy_norm[i]), sigma, tau, scale_a, n - int((cbrt((i + 1)*1.0 / pnum) * n)), n - int((cbrt(i * 1.0 / pnum) * 1.0 * n))));
	}
	for (int i = 0; i < pnum; i++) {
		v[i].join();
		atu += vatu[i];
		axby += vaxby[i];
		y_norm += vy_norm[i];
	}



}




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
void m_atax(int pnum,int n, Eigen::VectorXd& x, Eigen::VectorXd& res) {

	std::vector<std::thread>v;
	std::vector<Eigen::VectorXd>vatu;
	for (int i = 0; i < pnum; i++) {
		vatu.push_back(Eigen::VectorXd::Zero(n * (n - 1) / 2));
	}
	for (int i = 0; i < pnum; i++) {

		v.push_back(std::thread(atax, n, ref(x),ref(vatu[i]), n - int((cbrt((i + 1) * 1.0 / pnum) * n)), n - int((cbrt(i * 1.0 / pnum) * 1.0 * n))));
	}
	for (int i = 0; i < pnum; i++) {
		v[i].join();
		res += vatu[i];
	}
}
void atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd&b,std::map<int64_t,double>*u,Eigen::VectorXd& res,double sigma,double scale_a=1, int start = 0, int end = -1) {

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

				axbu = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]])/scale_a;
				flag = (u->find(longindex[0]) != u->end());
				
				if (flag)
					axbu +=(*u)[longindex[0]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu - (b[idx[0]] - b[idx[1]] - b[idx[2]])/scale_a;
				if (flag) {
					tem -= (*u)[longindex[0]] / sigma;
				}

				res[idx[0]] += tem/scale_a;
				res[idx[1]] -= tem/scale_a;
				res[idx[2]] -= tem/scale_a;
				
				//second scene
				axbu = (x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]])/scale_a;
				flag = (u->find(longindex[1]) != u->end());

				if (flag)
					axbu += (*u)[longindex[1]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu - (b[idx[1]] - b[idx[0]] - b[idx[2]])/scale_a;
				if (flag) {
					tem -= (*u)[longindex[1]] / sigma;
				}

				res[idx[1]] += tem/scale_a;
				res[idx[0]] -= tem/scale_a;
				res[idx[2]] -= tem/scale_a;
				//third scene
				axbu = (x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]])/scale_a;
				flag = (u->find(longindex[2]) != u->end());

				if (flag)
					axbu += (*u)[longindex[2]] / sigma;

				if (axbu > 0) {
					axbu = 0;
				}
				tem = axbu - (b[idx[2]] - b[idx[1]] - b[idx[0]])/scale_a;
				if (flag) {
					tem -= (*u)[longindex[2]] / sigma;
				}

				res[idx[2]] += tem/scale_a;
				res[idx[1]] -= tem/scale_a;
				res[idx[0]] -= tem/scale_a;



			}
		}
	}
}
void m_atyu(int pnum, int n, Eigen::VectorXd& x, Eigen::VectorXd& b, std::map<int64_t, double>* u, Eigen::VectorXd& res, double sigma, double scale_a = 1) {
	std::vector<std::thread>v;
	std::vector<Eigen::VectorXd>vatu;
	for (int i = 0; i < pnum; i++) {
		vatu.push_back(Eigen::VectorXd::Zero(n * (n - 1) / 2));
	}
	for (int i = 0; i < pnum; i++) {

		v.push_back(std::thread(atyu, n, ref(x), ref(b), u, ref(vatu[i]), sigma, scale_a, n - int((cbrt((i + 1) * 1.0 / pnum) * n)), n - int((cbrt(i * 1.0 / pnum) * 1.0 * n))));
	}
	for (int i = 0; i < pnum; i++) {
		v[i].join();
		res += vatu[i];

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

std::pair<int,int> admm_basic::run(double tau,double sigma,double scale_b,double scale_c) {

	double scale_a;
	//scale_a = 3;
	//scale_b = b.norm();
	//b = b / b.norm();
	//w = w / 9;
	//scale_b = b.norm();

	b = b * 3/scale_b;

	// = b / b.norm();
	scale_a = 3;

	w = w / scale_c;
	Eigen::VectorXd rhs(var_num);

	int dual_win(0), prim_win(0);
	int maxsize = 0;
	for(int k=0;;k++)
	{
		rhs = Eigen::VectorXd::Zero(var_num);
		//atyu(n, x, b, &u, rhs, sigma,scale_a);
		m_atyu(24, n, x, b, &u, rhs, sigma, scale_a);
		rhs += prox_p(x + v / sigma, sigma) - v / sigma;
		//cg
		//cout << "atyu\n";

		double tol = 1e-10;
		Eigen::VectorXd r = rhs;
		Eigen::VectorXd p(var_num);
		Eigen::VectorXd res = Eigen::VectorXd::Zero(var_num);
		double row, row_old;
		row = r.squaredNorm();
		double delta = tol * rhs.norm();
		//cout << "cgstart\n";
		{

			int i = 0;
			for (; sqrt(row) > delta; i++) {
				if (i == 0) {
					p = r;
				}
				else {
					p = r + row / row_old * p;
				}

				Eigen::VectorXd ap = Eigen::VectorXd::Zero(var_num);
				m_atax(24,n, p, ap);
				ap /= (scale_a * scale_a);
				ap += p;

				double mu = row / p.dot(ap);
				res += mu * p;
				r -= mu * ap;
				row_old = row;
				row = r.squaredNorm();
			}
			//cout << "  i:" << i << "  ";
		}
		//cout << "cgend\n";
		//cg end

 		//cout << "res:" << (rhs - test).norm() << endl;



		x = res;
		double y_norm, kkt_u;
		y_norm = kkt_u = 0;
		Eigen::VectorXd atu = Eigen::VectorXd::Zero(var_num);
		m_all_in_all(24,n, u, x, b, atu, kkt_u, y_norm, sigma, tau, scale_a);
		Eigen::VectorXd z = prox_p(x + v / sigma, sigma);
		v += sigma *tau* (x - z);

		double kkt_x = (atu + v).norm() / (1 + v.norm());
		double true_ku = scale_b * sqrt(kkt_u) / (1 + scale_b * sqrt(y_norm));
		kkt_u = sqrt(kkt_u) / (1 + sqrt(y_norm));
		double kkt_v = (x - z).norm() / (1 + z.norm());

		//cout << "allinall\n";
		if (std::fmax(kkt_u, kkt_v) > kkt_x) {
			dual_win++;
		}
		else {
			prim_win++;
		}
		if ((k + 1) % itertime(k) == 0) {
			if (prim_win > dual_win * 1.2 && prim_win > 1) {
				prim_win = 0;
				sigma /= 3;
			}
			else if (dual_win > prim_win * 1.2 && dual_win > 1) {
				dual_win = 0;
				sigma *= 3;
			}
		}
		kkt_x = scale_a * scale_a * (atu + v).norm() / (1 + scale_a * scale_a * v.norm());

		kkt_v = scale_b / scale_a * (x - z).norm() / (1 + scale_b / scale_a * z.norm());


		double kkt = std::fmax(kkt_x, std::fmax(true_ku, kkt_v));
		maxsize = fmax(maxsize, u.size());
		cout << k << "\t" << kkt << "\t"<<func_p(x)*scale_b<<"\t"<< sigma<<"\t"<<u.size()<<endl;
		if (kkt < 1e-3)
			return std::pair<int,int>(k,maxsize);


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

	double scale_a, scale_b, scale_c;
	scale_a = 3; scale_b = vb.norm(), scale_c = 27;
	A = A / scale_a; vb = vb / scale_b; w = w / scale_c;

	for (int k = 0; k < 2000; k++)
	{
		y = A * (x+b) + vu / sigma;
		y = y.cwiseMin(0.0);


		z = prox_p(x + v / sigma, sigma);

		rhs = A.transpose() * (y  +vb-vu / sigma)+z-v/sigma;

		x = ll.solve(rhs);
		y = A * x - vb + vu / sigma;
		y = y.cwiseMin(0.0);

		vu += sigma * (A * x - y - vb);
		z = prox_p(x + v / sigma, sigma);
		v += sigma * (x - z);
		
		
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