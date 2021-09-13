#include"admm_class.h"
#include"mytime.h"
#include <Eigen/PardisoSupport>
admm::admm(string graph_name, int type) :metric(graph_name, type)
{
}

admm::~admm()
{
}
int admm::itertime(int x) {
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
void all_in_all(mutex* mu, int n, vector<pair<int64_t, double>>&u , vector<pair<int64_t, double>>& unew, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, Eigen::VectorXd& old_x, double& axby, double& y_norm, double sigma, std::vector<int64_t>& cons, vector<int>&colidx, double tau = 1.0, double scale_a = 1.0, int start = 0, int end = -1) {
	auto square = [](double x) {return x * x; };
	std::vector<std::pair<int64_t, double> >ins;
	std::vector<std::map<int64_t, double>::iterator>del;
	if (end == -1)end = n;
	int uidx = 0;
	for (int istart = start; istart < end; istart++) {

		int64_t lindex = cons[istart];
		int i, j, k;
		k = lindex % n;
		j = lindex / n;
		i = j / n;
		j = j % n;

		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];


		int64_t longindex[3];
		longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
		longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
		longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;


		{
			double u_value = 0;
			
			if (uidx<u.size()&&longindex[0]==u[uidx].first) {
				u_value = u[uidx].second;
				uidx++;
			}

			double axb, axbu;

			axbu = (old_x[idx[0]] + b[idx[0]] - old_x[idx[1]] - b[idx[1]] - old_x[idx[2]] - b[idx[2]]) / scale_a;
			//calc y,ax
			axb = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;


			axbu += u_value / sigma;
			if (axbu > 0)
				axbu = 0;
			y_norm += square(axbu);
			axb -= axbu;
			axby += square(axb);
			double u_new = axb * tau * sigma + u_value;

			if (u_new!= 0) {
				
				unew.push_back(std::pair<int64_t, double>(longindex[0], u_new));
			}

			atu[idx[0]] += u_new / scale_a;
			atu[idx[1]] -= u_new / scale_a;
			atu[idx[2]] -= u_new / scale_a;

		}

		{
			double u_value = 0;
			if (uidx < u.size() && longindex[1] == u[uidx].first) {
				u_value = u[uidx].second;
				uidx++;
			}

			double axb, axbu;


			//calc y,ax

			axbu = (old_x[idx[1]] + b[idx[1]] - old_x[idx[0]] - b[idx[0]] - old_x[idx[2]] - b[idx[2]]) / scale_a;
			//calc y,ax
			axb = (x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]]) / scale_a;


			axbu += u_value / sigma;
			if (axbu > 0)
				axbu = 0;
			y_norm += square(axbu);
			axb -= axbu;
			axby += square(axb);
			double u_new = axb * tau * sigma + u_value;


			if (u_new != 0) {

				unew.push_back(std::pair<int64_t, double>(longindex[1], u_new));
			}
			atu[idx[1]] += u_new / scale_a;
			atu[idx[0]] -= u_new / scale_a;
			atu[idx[2]] -= u_new / scale_a;

		}

		{
			double u_value = 0;
			if (uidx < u.size() && longindex[2] == u[uidx].first) {
				u_value = u[uidx].second;
				uidx++;
			}

			double axb, axbu;


			//calc y,ax
			axbu = (old_x[idx[2]] + b[idx[2]] - old_x[idx[1]] - b[idx[1]] - old_x[idx[0]] - b[idx[0]]) / scale_a;
			//calc y,ax
			axb = (x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]]) / scale_a;


			axbu += u_value / sigma;
			if (axbu > 0)
				axbu = 0;
			y_norm += square(axbu);
			axb -= axbu;
			axby += square(axb);
			double u_new = axb * tau * sigma + u_value;

			if (u_new != 0) {

				unew.push_back(std::pair<int64_t, double>(longindex[2], u_new));
			}
			atu[idx[2]] += u_new / scale_a;
			atu[idx[1]] -= u_new / scale_a;
			atu[idx[0]] -= u_new / scale_a;

		}



	}
	u.clear();
}
void admm::m_all_in_all( int n, vector<vector<pair<int64_t, double>>>& u, vector<vector<pair<int64_t, double>>>& u_new, Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, Eigen::VectorXd& old_x, double& axby, double& y_norm, double sigma, std::vector<int64_t>& cons, vector<int>&colidx, double tau = 1.0, double scale_a = 1.0) {
	int pnum = mthread;
	std::vector<std::thread>v;
	std::vector<Eigen::VectorXd>vatu;
	std::vector<double>vaxby(pnum), vy_norm(pnum);
	std::mutex mu;
	for (int i = 0; i < pnum; i++) {
		vatu.push_back(Eigen::VectorXd::Zero(n * (n - 1) / 2));
	}
	for (int i = 0; i < pnum; i++) {


		int va, vb;
		va = cons.size() / pnum * i;
		if (i != pnum - 1) {
			vb = cons.size() / pnum * (i + 1);
		}
		else {
			vb = cons.size();
		}






		v.push_back(std::thread(all_in_all, &mu, n, ref(u[i]),ref(u_new[i]), ref(x), ref(b), ref(vatu[i]), ref(old_x), ref(vaxby[i]), ref(vy_norm[i]), sigma, ref(cons),ref(colidx), tau, scale_a, va, vb));
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

void calc_u(int n, Eigen::VectorXd& x, Eigen::VectorXd& b, std::map<int64_t, double>& u, double sigma, int start = 0, int end = -1) {

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
					u[longindex[0]] = sigma * axbu2;
				}
				else if (flag) {
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

void atax(int n, Eigen::VectorXd& x, Eigen::VectorXd& res, std::vector<int64_t>& cons, int start = 0, int end = -1) {

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
void m_atax(int pnum, int n, Eigen::VectorXd& x, Eigen::VectorXd& res, std::vector<int64_t>& cons) {

	std::vector<std::thread>v;
	std::vector<Eigen::VectorXd>vatu;
	for (int i = 0; i < pnum; i++) {
		vatu.push_back(Eigen::VectorXd::Zero(n * (n - 1) / 2));
	}
	for (int i = 0; i < pnum; i++) {

		v.push_back(std::thread(atax, n, ref(x), ref(vatu[i]), cons, n - int((cbrt((i + 1) * 1.0 / pnum) * n)), n - int((cbrt(i * 1.0 / pnum) * 1.0 * n))));
	}
	for (int i = 0; i < pnum; i++) {
		v[i].join();
		res += vatu[i];
	}
}
void atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd& b, vector<pair<int64_t, double>>& u, Eigen::VectorXd& res, double sigma, std::vector<int64_t>& cons, vector<int>&colidx,double scale_a = 1, int start = 0, int end = -1) {

	if (end < 0)end = cons.size();
	int uidx = 0;
	for (int istart = start; istart < end; istart++) {
		int64_t lindex = cons[istart];
		int i, j, k;
		k = lindex % n;
		j = lindex / n;
		i = j / n;
		j = j % n;

		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];

		double axbu;

		double flag;
		double tem;
		int64_t longindex[3];
		longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
		longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
		longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;

		axbu = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;
		if (uidx<u.size()&&u[uidx].first == longindex[0]) {
			flag = u[uidx].second;
			uidx++;
		}
		else {
			flag = 0;
		}

		


			axbu += flag / sigma;
			
		

		if (axbu > 0) {
			axbu = 0;
		}
		tem = axbu - (b[idx[0]] - b[idx[1]] - b[idx[2]]) / scale_a;

			tem -= flag / sigma;


		res[idx[0]] += tem / scale_a;
		res[idx[1]] -= tem / scale_a;
		res[idx[2]] -= tem / scale_a;

		//second scene
		axbu = (x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]]) / scale_a;
		if (uidx < u.size() && u[uidx].first == longindex[1]) {
			flag = u[uidx].second;
			uidx++;
		}
		else {
			flag = 0;
		}
			axbu += flag/ sigma;

		if (axbu > 0) {
			axbu = 0;
		}
		tem = axbu - (b[idx[1]] - b[idx[0]] - b[idx[2]]) / scale_a;
		
			tem -= flag/ sigma;
		

		res[idx[1]] += tem / scale_a;
		res[idx[0]] -= tem / scale_a;
		res[idx[2]] -= tem / scale_a;
		//third scene
		axbu = (x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]]) / scale_a;
		if (uidx < u.size() && u[uidx].first == longindex[2]) {
			flag = u[uidx].second;
			uidx++;
		}
		else {
			flag = 0;
		}
			axbu += flag/ sigma;

		if (axbu > 0) {
			axbu = 0;
		}
		tem = axbu - (b[idx[2]] - b[idx[1]] - b[idx[0]]) / scale_a;
		
			tem -= flag/ sigma;
		

		res[idx[2]] += tem / scale_a;
		res[idx[1]] -= tem / scale_a;
		res[idx[0]] -= tem / scale_a;
	}
}
void admm::m_atyu(int n, Eigen::VectorXd& x, Eigen::VectorXd& b, vector<vector<pair<int64_t, double>>>&u, Eigen::VectorXd& res, double sigma, std::vector<int64_t>& cons,vector<int>&colidx, double scale_a = 1) {
	int pnum = mthread;
	std::vector<std::thread>v;
	std::vector<Eigen::VectorXd>vatu;
	for (int i = 0; i < pnum; i++) {
		vatu.push_back(Eigen::VectorXd::Zero(n * (n - 1) / 2));
	}
	for (int i = 0; i < pnum; i++) {
		int va, vb;
		va = cons.size() / pnum * i;
		if (i != pnum - 1) {
			vb = cons.size() / pnum * (i + 1);
		}
		else {
			vb = cons.size();
		}
		v.push_back(std::thread(atyu, n, ref(x), ref(b), ref(u[i]), ref(vatu[i]), sigma, ref(cons),ref(colidx), scale_a, va, vb));
	}
	for (int i = 0; i < pnum; i++) {
		v[i].join();
		res += vatu[i];
	}


}
admm_basic::admm_basic(string graph_name,int type) :admm(graph_name,type)
{



}

admm_basic::~admm_basic()
{

}
 double admm_basic::func_p(Eigen::VectorXd x,Eigen::VectorXd w) {
	return w.cwiseProduct(x.cwiseAbs()).sum();
}
inline Eigen::VectorXd admm_basic::prox_p(Eigen::VectorXd x, double sigma,Eigen::VectorXd w) {
	return (x.cwiseAbs() - w / sigma).cwiseMax(0.0).cwiseProduct(x.cwiseSign());
}
double admm_basic::get_kkt(double sigma) {

	return 0.0;

	//cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;
}

void admm::constrain_generate(int istart,int iend,mutex&mut,vector<int>&colidx,int&colidx_num, std::vector<Eigen::Triplet<double> >&acons, std::vector<int64_t>&cons,std::set<int64_t>&viscons){

	for (int i = istart; i < iend; i++) {
		int64_t index;
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				index = static_cast<int64_t>(i) * n * n + j * n + k;

				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 3e-2) {

					if (viscons.find(index) == viscons.end()) {
						mut.lock();
						for (int tp = 0; tp < 3; tp++) {
							if (colidx[idx[tp]] == -1) {

								colidx[idx[tp]] = colidx_num;
								colidx_num++;

							}
						}
						cons.push_back(index);

						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[0]], 1));
						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[1]], -1));
						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[2]], -1));

						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[0]], -1));
						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[1]], 1));
						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[2]], -1));

						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[0]], -1));
						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[1]], -1));
						acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[2]], 1));



						viscons.insert(index);
						mut.unlock();
					}
				}
			}
		}
	}


}




double admm_basic::run(double tau, double sigma, double scale_b, double scale_c) {
	vector<vector<pair<int64_t, double>>>u1(mthread),u2(mthread);
	vector<vector<pair<int64_t, double>>>&u = u1;
	vector<vector<pair<int64_t, double>>>&u_new = u2;
	double scale_a;
	//scale_a = 3;
	//scale_b = b.norm();
	//b = b / b.norm();
	//w = w / 9;
	//scale_b = b.norm();
	//scale_c = 3;

	std::vector<int> colidx(var_num,-1);
	int colidx_num=0;
	// = b / b.norm();
	scale_a = 3;
	b = b * 3 / scale_b;
	w = w / scale_c;
	Eigen::VectorXd rhs(var_num);

	int dual_win(0), prim_win(0);
	int maxsize = 0;
	std::vector<int64_t>cons;
	std::set<int64_t>viscons;
	int old_size = viscons.size();
	std::vector<Eigen::Triplet<double> >acons;

	timer madoka,homula;
	double finalres;
	int loopflag = 0;

	while (true) {

		
		if (true) {


			for (int i = 0; i < n; i++) {
				int64_t index;
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						index = static_cast<int64_t>(i) * n * n + j * n + k;

						if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
							x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
							x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 3e-2) {

							if (viscons.find(index) == viscons.end()) {
								for (int tp = 0; tp < 3; tp++) {
									if (colidx[idx[tp]] == -1) {

										colidx[idx[tp]] = colidx_num;
										colidx_num++;

									}
								}
								cons.push_back(index);

								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[0]], 1));
								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[1]], -1));
								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[2]], -1));

								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[0]], -1));
								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[1]], 1));
								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[2]], -1));

								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[0]], -1));
								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[1]], -1));
								acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[2]], 1));



								viscons.insert(index);
							}
						}
					}
				}
			}
		}
		else {
			vector<thread>vt;
			mutex mut;
			for (int i = 0; i < mthread; i++) {

				int ta, tb;
				ta = n - (int)(cbrt((i + 1) * 1.0 / mthread) * n);
				tb = n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n);
				vt.push_back(thread(&admm::constrain_generate, this,ta,tb, ref(mut), ref(colidx), ref(colidx_num), ref(acons), ref(cons), ref(viscons)));
			}
			for (int i = 0; i < mthread; i++) {
				vt[i].join();
			}
		}
		cout << "cons:" << viscons.size() << endl;

		if (loopflag==2&&old_size == viscons.size()) {


			cout << "timeis:" << madoka.time() << " homulais:" << homula.time() << endl;

			x = x / scale_a * scale_b;
			b = b / 3;
			return finalres;
		}
		else if (loopflag == 2||old_size == viscons.size()) {
			loopflag=1;
		}

		old_size = viscons.size();
		Eigen::SparseMatrix<double> hst(viscons.size() * 3, colidx_num);
		Eigen::SparseMatrix<double> hs(colidx_num, colidx_num);
		hst.setFromTriplets(acons.begin(), acons.end());
		//acons.clear();

		hs = hst.transpose() * hst / (scale_a * scale_a);

		hs += Eigen::VectorXd::Ones(colidx_num).asDiagonal();
		Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
		//Eigen::PardisoLDLT<SparseMatrix<double>> cg;
		cg.compute(hs);

		sub_x = Eigen::VectorXd::Zero(colidx_num);
		sub_b = Eigen::VectorXd::Zero(colidx_num);
		sub_w = Eigen::VectorXd::Zero(colidx_num);
		sub_v = Eigen::VectorXd::Zero(colidx_num);



		for (int i = 0; i < var_num;i++) {
			if (colidx[i] != -1) {


				sub_x[colidx[i]] = x[i];
				sub_b[colidx[i]] = b[i];
				sub_w[colidx[i]] = w[i];
				sub_v[colidx[i]] = v[i];
			}
		} 


		int sub_var = colidx_num;

		for (int iter_I = 0;; iter_I++)
		{
			
			

			madoka.tic();
			Eigen::VectorXd old_x = sub_x;
			rhs = Eigen::VectorXd::Zero(sub_var);
			//atyu(n, x, b, &u, rhs, sigma,scale_a);
			//cout << u.size() << endl;
			m_atyu( n, sub_x, sub_b, u, rhs, sigma, cons, colidx, scale_a);
			rhs += prox_p(sub_x + sub_v / sigma, sigma,sub_w) - sub_v / sigma;
			//cg
			//cout << "atyu\n";
			madoka.toc();
			double tol = 1e-10;
			Eigen::VectorXd r = rhs;
			Eigen::VectorXd p(sub_var);
			Eigen::VectorXd res = Eigen::VectorXd::Zero(sub_var);
			double row, row_old;
			row = r.squaredNorm();
			double delta = tol * rhs.norm();
			//cout << "cgstart\n";
			//===============================CG===========================
			//{

			//	int i = 0;
			//	for (; sqrt(row) > delta; i++) {
			//		if (i == 0) {
			//			p = r;
			//		}
			//		else {
			//			p = r + row / row_old * p;
			//		}

			//		Eigen::VectorXd ap = Eigen::VectorXd::Zero(var_num);
			//		m_atax(24, n, p, ap,cons);
			//		ap /= (scale_a * scale_a);
			//		ap += p;

			//		double mu = row / p.dot(ap);
			//		res += mu * p;
			//		r -= mu * ap;
			//		row_old = row;
			//		row = r.squaredNorm();
			//	}
			//	//cout << "  i:" << i << "  ";
			//}
			//=================================CGEND=================================
			// 
			// 
			// 
			//cout << "cgend\n";
			//cg end

			//cout << "res:" << (rhs - test).norm() << endl;


			//double ita = std::fmin(0.5, 0.1 * rhs.norm());
			//cg.setTolerance(ita / rhs.norm());
			homula.tic();
			sub_x = cg.solve(rhs);
			homula.toc();

			double y_norm, kkt_u;
			y_norm = kkt_u = 0;
			Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_var);



			m_all_in_all( n, u,u_new, sub_x, sub_b, atu, old_x, kkt_u, y_norm, sigma, cons, colidx, tau, scale_a);
			std::swap(u, u_new);

			Eigen::VectorXd z = prox_p(sub_x + sub_v / sigma, sigma,sub_w);
			sub_v += sigma * tau * (sub_x - z);

			double kkt_x = (atu + sub_v).norm() / (1 + sub_v.norm());
			double true_ku = scale_b * sqrt(kkt_u) / (1 + scale_b * sqrt(y_norm));
			kkt_u = sqrt(kkt_u) / (1 + sqrt(y_norm));
			double kkt_v = (sub_x - z).norm() / (1 + z.norm());

			//cout << "allinall\n";
			if (std::fmax(kkt_u, kkt_v) > kkt_x) {
				dual_win++;
			}
			else {
				prim_win++;
			}
			if ((iter_I + 1) % itertime(iter_I) == 0) {
				if (prim_win > dual_win * 1.2 && prim_win > 1) {
					prim_win = 0;
					sigma /= 3;
				}
				else if (dual_win > prim_win * 1.2 && dual_win > 1) {
					dual_win = 0;
					sigma *= 3;
				}
			}
			kkt_x = scale_c*( atu + sub_v).norm() / (1 + scale_c*sub_v.norm());
			//double kkt_z=scale_b/scale_a*(z - prox_p(z + sub_v, 1.0, sub_w)).norm() / (1 + scale_b / scale_a * z.norm());
			kkt_v = scale_b / scale_a * (sub_x - z).norm() / (1 + scale_b / scale_a * z.norm());

			//cout << kkt_x << " " << kkt_v << " " << true_ku <<" !:"<<kkt_z<< endl;
			double kkt = std::fmax(kkt_x, std::fmax(true_ku, kkt_v));
			maxsize = fmax(maxsize, u.size());

			cout << iter_I << "\t" << kkt << "\t" << func_p(sub_x,sub_w) * scale_b*scale_c / 3 << "\t" << sigma << "\t" << u.size() << endl;
			if (kkt < 1e-3) {
				

				double pobj = func_p(sub_x,sub_w) * scale_b / 3;
				double dobj = 0;
#pragma omp parallel for reduction(+:dobj)
				for (int i = 0; i < u.size(); i++) {


					for (auto id : u[i]) {
						int64_t longindex = id.first;
						int i, j, k;
						k = longindex % n;
						j = longindex / n;
						i = j / n;
						j = j % n;
						int idx[3];
						idx[0] = colidx[ij(i, j, n)];
						idx[1] = colidx[ij(i, k, n)];
						idx[2] = colidx[ij(j, k, n)];
						dobj += id.second * (sub_b[idx[0]] - sub_b[idx[1]] - sub_b[idx[2]]);
					}
				}
				dobj /= scale_a * scale_a;
				dobj *= scale_b * scale_c;
				cout <<"pobj:"<<pobj<< " dobj:" << dobj << endl;
				if (fabs(pobj - dobj) / dobj < 1e-2) {

					if (loopflag == 0) {

						for (int i = 0; i < var_num; i++) {
							if (colidx[i] != -1) {


								 x[i]= sub_x[colidx[i]];
								 b[i]=sub_b[colidx[i]];
								 w[i]=sub_w[colidx[i]];
								 v[i]=sub_v[colidx[i]];
							}
						}

						//for (auto i = colidx.begin(); i != colidx.end(); i++) {
						//	x[i->first] = sub_x[i->second];
						//	b[i->first] = sub_b[i->second];
						//	w[i->first] = sub_w[i->second];
						//	v[i->first] = sub_v[i->second];
						//}
						finalres = func_p(sub_x, sub_w) * scale_b * scale_c / scale_a;
						break;
					}
					else {
						bool flg = true;
						for (auto iter : cons) {
							int64_t longindex = iter;
							int i, j, k;
							k = longindex % n;
							j = longindex / n;
							i = j / n;
							j = j % n;
							int idx[3];
							idx[0] = colidx[ij(i, j, n)];
							idx[1] = colidx[ij(i, k, n)];
							idx[2] = colidx[ij(j, k, n)];
							if (sub_x[idx[0]] + sub_b[idx[0]] - sub_x[idx[1]] - sub_b[idx[1]] - sub_x[idx[2]] - sub_b[idx[2]] > 3e-2 ||
								sub_x[idx[1]] + sub_b[idx[1]] - sub_x[idx[0]] - sub_b[idx[0]] - sub_x[idx[2]] - sub_b[idx[2]] > 3e-2 ||
								sub_x[idx[2]] + sub_b[idx[2]] - sub_x[idx[1]] - sub_b[idx[1]] - sub_x[idx[0]] - sub_b[idx[0]] > 3e-2) {
								flg = false;
								break;
							}
						}
						if (flg) {
							loopflag = 2;

							for (int i = 0; i < var_num; i++) {
								if (colidx[i] != -1) {


									x[i] = sub_x[colidx[i]];
									b[i] = sub_b[colidx[i]];
									w[i] = sub_w[colidx[i]];
									v[i] = sub_v[colidx[i]];
								}
							}
							finalres = func_p(sub_x, sub_w) * scale_b * scale_c / scale_a;
							break;
						}

					}
				}
			}
			//cout << "sig:"<<sigma << endl;
			//cout << rd << "\t" << rp << endl;

		}
	}
}

void admm_basic::run_dense() {
	Eigen::SparseMatrix<double>A(n * (n - 1) * (n - 2) / 2, var_num);

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
		ata += tool.asDiagonal();
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
		y = A * (x + b) + vu / sigma;
		y = y.cwiseMin(0.0);


		z = prox_p(x + v / sigma, sigma,w);

		rhs = A.transpose() * (y + vb - vu / sigma) + z - v / sigma;

		x = ll.solve(rhs);
		y = A * x - vb + vu / sigma;
		y = y.cwiseMin(0.0);

		vu += sigma * (A * x - y - vb);
		z = prox_p(x + v / sigma, sigma,w);
		v += sigma * (x - z);


		double rd = (A.transpose() * vu + v).norm() / (1 + v.norm());
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
/*
double admm_base::run(double tau, double sigma, double scale_b, double scale_c) {


	double scale_a;
	//scale_a = 3;
	//scale_b = b.norm();
	//b = b / b.norm();
	//w = w / 9;
	//scale_b = b.norm();

	b = b * 3 / scale_b;

	// = b / b.norm();
	scale_a = 3;

	w = w / scale_c;
	Eigen::VectorXd rhs(var_num);

	int dual_win(0), prim_win(0);
	int maxsize = 0;
	std::vector<int64_t>cons;
	std::set<int64_t>viscons;
	int old_size = viscons.size();
	std::vector<Eigen::Triplet<double> >acons;


	std::pair<int, int>finalres;

	while (true) {



		for (int i = 0; i < n; i++) {
			int64_t index;
			for (int j = i + 1; j < n; j++) {
				for (int k = j + 1; k < n; k++) {
					int idx[3];
					idx[0] = ij(i, j, n);
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					index = static_cast<int64_t>(i) * n * n + j * n + k;

					if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 1e-3 ||
						x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 1e-3 ||
						x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 1e-3) {

						if (viscons.find(index) == viscons.end()) {
							cons.push_back(index);

							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, idx[0], 1));
							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, idx[1], -1));
							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, idx[2], -1));

							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, idx[0], -1));
							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, idx[1], 1));
							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, idx[2], -1));

							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, idx[0], -1));
							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, idx[1], -1));
							acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, idx[2], 1));



							viscons.insert(index);
						}
					}
				}
			}
		}
		cout << "cons:" << viscons.size() << endl;

		if (old_size == viscons.size()) {
			break;
		}
		old_size = viscons.size();
		Eigen::SparseMatrix<double> hst(viscons.size() * 3, var_num);
		Eigen::SparseMatrix<double> hs(var_num, var_num);
		hst.setFromTriplets(acons.begin(), acons.end());
		hs = hst.transpose() * hst / (scale_a * scale_a);

		hs += w.cwiseProduct(w).asDiagonal();
		Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;

		cg.compute(hs);




		for (int k = 0;; k++)
		{
			Eigen::VectorXd z = w.cwiseProduct(x) + v / sigma - pi_c(w.cwiseProduct(x) + v / sigma, sigma);
			Eigen::VectorXd old_x = x;
			rhs = Eigen::VectorXd::Zero(var_num);
			//atyu(n, x, b, &u, rhs, sigma,scale_a);
			//m_atyu(32, n, x, b, &u, rhs, sigma, cons, scale_a);
			rhs += w.cwiseProduct(z -v / sigma);
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
			//===============================CG===========================
			//{

			//	int i = 0;
			//	for (; sqrt(row) > delta; i++) {
			//		if (i == 0) {
			//			p = r;
			//		}
			//		else {
			//			p = r + row / row_old * p;
			//		}

			//		Eigen::VectorXd ap = Eigen::VectorXd::Zero(var_num);
			//		m_atax(24, n, p, ap,cons);
			//		ap /= (scale_a * scale_a);
			//		ap += p;

			//		double mu = row / p.dot(ap);
			//		res += mu * p;
			//		r -= mu * ap;
			//		row_old = row;
			//		row = r.squaredNorm();
			//	}
			//	//cout << "  i:" << i << "  ";
			//}
			//=================================CGEND=================================
			// 
			// 
			// 
			//cout << "cgend\n";
			//cg end

			//cout << "res:" << (rhs - test).norm() << endl;


			//double ita = std::fmin(0.5, 0.1 * rhs.norm());
			//cg.setTolerance(ita / rhs.norm());

			x = cg.solve(rhs);


			double y_norm, kkt_u;
			y_norm = kkt_u = 0;
			Eigen::VectorXd atu = Eigen::VectorXd::Zero(var_num);
			//m_all_in_all(32, n, u, x, b, atu, old_x, kkt_u, y_norm, sigma, cons, tau, scale_a);
			z = w.cwiseProduct(x)+v/sigma-pi_c(w.cwiseProduct(x) + v / sigma, sigma);
			v += sigma * tau * (w.cwiseProduct(x) - z);

			double kkt_x = (atu + w.cwiseProduct(v)).norm() / (1 + w.cwiseProduct(v).norm());
			double true_ku = scale_b * sqrt(kkt_u) / (1 + scale_b * sqrt(y_norm));
			kkt_u = sqrt(kkt_u) / (1 + sqrt(y_norm));
			double kkt_v = (w.cwiseProduct(x) - z).norm() / (1 + z.norm());

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
			kkt_x = scale_a * scale_a * (atu + w.cwiseProduct(v)).norm() / (1 + scale_a * scale_a * w.cwiseProduct(v).norm());

			kkt_v = scale_b / scale_a * (w.cwiseProduct(x) - z).norm() / (1 + scale_b / scale_a * z.norm());


			double kkt = std::fmax(kkt_x, std::fmax(true_ku, kkt_v));
			maxsize = fmax(maxsize, u.size());
			cout << k << "\t" << kkt << "\t" << func_p(w.cwiseProduct(x)) * scale_b * 3 << "\t" << sigma << "\t" << u.size() << endl;
			if (kkt < 1e-3) {


				finalres = std::pair<int, int>(k, maxsize);
				break;
			}

			//cout << "sig:"<<sigma << endl;
			//cout << rd << "\t" << rp << endl;

		}

	}
	return func_p(w.cwiseProduct(x));

}
*/
inline double admm_base::func_p(Eigen::VectorXd x) {
	return x.cwiseAbs().sum();
}
inline Eigen::VectorXd admm_base::pi_c(Eigen::VectorXd x, double sigma) {
	return x.cwiseAbs().cwiseMin(1 / sigma).cwiseProduct(x.cwiseSign());
}

admm_base::admm_base(string graph_name) :metric(graph_name)
{



}

admm_base::~admm_base()
{

}

admm_2::admm_2(string graph_name):admm_base(graph_name)
{
	w = w.cwiseSqrt();
}

admm_2::~admm_2()
{
}

double admm_2::func_p(Eigen::VectorXd x)
{
	return x.norm();
}

Eigen::VectorXd admm_2::pi_c(Eigen::VectorXd x, double sigma)
{
	if (x.norm() <= 1 / sigma)
		return x;
	return x * (1 / sigma) / x.norm();
}
