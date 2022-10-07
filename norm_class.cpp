#include"mytime.h"
#include"norm_class.h"
#include<iomanip>
#include<atomic>
#include"cons_generate.h"
//#define EIGENCG
//#define SMCG
#define BDEG

#define DIAGMULT


constexpr int searchnum =10;
constexpr double tol = 1e-4;
int scale_a;

void new_atu(concurrency::concurrent_vector<cons_full>::iterator beg, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& res,int scale_a) {
	for (auto i = beg; i != end; i++) {
		int* idx = i->idx;
		double tem = i->u;
		res[idx[0]] += tem / scale_a;
		res[idx[1]] -= tem / scale_a;
		res[idx[2]] -= tem / scale_a;
	}
}


template <typename T>
void buildeigen(T&data,Eigen::SparseMatrix<double,Eigen::RowMajor>&hs,Eigen::VectorXd&diag) {
	timer bella;
	hs.reserve(2 * diag.cast<Eigen::Index>());
	auto dt = hs.valuePtr();
	auto nnz = hs.innerNonZeroPtr();
	auto inidx = hs.innerIndexPtr();
	auto oidx = hs.outerIndexPtr();




	for (int i = 0; i < data.size(); i++) {
		int* idx = data[i];
		//hs.coeffRef(idx[0], idx[1]) += -1;
		//hs.coeffRef(idx[1], idx[0]) += -1;
		//hs.coeffRef(idx[0], idx[2]) += -1;
		//hs.coeffRef(idx[2], idx[0]) += -1;
		//hs.coeffRef(idx[1], idx[2]) += 1;
		//hs.coeffRef(idx[2], idx[1]) += 1;
		dt[oidx[idx[0]] + (nnz[idx[0]])] = -1;
		inidx[oidx[idx[0]] + (nnz[idx[0]]++)] = idx[1];
		dt[oidx[idx[0]] + (nnz[idx[0]])] = -1;
		inidx[oidx[idx[0]] + (nnz[idx[0]]++)] = idx[2];
		dt[oidx[idx[1]] + (nnz[idx[1]])] = -1;
		inidx[oidx[idx[1]] + (nnz[idx[1]]++)] = idx[0];
		dt[oidx[idx[1]] + (nnz[idx[1]])] = 1;
		inidx[oidx[idx[1]] + (nnz[idx[1]]++)] = idx[2];
		dt[oidx[idx[2]] + (nnz[idx[2]])] = -1;
		inidx[oidx[idx[2]] + (nnz[idx[2]]++)] = idx[0];
		dt[oidx[idx[2]] + (nnz[idx[2]])] = 1;
		inidx[oidx[idx[2]] + (nnz[idx[2]]++)] = idx[1];
	}




}



template<typename T>
int psqmr_inf_diag(T& ata, Eigen::VectorXd& td1,Eigen::VectorXd& td2, Eigen::VectorXd& diag, Eigen::VectorXd& rnew, Eigen::VectorXd& rhs, Eigen::VectorXd& direction,double sigma, double delta = 1e-4) {
	//delta = 1e-4;
	Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd r = rhs;
	double rou = r.squaredNorm();
	//double delta = std::fmin(0.5, 0.0001 * std::fmax(rhs.norm(),1.0));
	//double delta = 1e-4;
	Eigen::VectorXd q = r.cwiseProduct(diag);
	Eigen::VectorXd res = r;
	Eigen::VectorXd aq, u;
	Eigen::VectorXd d = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd ad = Eigen::VectorXd::Zero(rhs.size());
	double tau;
	double rou_old = r.dot(q);
	double tau_old = q.norm(), theta_old = 0;
	double err = r.norm();
	int k;
	for (k = 0; k < 10000; k++) {
		aq = sigma*td1.cwiseProduct(ata * td1.cwiseProduct(q))/(scale_a*scale_a) + td2.cwiseProduct(q) - rnew.dot(q) * rnew;
		double _sigma = aq.dot(q);
		if (_sigma < 1e-30)break;
		double alpha = rou_old / _sigma;
		r -= alpha * aq;
		u = r.cwiseProduct(diag);
		double theta = u.norm() / tau_old;
		double c = 1 / (sqrt(1 + theta * theta));
		tau = tau_old * theta * c;
		double gam = (c * c * theta_old * theta_old);
		double eta = c * c * alpha;
		d = gam * d + eta * q;

		x += d;
		ad = gam * ad + eta * aq;
		res -= ad;
		err = res.norm();
		if (err < delta && k>5) {
			//cout << err << " " << delta << " md" << endl;
			break;
		}

		//if (k > 20) {    Õ£÷ÕºÏ≤È
		//	double ratio=
		//}
		if (abs(rou_old < 1e-30)) {
			break;
		}
		rou = r.dot(u);
		q = rou / rou_old * q + u;
		rou_old = rou;
		tau_old = tau;
		theta_old = theta;
	}
	direction = -x;
	return k;
}


template<typename T>
int psqmr_ori(T& ata, Eigen::VectorXd& diag, Eigen::VectorXd& rhs, Eigen::VectorXd& direction, double sigma, double delta = 1e-4) {
	//delta = 1e-4;
	Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd r = rhs;
	double rou = r.squaredNorm();
	//double delta = std::fmin(0.5, 0.0001 * std::fmax(rhs.norm(),1.0));
	//double delta = 1e-4;
	Eigen::VectorXd q = r.cwiseProduct(diag);
	Eigen::VectorXd res = r;
	Eigen::VectorXd aq, u;
	Eigen::VectorXd d = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd ad = Eigen::VectorXd::Zero(rhs.size());
	double tau;
	double rou_old = r.dot(q);
	double tau_old = q.norm(), theta_old = 0;
	double err = r.norm();
	int k;
	for (k = 0; k < 10000; k++) {
		aq = ata*q;
		double _sigma = aq.dot(q);
		if (_sigma < 1e-30)break;
		double alpha = rou_old / _sigma;
		r -= alpha * aq;
		u = r.cwiseProduct(diag);
		double theta = u.norm() / tau_old;
		double c = 1 / (sqrt(1 + theta * theta));
		tau = tau_old * theta * c;
		double gam = (c * c * theta_old * theta_old);
		double eta = c * c * alpha;
		d = gam * d + eta * q;

		x += d;
		ad = gam * ad + eta * aq;
		res -= ad;
		err = res.norm();
		if (err < delta && k>5) {
			//cout << err << " " << delta << " md" << endl;
			break;
		}

		//if (k > 20) {    Õ£÷ÕºÏ≤È
		//	double ratio=
		//}
		if (abs(rou_old < 1e-30)) {
			break;
		}
		rou = r.dot(u);
		q = rou / rou_old * q + u;
		rou_old = rou;
		tau_old = tau;
		theta_old = theta;
	}
	direction = -x;
	return k;
}
template<typename T>
int psqmr_ori_diag(T& ata, Eigen::VectorXd& td2, Eigen::VectorXd& diag, Eigen::VectorXd& rhs, Eigen::VectorXd& direction, double sigma, double delta = 1e-4) {
	//delta = 1e-4;
	Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd r = rhs;
	double rou = r.squaredNorm();
	//double delta = std::fmin(0.5, 0.0001 * std::fmax(rhs.norm(),1.0));
	//double delta = 1e-4;
	Eigen::VectorXd q = r.cwiseProduct(diag);
	Eigen::VectorXd res = r;
	Eigen::VectorXd aq, u;
	Eigen::VectorXd d = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd ad = Eigen::VectorXd::Zero(rhs.size());
	double tau;
	double rou_old = r.dot(q);
	double tau_old = q.norm(), theta_old = 0;
	double err = r.norm();
	int k;
	for (k = 0; k < 10000; k++) {
		aq = sigma * ata*q/ (scale_a * scale_a) + td2.cwiseProduct(q);
		double _sigma = aq.dot(q);
		if (_sigma < 1e-30)break;
		double alpha = rou_old / _sigma;
		r -= alpha * aq;
		u = r.cwiseProduct(diag);
		double theta = u.norm() / tau_old;
		double c = 1 / (sqrt(1 + theta * theta));
		tau = tau_old * theta * c;
		double gam = (c * c * theta_old * theta_old);
		double eta = c * c * alpha;
		d = gam * d + eta * q;

		x += d;
		ad = gam * ad + eta * aq;
		res -= ad;
		err = res.norm();
		if (err < delta && k>5) {
			//cout << err << " " << delta << " md" << endl;
			break;
		}

		//if (k > 20) {    Õ£÷ÕºÏ≤È
		//	double ratio=
		//}
		if (abs(rou_old < 1e-30)) {
			break;
		}
		rou = r.dot(u);
		q = rou / rou_old * q + u;
		rou_old = rou;
		tau_old = tau;
		theta_old = theta;
	}
	direction = -x;
	return k;
}



template<typename T>
int psqmr_inf(T& ata,Eigen::VectorXd&diag, Eigen::VectorXd& rnew, Eigen::VectorXd& rhs, Eigen::VectorXd& direction,double delta=1e-4) {
	//delta = 1e-4;
	Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd r = rhs;
	double rou = r.squaredNorm();
	//double delta = std::fmin(0.5, 0.0001 * std::fmax(rhs.norm(),1.0));
	//double delta = 1e-4;
	Eigen::VectorXd q = r.cwiseProduct(diag);
	Eigen::VectorXd res = r;
	Eigen::VectorXd aq, u;
	Eigen::VectorXd d = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd ad = Eigen::VectorXd::Zero(rhs.size());
	double tau;
	double rou_old = r.dot(q);
	double tau_old = q.norm(), theta_old = 0;
	double err = r.norm();
	int k;
	for (k = 0; k < 10000; k++) {
		aq = ata * q - rnew.dot(q) * rnew;
		double _sigma = aq.dot(q);
		if (_sigma < 1e-30)break;
		double alpha = rou_old / _sigma;
		r -= alpha * aq;
		u = r.cwiseProduct(diag);
		double theta = u.norm() / tau_old;
		double c = 1 / (sqrt(1 + theta * theta));
		tau = tau_old * theta * c;
		double gam = (c * c * theta_old * theta_old);
		double eta = c * c * alpha;
		d = gam * d + eta * q;

		x += d;
		ad = gam * ad + eta * aq;
		res -= ad;
		err = res.norm();
		if (err < delta && k>5) {
			//cout << err << " " << delta << " md" << endl;
			break;
		}

		//if (k > 20) {    Õ£÷ÕºÏ≤È
		//	double ratio=
		//}
		if (abs(rou_old < 1e-30)) {
			break;
		}
		rou = r.dot(u);
		q = rou / rou_old * q + u;
		rou_old = rou;
		tau_old = tau;
		theta_old = theta;
	}
	direction = -x;
	return k;
}

template<typename T>
int psqmr(T& ata,Eigen::VectorXd& rnew , Eigen::VectorXd& rhs, Eigen::VectorXd& direction) {

	Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd r = rhs;
	double rou = r.squaredNorm();
	//double delta = std::fmin(0.5, 0.000001 * rhs.norm());
	double delta = std::fmin(1e-3,1e-3*rhs.norm());
	Eigen::VectorXd q = r;
	Eigen::VectorXd res = r;
	Eigen::VectorXd aq, u;
	Eigen::VectorXd d = Eigen::VectorXd::Zero(rhs.size());
	Eigen::VectorXd ad = Eigen::VectorXd::Zero(rhs.size());
	double tau, rou_last, mu;
	double rou_old = r.squaredNorm();
	double tau_old = q.norm(), theta_old = 0;
	double err = r.norm();
	int k;
	for (k = 0; k < 1000; k++) {
		aq = ata * q-rnew.dot(q)*rnew;
		double _sigma = aq.dot(q);
		if (_sigma < 1e-30)break;
		double alpha = rou_old / _sigma;
		r -= alpha * aq;
		u = r;
		double theta = u.norm() / tau_old;
		double c = 1 / (sqrt(1 + theta * theta));
		double tau = tau_old * theta * c;
		double gam = (c * c * theta_old * theta_old);
		double eta = c * c * alpha;
		d = gam * d + eta * q;

		x += d;
		ad = gam * ad + eta * aq;
		res -= ad;
		err = res.norm();
		if (err < delta && k>5) {
			//cout << err << " " << delta << " md" << endl;
			break;
		}
		//if (k > 20) {    Õ£÷ÕºÏ≤È
		//	double ratio=
		//}
		if (abs(rou_old < 1e-30)) {
			break;
		}
		rou = r.dot(u);
		q = rou / rou_old * q + u;
		rou_old = rou;
		tau_old = tau;
		theta_old = theta;
	}
	direction = -x;
	return k;
}



norm_ori::norm_ori(string graph_name) :metric(graph_name)
{
	mthread = std::thread::hardware_concurrency();
	iter_all = 0;

	u.resize(mthread);
}

norm_ori::~norm_ori()
{
}

inline Eigen::VectorXd norm_ori::prox_p_star(Eigen::VectorXd x, double sigma,Eigen::VectorXd w) {
	return x.cwiseAbs().cwiseMin(w / sigma).cwiseProduct(x.cwiseSign());
}
inline Eigen::VectorXd norm_ori::prox_p(Eigen::VectorXd x, double sigma, Eigen::VectorXd w) {
	return (x.cwiseAbs() - w / sigma).cwiseMax(0.0).cwiseProduct(x.cwiseSign());
}
inline double norm_ori::func_p(Eigen::VectorXd x, Eigen::VectorXd w) {
	return w.cwiseProduct(x.cwiseAbs()).sum();
}

inline Eigen::VectorXd norm_ori::prox_p_star(Eigen::VectorXd x, double sigma) {
	return x.cwiseAbs().cwiseMin(w / sigma).cwiseProduct(x.cwiseSign());
}
inline Eigen::VectorXd norm_ori::prox_p(Eigen::VectorXd x, double sigma) {
	return (x.cwiseAbs() - w / sigma).cwiseMax(0.0).cwiseProduct(x.cwiseSign());
}
inline double norm_ori::func_p(Eigen::VectorXd x) {
	return w.cwiseProduct(x.cwiseAbs()).sum();
}






double norm_inf::run_axbu_sparse(double cst_const, int con_const, int sub_iter) {
	
	scale_a = 6;
	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();
	sv_b.default_value = scale_a;

	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 0.01;

	double c_t = 0;
	//std::unordered_map<int, int>colidx;

	timer time1, time2, time3, time4, time5;

	cout << "m:" << mthread << endl;
	int mxtnum = thread::hardware_concurrency();
	cons_generator_sparse generator(n, mxtnum);
	generator.consvio_limit = scale_a * 0.01;
	double total_kkt_u = 0, total_y_norm = 0;

	double sigini = 1e-2;

	bool consflag = false;
	while (true) {

		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(sv_x, sv_b, scale_a * tol,w_value_map);
		}

		double totalkkt = sqrt(generator.k1 + total_kkt_u) / (1 + sqrt(generator.k2 + total_y_norm));
		cout << "totalkkt:" << totalkkt << endl;
		if (generator.consvio)
			break;

		concurrency::concurrent_unordered_map<int,int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = sigini;
		cst_const = sigma * 1e-3;
		cout << "cst:" << cst_const << endl;
		int sub_var = colidx_num;
		std::cout << "total num of variables:" << var_num << " num of current variables:" << sub_var << endl;

		auto c_end = std::chrono::high_resolution_clock::now();
		auto c_eld = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start);
		c_t += double(c_eld.count()) * 1.0 * 1e-9;

		double dobj = 0;



		double	primal_win = 0;
		double dual_win = 0;
		double kkt0 = 1;
		kkt = 0.1;

		Eigen::VectorXd sub_x = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_b = Eigen::VectorXd::Ones(colidx_num);
		Eigen::VectorXd sub_w = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_v = Eigen::VectorXd::Zero(colidx_num);



		for (auto& i : colidx) {

			sub_x[i.second] = sv_x[i.first] * sv_w[i.first];
			sub_b[i.second] = sv_b[i.first];
			sub_w[i.second] = sv_w[i.first];
			sub_v[i.second] = sv_v[i.first];
		}
		Eigen::VectorXd inv_w = Eigen::VectorXd::Ones(sub_w.size()).cwiseQuotient(sub_w);
		Eigen::VectorXd ori_x = inv_w.cwiseProduct(sub_x);
		double kkt = 0.1;
		cout << sub_var << endl;
		double tol_cg = 5e-2;

		double min_kkt = 0.1;
		rlim = 1e-1;
		for (iter_all = 0; iter_all < 1000; iter_all++) {


			x0 = sub_x;
			cst = cst_const / sigma;



			//if (pwd != "skip" && iter_all > 0) {
			//	ofstream model(string("model_save/") + pwd + string("_") + graph_name + string("_") + to_string(iter_all) + string(".dat"), std::ios::out | std::ios::binary);
			//	auto i_end = std::chrono::high_resolution_clock::now();
			//	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
			//	double tm = double(eld.count()) * 1.0 * 1e-9;
			//	model << setprecision(searchnum);
			//	model << graph_name << endl;
			//	model << iter_all << endl;
			//	model << n << "\t" << var_num << "\t" << sigma << "\t" << tm + time_record << endl;



			//	for (int i = 0; i < var_num; i++)
			//		model << x[i] << endl;
			//	for (int i = 0; i < var_num; i++)
			//		model << v[i] << endl;
			//	model << u.size() << endl;
			//	for (auto i = u.begin(); i != u.end(); i++) {
			//		model << i->first << " " << i->second << endl;
			//	}
			//}

			std::cout << "start loop " << iter_all + 1 << ",sigma=" << sigma << endl;


			bool flaga = true;
			int iter_1;
			bool rhsflag = false;
			bool lxflag = false;
			double old_lx = 0;
			ori_x = inv_w.cwiseProduct(sub_x);

#pragma omp parallel for
			for (int i = 0; i < generator.data.size(); i++) {

				generator.data[i].y = (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
					+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a + generator.data[i].u / sigma;
			}


			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {

				ori_x = inv_w.cwiseProduct(sub_x);
				//				double bias = 0;
				//#pragma omp parallel  for reduction(+:bias)
				//				for (int i = 0; i < generator.data.size(); i++) {
				//
				//					bias+=fabs(generator.data[i].y - (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
				//						+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a- generator.data[i].u / sigma);
				//				}
				//				cout << "bias:" << bias << endl;
								//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
								//std::cout << "x:" << x.norm() << endl;

				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double, Eigen::RowMajor>ata(sub_var, sub_var);

				{
#ifdef BDEG
					concurrency::concurrent_vector<int*>mtriv;
#else

					concurrency::concurrent_vector<Eigen::Triplet<double>>mtriv;
#endif
					diag = Eigen::VectorXd::Zero(sub_var);
					if (!rhsflag)
					{


						rhs = Eigen::VectorXd::Zero(sub_var);

						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mdiag[idx[1]] += 1;
									mrhs(idx[2]) -= i->y / scale_a;
									mdiag[idx[2]] += 1;
									mrhs(idx[0]) += i->y / scale_a;
									mdiag[idx[0]] += 1;
#ifdef BDEG
									mtriv.push_back(idx);
#else
									//fstream tst("t123.txt", std::ios::app);


									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
									//tst << idx[0] << "\t" << idx[1]<<"\t"<<idx[2] << endl;
									//tst.close();
#endif
								}
							}
						};
						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var)), vdiag(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]), std::ref(vdiag[i])

							));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							rhs += vrhs[i];
							diag += vdiag[i];
						}
						rhs = sigma * rhs.cwiseProduct(inv_w);
						rhs += sigma * pi_c(sub_x + sub_v / sigma, sigma) + cst * (sub_x - x0);




					}
					else {
						vector<Eigen::VectorXd>vdiag(mthread, Eigen::VectorXd::Zero(sub_var));


						auto calc_diag = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								if (i->y > 0) {


									auto idx = i->idx;
									mdiag[idx[1]] += 1;
									mdiag[idx[2]] += 1;
									mdiag[idx[0]] += 1;


									mtriv.push_back(idx);
								}
							}
						};
						vector<thread>vt;
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calc_diag, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vdiag[i])));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							diag += vdiag[i];
						}

					}

					//
#ifdef BDEG
					buildeigen(mtriv, ata, diag);
#else
					ata.setFromTriplets(mtriv.begin(), mtriv.end());
#endif
					auto dt = ata.valuePtr();
					auto nnz = ata.innerNonZeroPtr();
					auto inidx = ata.innerIndexPtr();
					auto oidx = ata.outerIndexPtr();

					//cout << ata.norm() << " " << diag.norm() << endl;
#ifdef DIAGMULT

					diag /= (scale_a * scale_a);
					diag *= sigma;
					diag = diag.cwiseProduct(inv_w).cwiseProduct(inv_w);
#else

					ata += diag.asDiagonal();
					ata /= (scale_a * scale_a);
					ata = inv_w.asDiagonal() * ata * inv_w.asDiagonal();
					ata *= sigma;
#endif


				}


				time1.toc();
				//for (auto i : aj1) {
				//	int ci, cj, ck;
				//	ck = i % n;
				//	cj = i / n;
				//	ci = cj / n;
				//	cj = cj % n;

				//	int idx[3];
				//	idx[0] = ij(ci, cj, n);
				//	idx[1] = ij(ci, ck, n);
				//	idx[2] = ij(cj, ck, n);
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[1], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[2], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[2], sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[1], sigma));
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_v[i].row()<<" "<<test_v[i].col()<<" "<<test_v[i].value() << endl;
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_vb[i].row() << " " << test_vb[i].col() << " " << test_vb[i].value() << endl;
				//}
				//return;





				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;
				//cout << "rhsnorm1:" << rhs.norm() << endl;

				//cout << "rhsnorm2:" << rhs.norm() << endl;
				prt_res << iter_1 << "  " << rhs.norm() << "  ";
				double rhsnm = rhs.norm();

				if (iter_1 > 0 && rhs.norm() < rlim) {
					cout << "rhsout:" << rhs.norm() << endl;
					break;
				}
				//std::cout << "rhs:" << rhs.norm() << endl;
				Eigen::VectorXd direction;
				//if (flaga) {
				//	direction = get_direction(sigma, rhs, &aj1, &diag);
				//}
				//else {
				//	direction= get_direction(sigma, rhs, &aj1);
				//}
				int cgnum;
				time2.tic();
				int kcounter = 0;
				{
					double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);

					Eigen::VectorXd wx = sub_x + sub_v / sigma;

					Eigen::VectorXd r = pi_c(wx, sigma);
					r = Eigen::VectorXd::Ones(r.size()) - r.cwiseEqual(0).cast<double>();

					//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					//double ita = tol;
					//cg.setTolerance(ita / rhs.norm());
#ifdef DIAGMULT
					//ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
					diag += (a + cst) * Eigen::VectorXd::Ones(rhs.size());
#else
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
#endif
					if (wx.cwiseAbs().sum() < 1 / sigma) {
						//ata += sigma * sub_w.cwiseProduct(sub_w).asDiagonal();

						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						//cg.compute(ata);

						//direction = -cg.solve(rhs);
						//kcounter = cg.iterations();
#ifdef DIAGMULT
						diag += sigma * Eigen::VectorXd::Ones(rhs.size());
						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = diag;



#else
						ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = ata.diagonal();
#endif

						Eigen::VectorXd rnew = Eigen::VectorXd::Zero(newdiag.size());


						newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);



#ifdef DIAGMULT
						kcounter = psqmr_inf_diag(ata, inv_w, diag, newdiag, rnew, rhs, direction, sigma, tol_cg);
#else
						kcounter = psqmr_inf(ata, newdiag, rnew, rhs, direction, tol_cg);
#endif






					}
					else {
						if (true) {
#ifdef DIAGMULT
							diag += sigma * r;
							//ata += sigma * r.asDiagonal();
							//cout << "ata:" << ata.diagonal().norm() << "\tdiag:" << diag.norm() << endl;
#else
							ata += sigma * r.asDiagonal();
#endif


							double rnum = r.sum();
							Eigen::VectorXd rnew = r.cwiseProduct(wx.cwiseSign()) * sqrt(sigma / rnum);

							timer diana;
							diana.tic();
							timer ava;

#ifdef SMCG

							cg.compute(ata);
							ava.tic();
							direction = -cg.solve(rhs);
							ava.toc();
							kcounter = cg.iterations();
							ita = std::fmin(0.5, 0.000001 * rnew.norm());
							cg.setTolerance(ita / rnew.norm());
							ava.tic();
							Eigen::VectorXd ar = cg.solve(rnew);
							ava.toc();
							direction += (ar * rnew.dot(direction)) / (1 - rnew.dot(ar));
							kcounter += cg.iterations();
#else
#ifdef EIGENCG
							MatrixReplacement_inf AX(&ata, sub_var, &rnew, sigma);
							Eigen::ConjugateGradient<MatrixReplacement_inf, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg2;
							cg2.compute(AX);
							ita = std::fmin(0.5, 0.000001 * rhs.norm());
							cg2.setTolerance(ita / rhs.norm());
							direction = -cg2.solve(rhs);
							kcounter = cg2.iterations();
#else


							//diag = sigma * diag + (a + cst) * Eigen::VectorXd::Ones(diag.size());
							//diag += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r);
							//diag = Eigen::VectorXd::Ones(diag.size()).cwiseQuotient(diag);





#ifdef DIAGMULT
							Eigen::VectorXd newdiag = diag;
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);
							kcounter = psqmr_inf_diag(ata, inv_w, diag, newdiag, rnew, rhs, direction, sigma, tol_cg);
#else
							Eigen::VectorXd newdiag = ata.diagonal();
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);
							kcounter = psqmr_inf(ata, newdiag, rnew, rhs, direction, tol_cg);
#endif

							//cout << "cg:" << (ata * direction - rnew.dot(direction) * rnew + rhs).norm() << endl;
							//kcounter=psqmr(ata, rnew, rhs, direction);

#endif
#endif
							diana.toc();
							//cout << "cg:" << diana.blink * 100 << " " << kcounter << " " << diana.blink / kcounter * 10000 << endl;
							//cout << ava.time() * 100 << endl;


							//cout << "di:" << direction.norm() << endl;
							//return 0;
							//cout << "woodbury:" << direction.norm() << endl;
						}
						else {



						}
					}



					//std::cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

					//Eigen::PardisoLDLT<SparseMatrix<double> > cg;



					//cgnum = cg.iterations();
					//std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
					//prt_res << cg.iterations() << "\t";
					//std::cout << "pro_res:" << (A * direction+rhs).norm() << endl;
					//std::cout << "CG:" << cg.iterations() << std::endl;
					//std::cout << "di:" << direction.norm() << endl;

				}
				time2.toc();
				//cout << "di:" << setprecision(10) << direction.norm() << endl;
				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;
				Eigen::VectorXd wdirection = direction.cwiseProduct(inv_w);
				time3.tic();
				string tst(":");
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;

					double galp = 0;
					double g0 = direction.dot(rhs);
					if (g0 >= 0) {
						std::cout << "Need a descent direction\n";
						exit(-1);
					}
					//if (mthread == 0) {
					//	axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
					//	(*axbu)();
					//}
					//else {
					//	m_calc_phi(mthread, lx0, n, *axbu, 0, &direction);
					//}


					auto square = [](double x) {return x * x; };



					Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(sub_x.size());




					{


						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								i->y += (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mrhs(idx[2]) -= i->y / scale_a;
									mrhs(idx[0]) += i->y / scale_a;
								}
							}
						};

						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]))
							);
						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							nrhs += vrhs[i];
						}
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);



					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs.cwiseProduct(inv_w);
					nrhs += sigma * pi_c(sub_x + direction + sub_v / sigma, sigma) + cst * (sub_x + direction - x0);
					rhs = nrhs;
					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";

						rhsflag = true;
						lxflag = false;
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					Eigen::VectorXd lxnorm = Eigen::VectorXd::Zero(searchnum);
					//omp_set_num_threads(mthread);
					double lx0 = 0;
					if (lxflag == false) {
#pragma omp parallel for reduction(+:lx0)
						for (int i = 0; i < generator.data.size(); i++) {
							auto idx = generator.data[i].idx;
							float y = generator.data[i].y - (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;

							if (y > 0) {
								lx0 += square(y);
							}
						}
						lx0 = sigma / 2 * lx0 + sigma / 2 * pi_c(sub_x + sub_v / sigma, sigma).squaredNorm();

						lx0 += func_p(sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma)) + cst / 2 * (sub_x - x0).squaredNorm();

					}
					else {





						//	cout << "lxbias:" << lx0 - old_lx << endl;

						lx0 = old_lx;
					}










					//lx0 = lxnorm[0];




					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();
					alpha[0] = 1;

					for (int i = 1; i < searchnum; i++) {

						double lx = 0;

#pragma omp parallel for reduction(+:lx)
						for (int ii = 0; ii < generator.data.size(); ii++) {
							auto idx = generator.data[ii].idx;
							generator.data[ii].y -= (alpha[i - 1] - alpha[i]) * (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;

							if (generator.data[ii].y > 0) {
								lx += square(generator.data[ii].y);
							}

						}


						double alp = alpha[i];



						lx = sigma / 2 * lx;
						lx += sigma / 2 * pi_c(sub_x + alp * direction + sub_v / sigma, sigma).squaredNorm();
						lx += func_p(sub_x + alp * direction + sub_v / sigma - pi_c(sub_x + alp * direction + sub_v / sigma, sigma)) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();
						lxflag = true;
						old_lx = lx;
						if (lx <= down) {

							if (i == 1) {
								rhsflag = true;

							}
							else {
								rhsflag = false;
							}
							return alp;
						}

					}


					rhsflag = false;
					return alpha[searchnum - 1];





				}();
				time3.toc();
				//double alpha = 1;



				sub_x += dalpha * direction;
				//cout << "x:" << setprecision(10) << x.norm() << endl;
				//cout << "x:" << sub_x.norm() << endl;
				prt_res << dalpha << std::endl;
				std::cout << setprecision(10) << iter_1 << "\t" << rhsnm << "\t" << kcounter << "\t" << dalpha << endl;
			}


			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				consflag = true;

#pragma omp parallel for reduction(+:kkt_y,kkt_u,y_norm)
				for (int ii = 0; ii < generator.data.size(); ii++) {
					auto idx = generator.data[ii].idx;
					double y = generator.data[ii].y;
					float u = generator.data[ii].u;
					if (y < 0) {
						kkt_u += (u / sigma) * (u / sigma);
					}
					else {
						kkt_u += (y - u / sigma) * (y - u / sigma);
					}
					double ny = y < 0 ? y : 0;
					double nu = sigma * (y - ny);
					if (y - u / sigma > 0.01)consflag = false;

					//double y = axbu < 0 ? axbu : 0;
					y_norm += ny * ny;
					//	cout << "wtf\n";
						//*kkt_u += (axb - y) * (axb - y);


					double y_plus_u = (ny + nu > 0 ? ny : -nu);
					kkt_y += y_plus_u * y_plus_u;

					generator.data[ii].u = nu;

				}


				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();

				total_kkt_u = kkt_u;
				total_y_norm = y_norm;


				//std::cout << "u:" << u.size() << endl;
				Eigen::VectorXd z = sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma);
				sub_v = sigma * pi_c(sub_x + sub_v / sigma, sigma);

				y_norm = sqrt(y_norm);
				kkt_u = sqrt(kkt_u) / (1 + y_norm);
				kkt_y = sqrt(kkt_y) / (1 + y_norm);

				std::cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
				Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(sub_x.size());


				Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_x.size());



				//m_calc_atu(mthread, &atu, n, au);

				{
					int sz = generator.data.size();
					int pnum = mthread;


					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(sub_var);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i]), scale_a));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu.cwiseProduct(inv_w) + sub_v).norm() / (1 + scale_a * scale_a * sub_v.norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_x - z).norm() / scale_a / (1 + sub_x.norm() / scale_a);
				kkt_z = ((pi_c(z + sub_v, 1.0) - sub_v).norm() / scale_a) / (1 + z.norm() / scale_a);
				//for (int i = 0; i < var_num; i++)
				//    t_res += abs(x[i] * w[i]);
				//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
				double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
				std::cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
				std::cout << "kkt:" << kkt_final << endl;

				rp = std::fmax(kkt_u, kkt_v);
				rd = kkt_x;
				kkt = kkt_final;

				min_kkt = min(kkt, min_kkt);
				tol_cg = min(min_kkt, 5.0e-3);
				rlim = max(1e-4, min(kkt, 1.0) * 0.9);
				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;


#pragma omp parallel for reduction(+:dobj)
				for (int i = 0; i < generator.data.size(); i++) {
					int* idx = generator.data[i].idx;
					dobj += generator.data[i].u * (sub_b[idx[0]] - sub_b[idx[1]] - sub_b[idx[2]]);
				}
				dobj /= scale_a * scale_a;




			}
			time4.toc();

			double pobj = func_p(sub_x) / scale_a;
			cout << "pobj:" << pobj << ",dobj:" << dobj << endl;


			prt_res << "==================\n";
			if (kkt < tol && consflag && fabs((pobj - dobj) / (1 + dobj)) < 1e-2) {

				sigini = max(sigini / 10, 1e-4);

				cst_const = sigini * 1e-4;

				//int pnum = mthread;


				//vector<thread>tv;

				//std::vector<almfunc_rhs>vf;
				//for (int i = 0; i < mthread; i++) {
				//	vf.push_back(almfunc_rhs(0, NULL, mrhs + i));
				//}
				//int sz = generator.data.size();
				//for (int i = 0; i < pnum; i++) {
				//	tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
				//		n, std::ref(sub_x), std::ref(sub_b), sigma));

				//}
				//for (int i = 0; i < pnum; i++)
				//	tv[i].join();
				//for (int i = 0; i < pnum; i++)
				//	nrhs += mrhs[i];
				//delete[]mrhs;
				break;
			}

			/*
			if (iter_1 == sub_iter) {
				sigma = 1;
				primal_win = 0;
				dual_win = 0;
			}
			else {
				if (5 * rp < rd)
					primal_win += 1;
				else if (rd < rp)
					dual_win += 1;


				if (iter_all % 3 == 2) {

					if (primal_win > std::fmax(1, 1.2 * dual_win)) {
						primal_win = 0;
						sigma = min(sigma / 1.2, 1.0e7);
					}
					else if ((dual_win > std::fmax(1, 1.2 * primal_win)) && (kkt > kkt0 * 0.7)) {

						dual_win = 0;
						sigma = max(sigma * 1.5, 1.0e0);
					}

				}
				else if (kkt > kkt0 * 0.7) {
					sigma = std::fmin(sigma * 1.5, 1.0e7);
				}
			}
			*/
			if (iter_1 == sub_iter && iter_all > 6) {


				sigma = max(sigma / 10, 1.0e-4);
			}

			else if (iter_1 != sub_iter) {
				if (kkt > 4.0e-4)
					sigma = min(sigma * 1.5, 1.0e7);
				else
					sigma = min(sigma * 1.2, 1.0e7);
			}






			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
			kkt0 = kkt;

		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
		itervec.push_back(iter_all);
		for (auto& i:colidx) {
			


				sparse_x[i.first] = sub_x[i.second] / sub_w[i.second];
				//b[i] = sub_b[colidx[i]];
				//w[i] = sub_w[colidx[i]];
				sparse_v[i.first] = sub_v[i.second];
			
		}
	}


	//for (auto iter = u.begin(); iter != u.end(); iter++) {
	//	int64_t longindex = iter->first;
	//	int i, j, k;
	//	k = longindex % n;
	//	j = longindex / n;
	//	i = j / n;
	//	j = j % n;
	//	int idx[3];
	//	idx[0] = ij(i, j, n);
	//	idx[1] = ij(i, k, n);
	//	idx[2] = ij(j, k, n);
	//	dobj += iter->second * (b[idx[0]]-b[idx[1]]-b[idx[2]]);


	//}


	std::cout << "func_value:" << std::setprecision(8) << 0 << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	cout << "diana:" << time5.time() << endl;
	
	//x = x.cwiseQuotient(w);
	return -1;




}

double norm_inf::run_axbu(double cst_const, int con_const, int sub_iter) {
	scale_a = 6;
	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();
	b = b * scale_a;

	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 0.01;

	double c_t = 0;
	//std::unordered_map<int, int>colidx;

	timer time1, time2, time3, time4, time5;

	cout << "m:" << mthread << endl;
	int mxtnum = thread::hardware_concurrency();
	cons_generator generator(n, mxtnum, cgflg);
	generator.consvio_limit = scale_a * 0.01;
	double total_kkt_u = 0, total_y_norm = 0;

	double sigini = 1e-2;
	runhist[string("iter")] = -1;

	bool consflag = false;
	while (true) {
		runhist[string("iter")]+=1;
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(x, b, scale_a * tol);
		}
		sort(generator.data.begin(), generator.data.end(), [](cons_full &a, cons_full &b) {
			if (min(a.idx[0], min(a.idx[1], a.idx[2])) != min(b.idx[0], min(b.idx[1], b.idx[2]))) {
				return min(a.idx[0], min(a.idx[1], a.idx[2])) != min(b.idx[0], min(b.idx[1], b.idx[2]));
			}
			return max(a.idx[0], max(a.idx[1], a.idx[2])) != max(b.idx[0], max(b.idx[1], b.idx[2]));

			});
		double totalkkt = sqrt(generator.k1 + total_kkt_u) / (1 + sqrt(generator.k2 + total_y_norm));
		cout << "totalkkt:" << totalkkt << endl;
		if (generator.consvio)
			break;


		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = sigini;
		cst_const = sigma * 1e-3;
		cout << "cst:" << cst_const << endl;
		int sub_var = colidx_num;
		std::cout << "total num of variables:" << var_num << " num of current variables:" << sub_var << endl;

		auto c_end = std::chrono::high_resolution_clock::now();
		auto c_eld = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start);
		c_t += double(c_eld.count()) * 1.0 * 1e-9;

		double dobj = 0;



		double	primal_win = 0;
		double dual_win = 0;
		double kkt0 = 1;
		kkt = 0.1;

		Eigen::VectorXd sub_x = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_b = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_w = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_v = Eigen::VectorXd::Zero(colidx_num);



		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				sub_x[colidx[i]] = x[i] * w[i];
				sub_b[colidx[i]] = b[i];
				sub_w[colidx[i]] = w[i];
				sub_v[colidx[i]] = v[i];
			}
		}
		Eigen::VectorXd inv_w = Eigen::VectorXd::Ones(sub_w.size()).cwiseQuotient(sub_w);
		Eigen::VectorXd ori_x = inv_w.cwiseProduct(sub_x);
		double kkt = 0.1;
		cout << sub_var << endl;
		double tol_cg = 5e-2;

		double min_kkt = 0.1;
		rlim = 1e-1;
		for (iter_all = 0; iter_all < 1000; iter_all++) {


			x0 = sub_x;
			cst = cst_const / sigma;
			


			//if (pwd != "skip" && iter_all > 0) {
			//	ofstream model(string("model_save/") + pwd + string("_") + graph_name + string("_") + to_string(iter_all) + string(".dat"), std::ios::out | std::ios::binary);
			//	auto i_end = std::chrono::high_resolution_clock::now();
			//	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
			//	double tm = double(eld.count()) * 1.0 * 1e-9;
			//	model << setprecision(searchnum);
			//	model << graph_name << endl;
			//	model << iter_all << endl;
			//	model << n << "\t" << var_num << "\t" << sigma << "\t" << tm + time_record << endl;



			//	for (int i = 0; i < var_num; i++)
			//		model << x[i] << endl;
			//	for (int i = 0; i < var_num; i++)
			//		model << v[i] << endl;
			//	model << u.size() << endl;
			//	for (auto i = u.begin(); i != u.end(); i++) {
			//		model << i->first << " " << i->second << endl;
			//	}
			//}

			std::cout << "start loop " << iter_all + 1 << ",sigma=" << sigma << endl;


			bool flaga = true;
			int iter_1;
			bool rhsflag = false;
			bool lxflag = false;
			double old_lx = 0;
			ori_x = inv_w.cwiseProduct(sub_x);

#pragma omp parallel for
			for (int i = 0; i < generator.data.size(); i++) {

				generator.data[i].y = (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
					+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a + generator.data[i].u / sigma;
			}


			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {
			
				ori_x = inv_w.cwiseProduct(sub_x);
//				double bias = 0;
//#pragma omp parallel  for reduction(+:bias)
//				for (int i = 0; i < generator.data.size(); i++) {
//
//					bias+=fabs(generator.data[i].y - (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
//						+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a- generator.data[i].u / sigma);
//				}
//				cout << "bias:" << bias << endl;
				//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
				//std::cout << "x:" << x.norm() << endl;
				
				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double,Eigen::RowMajor>ata(sub_var, sub_var);
							
				{
#ifdef BDEG
					concurrency::concurrent_vector<int*>mtriv;
#else

					concurrency::concurrent_vector<Eigen::Triplet<double>>mtriv;
#endif
					diag = Eigen::VectorXd::Zero(sub_var);	
					if (!rhsflag)
					{


						rhs = Eigen::VectorXd::Zero(sub_var);
						
						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end,Eigen::VectorXd&mrhs,Eigen::VectorXd&mdiag) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mdiag[idx[1]] += 1;
									mrhs(idx[2]) -= i->y / scale_a;
									mdiag[idx[2]] += 1;
									mrhs(idx[0]) += i->y / scale_a;
									mdiag[idx[0]] += 1;
#ifdef BDEG
									mtriv.push_back(idx);
#else
									//fstream tst("t123.txt", std::ios::app);


									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
									//tst << idx[0] << "\t" << idx[1]<<"\t"<<idx[2] << endl;
									//tst.close();
#endif
								}
							}
						};
						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var)), vdiag(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]), std::ref(vdiag[i])

							));
								
						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							rhs += vrhs[i];
							diag += vdiag[i];
						}
						rhs = sigma * rhs.cwiseProduct(inv_w);
						rhs+=sigma * pi_c(sub_x + sub_v / sigma, sigma) + cst * (sub_x - x0);




					}
					else {
						vector<Eigen::VectorXd>vdiag(mthread, Eigen::VectorXd::Zero(sub_var));


						auto calc_diag = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								if (i->y > 0) {


									auto idx = i->idx;
									mdiag[idx[1]] += 1;
									mdiag[idx[2]] += 1;
									mdiag[idx[0]] += 1;


									mtriv.push_back(idx);
								}
							}
						};
						vector<thread>vt;
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calc_diag, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vdiag[i])));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							diag += vdiag[i];
						}

					}
					
					//
#ifdef BDEG
					buildeigen(mtriv, ata, diag);
#else
					ata.setFromTriplets(mtriv.begin(), mtriv.end());
#endif
					auto dt = ata.valuePtr();
					auto nnz = ata.innerNonZeroPtr();
					auto inidx = ata.innerIndexPtr();
					auto oidx = ata.outerIndexPtr();

					//cout << ata.norm() << " " << diag.norm() << endl;
#ifdef DIAGMULT

					diag /= (scale_a * scale_a);
					diag *= sigma;
					diag = diag.cwiseProduct(inv_w).cwiseProduct(inv_w);
#else

					ata += diag.asDiagonal();
					ata/= (scale_a * scale_a);
					ata = inv_w.asDiagonal() * ata * inv_w.asDiagonal();
					ata *= sigma;
#endif


				}


				time1.toc();
				//for (auto i : aj1) {
				//	int ci, cj, ck;
				//	ck = i % n;
				//	cj = i / n;
				//	ci = cj / n;
				//	cj = cj % n;

				//	int idx[3];
				//	idx[0] = ij(ci, cj, n);
				//	idx[1] = ij(ci, ck, n);
				//	idx[2] = ij(cj, ck, n);
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[1], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[2], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[2], sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[1], sigma));
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_v[i].row()<<" "<<test_v[i].col()<<" "<<test_v[i].value() << endl;
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_vb[i].row() << " " << test_vb[i].col() << " " << test_vb[i].value() << endl;
				//}
				//return;



				

				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;
				//cout << "rhsnorm1:" << rhs.norm() << endl;
				
				//cout << "rhsnorm2:" << rhs.norm() << endl;
				prt_res << iter_1 << "  " << rhs.norm() << "  ";
				double rhsnm = rhs.norm();

				if (iter_1>0&&rhs.norm() < rlim) {
					cout << "rhsout:" << rhs.norm() << endl;
					break;
				}
				//std::cout << "rhs:" << rhs.norm() << endl;
				Eigen::VectorXd direction;
				//if (flaga) {
				//	direction = get_direction(sigma, rhs, &aj1, &diag);
				//}
				//else {
				//	direction= get_direction(sigma, rhs, &aj1);
				//}
				int cgnum;
				time2.tic();
				int kcounter = 0;
				{
					double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);

					Eigen::VectorXd wx = sub_x + sub_v / sigma;

					Eigen::VectorXd r = pi_c(wx, sigma);
					r = Eigen::VectorXd::Ones(r.size()) - r.cwiseEqual(0).cast<double>();

					//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					//double ita = tol;
					//cg.setTolerance(ita / rhs.norm());
#ifdef DIAGMULT
					//ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
					diag += (a + cst) * Eigen::VectorXd::Ones(rhs.size());
#else
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
#endif
					if (wx.cwiseAbs().sum() < 1 / sigma) {
						//ata += sigma * sub_w.cwiseProduct(sub_w).asDiagonal();

						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						//cg.compute(ata);

						//direction = -cg.solve(rhs);
						//kcounter = cg.iterations();
#ifdef DIAGMULT
						diag += sigma * Eigen::VectorXd::Ones(rhs.size());
						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = diag;



#else
						ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = ata.diagonal();
#endif
						
						Eigen::VectorXd rnew = Eigen::VectorXd::Zero(newdiag.size());

						
						newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);



#ifdef DIAGMULT
						kcounter = psqmr_inf_diag(ata,inv_w,diag, newdiag, rnew,rhs, direction,sigma, tol_cg);
#else
						kcounter = psqmr_inf(ata, newdiag, rnew, rhs, direction, tol_cg);
#endif

						




					}
					else {
						if (true) {
#ifdef DIAGMULT
							diag += sigma * r;
							//ata += sigma * r.asDiagonal();
							//cout << "ata:" << ata.diagonal().norm() << "\tdiag:" << diag.norm() << endl;
#else
							ata += sigma * r.asDiagonal();
#endif

							
							double rnum = r.sum();
							Eigen::VectorXd rnew = r.cwiseProduct(wx.cwiseSign()) * sqrt(sigma / rnum);

							timer diana;
							diana.tic();
							timer ava;

#ifdef SMCG

							cg.compute(ata);
							ava.tic();
							direction = -cg.solve(rhs);
							ava.toc();
							kcounter = cg.iterations();
							ita = std::fmin(0.5, 0.000001 * rnew.norm());
							cg.setTolerance(ita / rnew.norm());
							ava.tic();
							Eigen::VectorXd ar = cg.solve(rnew);
							ava.toc();
							direction += (ar * rnew.dot(direction)) / (1 - rnew.dot(ar));
							kcounter += cg.iterations();
#else
#ifdef EIGENCG
							MatrixReplacement_inf AX(&ata, sub_var, &rnew, sigma);
							Eigen::ConjugateGradient<MatrixReplacement_inf, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg2;
							cg2.compute(AX);
							ita = std::fmin(0.5, 0.000001 * rhs.norm());
							cg2.setTolerance(ita / rhs.norm());
							direction = -cg2.solve(rhs);
							kcounter = cg2.iterations();
#else


							//diag = sigma * diag + (a + cst) * Eigen::VectorXd::Ones(diag.size());
							//diag += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r);
							//diag = Eigen::VectorXd::Ones(diag.size()).cwiseQuotient(diag);





#ifdef DIAGMULT
							Eigen::VectorXd newdiag =diag;
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);
							kcounter = psqmr_inf_diag(ata, inv_w, diag, newdiag, rnew, rhs, direction, sigma, tol_cg);
#else
							Eigen::VectorXd newdiag = ata.diagonal();
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);
							kcounter = psqmr_inf(ata, newdiag, rnew, rhs, direction, tol_cg);
#endif

							//cout << "cg:" << (ata * direction - rnew.dot(direction) * rnew + rhs).norm() << endl;
							//kcounter=psqmr(ata, rnew, rhs, direction);

#endif
#endif
							diana.toc();
							//cout << "cg:" << diana.blink * 100 << " " << kcounter << " " << diana.blink / kcounter * 10000 << endl;
							//cout << ava.time() * 100 << endl;


							//cout << "di:" << direction.norm() << endl;
							//return 0;
							//cout << "woodbury:" << direction.norm() << endl;
						}
						else {



						}
					}



					//std::cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

					//Eigen::PardisoLDLT<SparseMatrix<double> > cg;



					//cgnum = cg.iterations();
					//std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
					//prt_res << cg.iterations() << "\t";
					//std::cout << "pro_res:" << (A * direction+rhs).norm() << endl;
					//std::cout << "CG:" << cg.iterations() << std::endl;
					//std::cout << "di:" << direction.norm() << endl;

				}
				time2.toc();
				//cout << "di:" << setprecision(10) << direction.norm() << endl;
				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;
				Eigen::VectorXd wdirection = direction.cwiseProduct(inv_w);
				time3.tic();
				string tst(":");
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;
					
					double galp = 0;
					double g0 = direction.dot(rhs);
					if (g0 >= 0) {
						std::cout << "Need a descent direction\n";
						exit(-1);
					}
					//if (mthread == 0) {
					//	axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
					//	(*axbu)();
					//}
					//else {
					//	m_calc_phi(mthread, lx0, n, *axbu, 0, &direction);
					//}


					auto square = [](double x) {return x * x; };



					Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(sub_x.size());


					

					{
			

						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								i->y += (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]])/scale_a;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mrhs(idx[2]) -= i->y / scale_a;
									mrhs(idx[0]) += i->y / scale_a;
								}
							}
						};

						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]))
							);
						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							nrhs += vrhs[i];
						}
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);
					


					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs.cwiseProduct(inv_w);
					nrhs += sigma * pi_c(sub_x + direction + sub_v / sigma, sigma) + cst * (sub_x + direction - x0);
					rhs = nrhs;
					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";
						
						rhsflag = true;
						lxflag = false;
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					Eigen::VectorXd lxnorm = Eigen::VectorXd::Zero(searchnum);
					//omp_set_num_threads(mthread);
					double lx0 = 0;
					if (lxflag==false) {
#pragma omp parallel for reduction(+:lx0)
						for (int i = 0; i < generator.data.size(); i++) {
							auto idx = generator.data[i].idx;
							float y = generator.data[i].y - (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;
							
							if (y > 0) {
								lx0 += square(y);
							}
						}
						lx0 = sigma / 2 * lx0 + sigma / 2 * pi_c(sub_x + sub_v / sigma, sigma).squaredNorm();

						lx0 += func_p(sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma)) + cst / 2 * (sub_x - x0).squaredNorm();

					}
					else {




						
					//	cout << "lxbias:" << lx0 - old_lx << endl;

						lx0 = old_lx;
					}
					









					//lx0 = lxnorm[0];




					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();
					alpha[0] = 1;

					for (int i = 1; i < searchnum; i++) {

						double lx = 0;

#pragma omp parallel for reduction(+:lx)
						for (int ii = 0; ii < generator.data.size(); ii++) {
							auto idx = generator.data[ii].idx;
							generator.data[ii].y-=(alpha[i-1]-alpha[i])* (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;

								if (generator.data[ii].y > 0) {
									lx += square(generator.data[ii].y);
								}
							
						}


						double alp = alpha[i];

						

						lx = sigma / 2 * lx;
						lx += sigma / 2 * pi_c(sub_x + alp * direction + sub_v / sigma, sigma).squaredNorm();
						lx += func_p(sub_x + alp * direction + sub_v / sigma - pi_c(sub_x + alp * direction + sub_v / sigma, sigma)) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();
						lxflag = true;
						old_lx = lx;
						if (lx <= down) {
							
							if (i == 1) {
								rhsflag = true;
								
							}
							else {
								rhsflag = false;
							}
							return alp;
						}

					}


					rhsflag = false;
					return alpha[searchnum - 1];





				}();
				time3.toc();
				//double alpha = 1;



				sub_x += dalpha * direction;
				//cout << "x:" << setprecision(10) << x.norm() << endl;
				//cout << "x:" << sub_x.norm() << endl;
				prt_res << dalpha << std::endl;
				std::cout << setprecision(10) << iter_1 << "\t" << rhsnm << "\t" << kcounter << "\t" << dalpha << endl;
			}


			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				consflag = true;

#pragma omp parallel for reduction(+:kkt_y,kkt_u,y_norm)
				for (int ii = 0; ii < generator.data.size(); ii++) {
					auto idx = generator.data[ii].idx;
					double y = generator.data[ii].y;
					float u = generator.data[ii].u;
					if (y < 0) {
						kkt_u += (u / sigma) * (u / sigma);
					}
					else {
						kkt_u += (y - u / sigma) * (y - u / sigma);
					}
					double ny = y < 0 ? y : 0;
					double nu=  sigma*(y-ny) ;
					if (y - u / sigma > 0.01)consflag = false;

					//double y = axbu < 0 ? axbu : 0;
					y_norm += ny * ny;
					//	cout << "wtf\n";
						//*kkt_u += (axb - y) * (axb - y);

					
					double y_plus_u = (ny + nu > 0 ? ny : -nu);
					kkt_y += y_plus_u * y_plus_u;

					generator.data[ii].u =nu ;

				}


				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();

				total_kkt_u = kkt_u;
				total_y_norm = y_norm;


				//std::cout << "u:" << u.size() << endl;
				Eigen::VectorXd z = sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma);
				sub_v = sigma * pi_c(sub_x + sub_v / sigma, sigma);

				y_norm = sqrt(y_norm);
				kkt_u = sqrt(kkt_u) / (1 + y_norm);
				kkt_y = sqrt(kkt_y) / (1 + y_norm);

				std::cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
				Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(sub_x.size());


				Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_x.size());



				//m_calc_atu(mthread, &atu, n, au);

				{
					int sz = generator.data.size();
					int pnum = mthread;

					
					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(sub_var);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i]), scale_a));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu.cwiseProduct(inv_w) + sub_v).norm() / (1 + scale_a*scale_a*sub_v.norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_x - z).norm()/scale_a / (1 + sub_x.norm()/scale_a);
				//kkt_z = ((pi_c(z + sub_v, 1.0) - sub_v).norm() / scale_a) / (1 + z.norm() / scale_a);
				kkt_z = ((z -( sub_x + sub_v- pi_c(sub_x + sub_v, 1.0))).norm() / scale_a) / (1 + z.norm() / scale_a);

				//for (int i = 0; i < var_num; i++)
				//    t_res += abs(x[i] * w[i]);
				//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
				double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
				std::cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
				std::cout << "kkt:" << kkt_final << endl;
				runhist[string("kkt")] = kkt_final;
				rp = std::fmax(kkt_u, kkt_v);
				rd = kkt_x;
				kkt = kkt_final;

				min_kkt = min(kkt, min_kkt);
				tol_cg = min(min_kkt, 5.0e-3);
				rlim = max(1e-4, min(kkt, 1.0) * 0.9);
				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;

				dobj = 0;
#pragma omp parallel for reduction(+:dobj)
				for (int i = 0; i < generator.data.size(); i++) {
					int* idx = generator.data[i].idx;
					dobj += generator.data[i].u * (sub_b[idx[0]] - sub_b[idx[1]] - sub_b[idx[2]]);
				}
				dobj /= scale_a * scale_a;




			}
			time4.toc();

			double pobj = func_p(sub_x) / scale_a;
			cout << "pobj:" << pobj << ",dobj:" << dobj << endl;

			runhist[string("R_g")] = fabs((pobj - dobj) / (1 + dobj));
			prt_res << "==================\n";
			if (kkt < tol&&consflag&&fabs((pobj-dobj)/(1+dobj))<1e-4) {

				sigini = max(sigini / 10, 1e-4);
				
				cst_const = sigini * 1e-4;

				//int pnum = mthread;


				//vector<thread>tv;

				//std::vector<almfunc_rhs>vf;
				//for (int i = 0; i < mthread; i++) {
				//	vf.push_back(almfunc_rhs(0, NULL, mrhs + i));
				//}
				//int sz = generator.data.size();
				//for (int i = 0; i < pnum; i++) {
				//	tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
				//		n, std::ref(sub_x), std::ref(sub_b), sigma));

				//}
				//for (int i = 0; i < pnum; i++)
				//	tv[i].join();
				//for (int i = 0; i < pnum; i++)
				//	nrhs += mrhs[i];
				//delete[]mrhs;
				break;
			}

			/*
			if (iter_1 == sub_iter) {
				sigma = 1;
				primal_win = 0;
				dual_win = 0;
			}
			else {
				if (5 * rp < rd)
					primal_win += 1;
				else if (rd < rp)
					dual_win += 1;


				if (iter_all % 3 == 2) {

					if (primal_win > std::fmax(1, 1.2 * dual_win)) {
						primal_win = 0;
						sigma = min(sigma / 1.2, 1.0e7);
					}
					else if ((dual_win > std::fmax(1, 1.2 * primal_win)) && (kkt > kkt0 * 0.7)) {

						dual_win = 0;
						sigma = max(sigma * 1.5, 1.0e0);
					}

				}
				else if (kkt > kkt0 * 0.7) {
					sigma = std::fmin(sigma * 1.5, 1.0e7);
				}
			}
			*/
			if (iter_1 == sub_iter && iter_all > 6) {


				sigma = max(sigma / 10, 1.0e-4);
			}

			else if(iter_1!=sub_iter) {
				if (kkt > 4.0e-4)
					sigma = min(sigma * 1.5, 1.0e7);
				else
					sigma = min(sigma * 1.2, 1.0e7);
			}

			




			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
			kkt0 = kkt;

		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
		itervec.push_back(iter_all);
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]] / sub_w[colidx[i]];
				//b[i] = sub_b[colidx[i]];
				//w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}


	}


	//for (auto iter = u.begin(); iter != u.end(); iter++) {
	//	int64_t longindex = iter->first;
	//	int i, j, k;
	//	k = longindex % n;
	//	j = longindex / n;
	//	i = j / n;
	//	j = j % n;
	//	int idx[3];
	//	idx[0] = ij(i, j, n);
	//	idx[1] = ij(i, k, n);
	//	idx[2] = ij(j, k, n);
	//	dobj += iter->second * (b[idx[0]]-b[idx[1]]-b[idx[2]]);


	//}


	std::cout << "func_value:" << std::setprecision(8) << func_p(x) <<  endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	cout << "diana:" << time5.time() << endl;
	x /= scale_a;
	b /= scale_a;
	//x = x.cwiseQuotient(w);
		runhist[string("time_of_cg")] = generator.ava.time();


	return func_p(w.cwiseProduct(x));




}
double norm_2::run_axbu(double cst_const, int con_const, int sub_iter) {
	scale_a = 6;
	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();
	b = b * scale_a;
	w = w.cwiseSqrt();
	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 0.01;
	double dobj = 0;
	double c_t = 0;
	//std::unordered_map<int, int>colidx;

	timer time1, time2, time3, time4, time5;

	cout << "m:" << mthread << endl;
	int mxtnum = thread::hardware_concurrency();
	cons_generator generator(n, mxtnum, cgflg);
	generator.consvio_limit = scale_a * 0.01;
	double total_kkt_u = 0, total_y_norm = 0;

	double sigini = 1e-2;
	runhist[string("iter")] = -1;
	bool consflag = false;
	while (true) {
		runhist[string("iter")] += 1;
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(x, b, scale_a * tol);
		}
		sort(generator.data.begin(), generator.data.end(), [](cons_full& a, cons_full& b) {
			if (min(a.idx[0], min(a.idx[1], a.idx[2])) != min(b.idx[0], min(b.idx[1], b.idx[2]))) {
				return min(a.idx[0], min(a.idx[1], a.idx[2])) != min(b.idx[0], min(b.idx[1], b.idx[2]));
			}
			return max(a.idx[0], max(a.idx[1], a.idx[2])) != max(b.idx[0], max(b.idx[1], b.idx[2]));

			});
		double totalkkt = sqrt(generator.k1 + total_kkt_u) / (1 + sqrt(generator.k2 + total_y_norm));
		cout << "totalkkt:" << totalkkt << endl;
		if (generator.consvio)
			break;


		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = sigini;
		cst_const = 0;
		cout << "cst:" << cst_const << endl;
		int sub_var = colidx_num;
		std::cout << "total num of variables:" << var_num << " num of current variables:" << sub_var << endl;

		auto c_end = std::chrono::high_resolution_clock::now();
		auto c_eld = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start);
		c_t += double(c_eld.count()) * 1.0 * 1e-9;





		double	primal_win = 0;
		double dual_win = 0;
		double kkt0 = 1;
		kkt = 0.1;

		Eigen::VectorXd sub_x = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_b = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_w = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_v = Eigen::VectorXd::Zero(colidx_num);



		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				sub_x[colidx[i]] = x[i] * w[i];
				sub_b[colidx[i]] = b[i];
				sub_w[colidx[i]] = w[i];
				sub_v[colidx[i]] = v[i];
			}
		}
		Eigen::VectorXd inv_w = Eigen::VectorXd::Ones(sub_w.size()).cwiseQuotient(sub_w);
		Eigen::VectorXd ori_x = inv_w.cwiseProduct(sub_x);
		double kkt = 0.1;
		cout << sub_var << endl;
		double tol_cg = 5e-2;

		double min_kkt = 0.1;
		rlim = 1e-1;
		for (iter_all = 0; iter_all < 1000; iter_all++) {


			x0 = sub_x;
			cst = cst_const / sigma;



			//if (pwd != "skip" && iter_all > 0) {
			//	ofstream model(string("model_save/") + pwd + string("_") + graph_name + string("_") + to_string(iter_all) + string(".dat"), std::ios::out | std::ios::binary);
			//	auto i_end = std::chrono::high_resolution_clock::now();
			//	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
			//	double tm = double(eld.count()) * 1.0 * 1e-9;
			//	model << setprecision(searchnum);
			//	model << graph_name << endl;
			//	model << iter_all << endl;
			//	model << n << "\t" << var_num << "\t" << sigma << "\t" << tm + time_record << endl;



			//	for (int i = 0; i < var_num; i++)
			//		model << x[i] << endl;
			//	for (int i = 0; i < var_num; i++)
			//		model << v[i] << endl;
			//	model << u.size() << endl;
			//	for (auto i = u.begin(); i != u.end(); i++) {
			//		model << i->first << " " << i->second << endl;
			//	}
			//}

			std::cout << "start loop " << iter_all + 1 << ",sigma=" << sigma << endl;


			bool flaga = true;
			int iter_1;
			bool rhsflag = false;
			bool lxflag = false;
			double old_lx = 0;
			ori_x = inv_w.cwiseProduct(sub_x);

#pragma omp parallel for
			for (int i = 0; i < generator.data.size(); i++) {

				generator.data[i].y = (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
					+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a + generator.data[i].u / sigma;
			}


			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {

				ori_x = inv_w.cwiseProduct(sub_x);
				//				double bias = 0;
				//#pragma omp parallel  for reduction(+:bias)
				//				for (int i = 0; i < generator.data.size(); i++) {
				//
				//					bias+=fabs(generator.data[i].y - (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
				//						+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a- generator.data[i].u / sigma);
				//				}
				//				cout << "bias:" << bias << endl;
								//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
								//std::cout << "x:" << x.norm() << endl;

				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double, Eigen::RowMajor>ata(sub_var, sub_var);

				{
#ifdef BDEG
					concurrency::concurrent_vector<int*>mtriv;
#else

					concurrency::concurrent_vector<Eigen::Triplet<double>>mtriv;
#endif
					diag = Eigen::VectorXd::Zero(sub_var);
					if (!rhsflag)
					{


						rhs = Eigen::VectorXd::Zero(sub_var);

						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mdiag[idx[1]] += 1;
									mrhs(idx[2]) -= i->y / scale_a;
									mdiag[idx[2]] += 1;
									mrhs(idx[0]) += i->y / scale_a;
									mdiag[idx[0]] += 1;
#ifdef BDEG
									mtriv.push_back(idx);
#else
									//fstream tst("t123.txt", std::ios::app);


									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
									//tst << idx[0] << "\t" << idx[1]<<"\t"<<idx[2] << endl;
									//tst.close();
#endif
								}
							}
						};
						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var)), vdiag(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]), std::ref(vdiag[i])

							));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							rhs += vrhs[i];
							diag += vdiag[i];
						}
						rhs = sigma * rhs.cwiseProduct(inv_w);
						rhs += sigma * pi_c(sub_x + sub_v / sigma, sigma) + cst * (sub_x - x0);



					}
					else {

						vector<Eigen::VectorXd>vdiag(mthread, Eigen::VectorXd::Zero(sub_var));


						auto calc_diag = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								if (i->y > 0) {


									auto idx = i->idx;
									mdiag[idx[1]] += 1;
									mdiag[idx[2]] += 1;
									mdiag[idx[0]] += 1;


									mtriv.push_back(idx);
								}
							}
						};
						vector<thread>vt;
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calc_diag, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vdiag[i])));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							diag += vdiag[i];
						}


					}

					//
#ifdef BDEG
					buildeigen(mtriv, ata, diag);
#else
					ata.setFromTriplets(mtriv.begin(), mtriv.end());
#endif
					auto dt = ata.valuePtr();
					auto nnz = ata.innerNonZeroPtr();
					auto inidx = ata.innerIndexPtr();
					auto oidx = ata.outerIndexPtr();

					//cout << ata.norm() << " " << diag.norm() << endl;
#ifdef DIAGMULT

					diag /= (scale_a * scale_a);
					diag *= sigma;
					diag = diag.cwiseProduct(inv_w).cwiseProduct(inv_w);
#else

					ata += diag.asDiagonal();
					ata /= (scale_a * scale_a);
					ata = inv_w.asDiagonal() * ata * inv_w.asDiagonal();
					ata *= sigma;
#endif


				}


				time1.toc();
				//for (auto i : aj1) {
				//	int ci, cj, ck;
				//	ck = i % n;
				//	cj = i / n;
				//	ci = cj / n;
				//	cj = cj % n;

				//	int idx[3];
				//	idx[0] = ij(ci, cj, n);
				//	idx[1] = ij(ci, ck, n);
				//	idx[2] = ij(cj, ck, n);
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[1], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[2], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[2], sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[1], sigma));
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_v[i].row()<<" "<<test_v[i].col()<<" "<<test_v[i].value() << endl;
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_vb[i].row() << " " << test_vb[i].col() << " " << test_vb[i].value() << endl;
				//}
				//return;





				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;
				//cout << "rhsnorm1:" << rhs.norm() << endl;

				//cout << "rhsnorm2:" << rhs.norm() << endl;
				prt_res << iter_1 << "  " << rhs.norm() << "  ";
				double rhsnm = rhs.norm();

				if (iter_1 > 0 && rhs.norm() < rlim) {
					cout << "rhsout:" << rhs.norm() << endl;
					break;
				}
				//std::cout << "rhs:" << rhs.norm() << endl;
				Eigen::VectorXd direction;
				//if (flaga) {
				//	direction = get_direction(sigma, rhs, &aj1, &diag);
				//}
				//else {
				//	direction= get_direction(sigma, rhs, &aj1);
				//}
				int cgnum;
				time2.tic();
				int kcounter = 0;
				{
					double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);

					Eigen::VectorXd wx = sub_x + sub_v / sigma;

					Eigen::VectorXd r = wx;


					//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					//double ita = tol;
					//cg.setTolerance(ita / rhs.norm());
#ifdef DIAGMULT
					//ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
					diag += (a + cst) * Eigen::VectorXd::Ones(rhs.size());
#else
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
#endif
					if (func_p(wx) < 1 / sigma) {

						//ata += sigma * sub_w.cwiseProduct(sub_w).asDiagonal();

						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						//cg.compute(ata);

						//direction = -cg.solve(rhs);
						//kcounter = cg.iterations();
#ifdef DIAGMULT
						diag += sigma * Eigen::VectorXd::Ones(rhs.size());
						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = diag;



#else
						ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = ata.diagonal();
#endif




						newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);

						Eigen::VectorXd rnew = Eigen::VectorXd::Zero(rhs.size());

#ifdef DIAGMULT
						kcounter = psqmr_inf_diag(ata, inv_w, diag, newdiag, rnew, rhs, direction, sigma, tol_cg);
#else
						kcounter = psqmr_inf(ata, newdiag, rnew, rhs, direction, tol_cg);
#endif






					}
					else {

						if (true) {
#ifdef DIAGMULT
							diag +=1.0 / wx.norm() * Eigen::VectorXd::Ones(rhs.size());
							//ata += sigma * r.asDiagonal();
							//cout << "ata:" << ata.diagonal().norm() << "\tdiag:" << diag.norm() << endl;
#else
							ata += sigma/wx.norm() * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
#endif



							Eigen::VectorXd rnew = wx/pow(wx.norm(), 1.5);

							timer diana;
							diana.tic();
							timer ava;

#ifdef SMCG

							cg.compute(ata);
							ava.tic();
							direction = -cg.solve(rhs);
							ava.toc();
							kcounter = cg.iterations();
							ita = std::fmin(0.5, 0.000001 * rnew.norm());
							cg.setTolerance(ita / rnew.norm());
							ava.tic();
							Eigen::VectorXd ar = cg.solve(rnew);
							ava.toc();
							direction += (ar * rnew.dot(direction)) / (1 - rnew.dot(ar));
							kcounter += cg.iterations();
#else
#ifdef EIGENCG
							MatrixReplacement_inf AX(&ata, sub_var, &rnew, sigma);
							Eigen::ConjugateGradient<MatrixReplacement_inf, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg2;
							cg2.compute(AX);
							ita = std::fmin(0.5, 0.000001 * rhs.norm());
							cg2.setTolerance(ita / rhs.norm());
							direction = -cg2.solve(rhs);
							kcounter = cg2.iterations();
#else


							//diag = sigma * diag + (a + cst) * Eigen::VectorXd::Ones(diag.size());
							//diag += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r);
							//diag = Eigen::VectorXd::Ones(diag.size()).cwiseQuotient(diag);





#ifdef DIAGMULT
							Eigen::VectorXd newdiag = diag;
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);
							kcounter = psqmr_inf_diag(ata, inv_w, diag, newdiag, rnew, rhs, direction, sigma, tol_cg);
#else
							Eigen::VectorXd newdiag = ata.diagonal();
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);
							kcounter = psqmr_inf(ata, newdiag, rnew, rhs, direction, tol_cg);
#endif

							//cout << "cg:" << (ata * direction - rnew.dot(direction) * rnew + rhs).norm() << endl;
							//kcounter=psqmr(ata, rnew, rhs, direction);

#endif
#endif
							diana.toc();
							//cout << "cg:" << diana.blink * 100 << " " << kcounter << " " << diana.blink / kcounter * 10000 << endl;
							//cout << ava.time() * 100 << endl;


							//cout << "di:" << direction.norm() << endl;
							//return 0;
							//cout << "woodbury:" << direction.norm() << endl;
						}
						else {



						}
					}



					//std::cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

					//Eigen::PardisoLDLT<SparseMatrix<double> > cg;



					//cgnum = cg.iterations();
					//std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
					//prt_res << cg.iterations() << "\t";
					//std::cout << "pro_res:" << (A * direction+rhs).norm() << endl;
					//std::cout << "CG:" << cg.iterations() << std::endl;
					//std::cout << "di:" << direction.norm() << endl;

				}
				time2.toc();
				//cout << "di:" << setprecision(10) << direction.norm() << endl;
				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;
				Eigen::VectorXd wdirection = direction.cwiseProduct(inv_w);
				time3.tic();
				string tst(":");
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;

					double galp = 0;
					double g0 = direction.dot(rhs);
					if (g0 >= 0) {
						std::cout << "Need a descent direction\n";
						exit(-1);
					}
					//if (mthread == 0) {
					//	axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
					//	(*axbu)();
					//}
					//else {
					//	m_calc_phi(mthread, lx0, n, *axbu, 0, &direction);
					//}


					auto square = [](double x) {return x * x; };



					Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(sub_x.size());




					{


						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								i->y += (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mrhs(idx[2]) -= i->y / scale_a;
									mrhs(idx[0]) += i->y / scale_a;
								}
							}
						};

						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]))
							);
						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							nrhs += vrhs[i];
						}
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);



					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs.cwiseProduct(inv_w);
					nrhs += sigma * pi_c(sub_x + direction + sub_v / sigma, sigma) + cst * (sub_x + direction - x0);
					rhs = nrhs;
					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";

						rhsflag = true;
						lxflag = false;
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					Eigen::VectorXd lxnorm = Eigen::VectorXd::Zero(searchnum);
					//omp_set_num_threads(mthread);
					double lx0 = 0;
					if (lxflag == false) {
#pragma omp parallel for reduction(+:lx0)
						for (int i = 0; i < generator.data.size(); i++) {
							auto idx = generator.data[i].idx;
							float y = generator.data[i].y - (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;

							if (y > 0) {
								lx0 += square(y);
							}
						}
						lx0 = sigma / 2 * lx0 + sigma / 2 * pi_c(sub_x + sub_v / sigma, sigma).squaredNorm();

						lx0 += func_p(sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma)) + cst / 2 * (sub_x - x0).squaredNorm();

					}
					else {





						//	cout << "lxbias:" << lx0 - old_lx << endl;

						lx0 = old_lx;
					}










					//lx0 = lxnorm[0];




					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();
					alpha[0] = 1;

					for (int i = 1; i < searchnum; i++) {

						double lx = 0;

#pragma omp parallel for reduction(+:lx)
						for (int ii = 0; ii < generator.data.size(); ii++) {
							auto idx = generator.data[ii].idx;
							generator.data[ii].y -= (alpha[i - 1] - alpha[i]) * (wdirection[idx[0]] - wdirection[idx[1]] - wdirection[idx[2]]) / scale_a;

							if (generator.data[ii].y > 0) {
								lx += square(generator.data[ii].y);
							}

						}


						double alp = alpha[i];



						lx = sigma / 2 * lx;
						lx += sigma / 2 * pi_c(sub_x + alp * direction + sub_v / sigma, sigma).squaredNorm();
						lx += func_p(sub_x + alp * direction + sub_v / sigma - pi_c(sub_x + alp * direction + sub_v / sigma, sigma)) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();
						lxflag = true;
						old_lx = lx;
						if (lx <= down) {

							if (i == 1) {
								rhsflag = true;

							}
							else {
								rhsflag = false;
							}
							return alp;
						}

					}


					rhsflag = false;
					return alpha[searchnum - 1];





				}();
				time3.toc();
				//double alpha = 1;



				sub_x += dalpha * direction;
				//cout << "x:" << setprecision(10) << x.norm() << endl;
				//cout << "x:" << sub_x.norm() << endl;
				prt_res << dalpha << std::endl;
				std::cout << setprecision(10) << iter_1 << "\t" << rhsnm << "\t" << kcounter << "\t" << dalpha << endl;
			}


			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				consflag = true;
				double u_norm = 0;
#pragma omp parallel for reduction(+:kkt_y,kkt_u,y_norm,u_norm)
				for (int ii = 0; ii < generator.data.size(); ii++) {
					auto idx = generator.data[ii].idx;
					double y = generator.data[ii].y;
					float u = generator.data[ii].u;
					if (y < 0) {
						kkt_u += (u / sigma) * (u / sigma);
					}
					else {
						kkt_u += (y - u / sigma) * (y - u / sigma);
					}
					double ny = y < 0 ? y : 0;
					double nu = sigma * (y - ny);
					if (y - u / sigma > 0.01)consflag = false;

					//double y = axbu < 0 ? axbu : 0;
					y_norm += ny * ny;
					//	cout << "wtf\n";
						//*kkt_u += (axb - y) * (axb - y);


					double y_plus_u = (ny + nu > 0 ? ny : -nu);
					kkt_y += y_plus_u * y_plus_u;

					generator.data[ii].u = nu;
					u_norm += nu * nu;
				}


				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();

				total_kkt_u = kkt_u;
				total_y_norm = y_norm;


				//std::cout << "u:" << u.size() << endl;
				Eigen::VectorXd z = sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma);
				sub_v = sigma * pi_c(sub_x + sub_v / sigma, sigma);

				y_norm = sqrt(y_norm);
				kkt_u = sqrt(kkt_u) / (1 + y_norm);
				kkt_y = sqrt(kkt_y) / (1 + y_norm);

				std::cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;

				Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(sub_x.size());


				Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_x.size());



				//m_calc_atu(mthread, &atu, n, au);

				{
					int sz = generator.data.size();
					int pnum = mthread;


					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(sub_var);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i]), scale_a));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu.cwiseProduct(inv_w) + sub_v).norm() / (1 + scale_a * scale_a * sub_v.norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_x - z).norm() / scale_a / (1 + sub_x.norm() / scale_a);
				kkt_z = ((z - (sub_x + sub_v - pi_c(sub_x + sub_v, 1.0))).norm() / scale_a) / (1 + z.norm() / scale_a);
				//for (int i = 0; i < var_num; i++)
				//    t_res += abs(x[i] * w[i]);
				//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
				double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
				std::cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
				std::cout << "kkt:" << kkt_final << endl;
				runhist[string("kkt")] = kkt_final;
				rp = std::fmax(kkt_u, kkt_v);
				rd = kkt_x;
				kkt = kkt_final;

				min_kkt = min(kkt, min_kkt);
				tol_cg = min(min_kkt, 5.0e-3);
				rlim = max(1e-4, min(kkt, 1.0) * 0.9);
				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;
				dobj = 0;
#pragma omp parallel for reduction(+:dobj)
				for (int i = 0; i < generator.data.size(); i++) {
					int* idx = generator.data[i].idx;
					dobj += generator.data[i].u * (sub_b[idx[0]] - sub_b[idx[1]] - sub_b[idx[2]]);
				}
				dobj /= scale_a * scale_a;





			}
			time4.toc();



			double pobj = func_p(sub_x)/scale_a;
			cout << "pobj:" << pobj << ",dobj:" << dobj << endl;
			runhist[string("R_g")] = fabs((pobj - dobj) / (1 + dobj));
			if (kkt < tol && consflag&&fabs((pobj-dobj)/(1+dobj))<1e-4) {

				sigini = 0.01;

				cst_const = 0;

				//int pnum = mthread;


				//vector<thread>tv;

				//std::vector<almfunc_rhs>vf;
				//for (int i = 0; i < mthread; i++) {
				//	vf.push_back(almfunc_rhs(0, NULL, mrhs + i));
				//}
				//int sz = generator.data.size();
				//for (int i = 0; i < pnum; i++) {
				//	tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
				//		n, std::ref(sub_x), std::ref(sub_b), sigma));

				//}
				//for (int i = 0; i < pnum; i++)
				//	tv[i].join();
				//for (int i = 0; i < pnum; i++)
				//	nrhs += mrhs[i];
				//delete[]mrhs;
				break;
			}

			/*
			if (iter_1 == sub_iter) {
				sigma = 1;
				primal_win = 0;
				dual_win = 0;
			}
			else {
				if (5 * rp < rd)
					primal_win += 1;
				else if (rd < rp)
					dual_win += 1;


				if (iter_all % 3 == 2) {

					if (primal_win > std::fmax(1, 1.2 * dual_win)) {
						primal_win = 0;
						sigma = min(sigma / 1.2, 1.0e7);
					}
					else if ((dual_win > std::fmax(1, 1.2 * primal_win)) && (kkt > kkt0 * 0.7)) {

						dual_win = 0;
						sigma = max(sigma * 1.5, 1.0e0);
					}

				}
				else if (kkt > kkt0 * 0.7) {
					sigma = std::fmin(sigma * 1.5, 1.0e7);
				}
			}
			*/
			if (iter_1 == sub_iter && iter_all > 6) {


				sigma = max(sigma / 10, 1.0e-4);
			}

			else if (iter_1 != sub_iter) {
				if (kkt > kkt0*0.7)
					sigma = min(sigma * 1.5, 1.0e7);

			}
			kkt0 = kkt;





			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
			kkt0 = kkt;

		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;

		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]] / sub_w[colidx[i]];
				//b[i] = sub_b[colidx[i]];
				//w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}
	}


	//for (auto iter = u.begin(); iter != u.end(); iter++) {
	//	int64_t longindex = iter->first;
	//	int i, j, k;
	//	k = longindex % n;
	//	j = longindex / n;
	//	i = j / n;
	//	j = j % n;
	//	int idx[3];
	//	idx[0] = ij(i, j, n);
	//	idx[1] = ij(i, k, n);
	//	idx[2] = ij(j, k, n);
	//	dobj += iter->second * (b[idx[0]]-b[idx[1]]-b[idx[2]]);


	//}


	std::cout << "func_value:" << std::setprecision(8) << func_p(x) << "\t" << dobj << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	cout << "diana:" << time5.time() << endl;
	x /= scale_a;
	b /= scale_a;
	//x = x.cwiseQuotient(w);
	runhist[string("time_of_cg")] = generator.ava.time();

	return func_p(w.cwiseProduct(x));




}
double norm_ori::run_axbu(double cst_const, int con_const, int sub_iter) {
	scale_a = 3;
	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();
	b = b * scale_a;
	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 0.01;

	double c_t = 0;
	//std::unordered_map<int, int>colidx;

	timer time1, time2, time3, time4, time5;

	cout << "m:" << mthread << endl;
	int mxtnum = thread::hardware_concurrency();
	cons_generator generator(n, mxtnum, cgflg);
	generator.consvio_limit = scale_a * 0.01;
	double total_kkt_u = 0, total_y_norm = 0;
	double dobj = 0;
	double sigini = 1e-2;
	runhist[string("iter")] = -1;
	bool consflag = false;
	while (true) {
		runhist[string("iter")] += 1;
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(x, b, scale_a * tol);
		}
		sort(generator.data.begin(), generator.data.end(), [](cons_full& a, cons_full& b) {
			if (min(a.idx[0], min(a.idx[1], a.idx[2])) != min(b.idx[0], min(b.idx[1], b.idx[2]))) {
				return min(a.idx[0], min(a.idx[1], a.idx[2])) != min(b.idx[0], min(b.idx[1], b.idx[2]));
			}
			return max(a.idx[0], max(a.idx[1], a.idx[2])) != max(b.idx[0], max(b.idx[1], b.idx[2]));

			});
		double totalkkt = sqrt(generator.k1 + total_kkt_u) / (1 + sqrt(generator.k2 + total_y_norm));
		cout << "totalkkt:" << totalkkt << endl;
		if (generator.consvio)
			break;


		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = sigini;
		cst_const = sigini*1e-3;
		cout << "cst:" << cst_const << endl;
		int sub_var = colidx_num;
		std::cout << "total num of variables:" << var_num << " num of current variables:" << sub_var << endl;

		auto c_end = std::chrono::high_resolution_clock::now();
		auto c_eld = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start);
		c_t += double(c_eld.count()) * 1.0 * 1e-9;





		double	primal_win = 0;
		double dual_win = 0;
		double kkt0 = 1;
		kkt = 0.1;

		Eigen::VectorXd sub_x = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_b = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_w = Eigen::VectorXd::Zero(colidx_num);
		Eigen::VectorXd sub_v = Eigen::VectorXd::Zero(colidx_num);



		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				sub_x[colidx[i]] = x[i];
				sub_b[colidx[i]] = b[i];
				sub_w[colidx[i]] = w[i];
				sub_v[colidx[i]] = v[i];
			}
		}

		double kkt = 0.1;
		cout << sub_var << endl;
		double tol_cg = 5e-2;

		double min_kkt = 0.1;
		rlim = 5e-2;
		for (iter_all = 0; iter_all < 1000; iter_all++) {


			x0 = sub_x;
			cst = cst_const / sigma;



			//if (pwd != "skip" && iter_all > 0) {
			//	ofstream model(string("model_save/") + pwd + string("_") + graph_name + string("_") + to_string(iter_all) + string(".dat"), std::ios::out | std::ios::binary);
			//	auto i_end = std::chrono::high_resolution_clock::now();
			//	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
			//	double tm = double(eld.count()) * 1.0 * 1e-9;
			//	model << setprecision(searchnum);
			//	model << graph_name << endl;
			//	model << iter_all << endl;
			//	model << n << "\t" << var_num << "\t" << sigma << "\t" << tm + time_record << endl;



			//	for (int i = 0; i < var_num; i++)
			//		model << x[i] << endl;
			//	for (int i = 0; i < var_num; i++)
			//		model << v[i] << endl;
			//	model << u.size() << endl;
			//	for (auto i = u.begin(); i != u.end(); i++) {
			//		model << i->first << " " << i->second << endl;
			//	}
			//}

			std::cout << "start loop " << iter_all + 1 << ",sigma=" << sigma<<",rlim="<<rlim << endl;


			bool flaga = true;
			int iter_1;
			bool rhsflag = false;
			bool lxflag = false;
			double old_lx = 0;


#pragma omp parallel for
			for (int i = 0; i < generator.data.size(); i++) {

				generator.data[i].y = (sub_x[generator.data[i].idx[0]] - sub_x[generator.data[i].idx[1]] - sub_x[generator.data[i].idx[2]]
					+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a + generator.data[i].u / sigma;
			}


			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {


				//				double bias = 0;
				//#pragma omp parallel  for reduction(+:bias)
				//				for (int i = 0; i < generator.data.size(); i++) {
				//
				//					bias+=fabs(generator.data[i].y - (ori_x[generator.data[i].idx[0]] - ori_x[generator.data[i].idx[1]] - ori_x[generator.data[i].idx[2]]
				//						+ sub_b[generator.data[i].idx[0]] - sub_b[generator.data[i].idx[1]] - sub_b[generator.data[i].idx[2]]) / scale_a- generator.data[i].u / sigma);
				//				}
				//				cout << "bias:" << bias << endl;
								//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
								//std::cout << "x:" << x.norm() << endl;

				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double, Eigen::RowMajor>ata(sub_var, sub_var);

				{
#ifdef BDEG
					concurrency::concurrent_vector<int*>mtriv;
#else

					concurrency::concurrent_vector<Eigen::Triplet<double>>mtriv;
#endif
					diag = Eigen::VectorXd::Zero(sub_var);
					if (!rhsflag)
					{


						rhs = Eigen::VectorXd::Zero(sub_var);

						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mdiag[idx[1]] += 1;
									mrhs(idx[2]) -= i->y / scale_a;
									mdiag[idx[2]] += 1;
									mrhs(idx[0]) += i->y / scale_a;
									mdiag[idx[0]] += 1;
#ifdef BDEG
									mtriv.push_back(idx);
#else
									//fstream tst("t123.txt", std::ios::app);


									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[1], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[0], idx[2], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[0], -1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[1], idx[2], 1.0));
									mtriv.push_back(Eigen::Triplet<double>(idx[2], idx[1], 1.0));
									//tst << idx[0] << "\t" << idx[1]<<"\t"<<idx[2] << endl;
									//tst.close();
#endif
								}
							}
						};
						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var)), vdiag(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]), std::ref(vdiag[i])
							));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							rhs += vrhs[i];
							diag += vdiag[i];
						}
						rhs = sigma * rhs;
						rhs += sigma * prox_p_star(sub_x + sub_v / sigma, sigma, sub_w) + cst * (sub_x - x0);




					}
					else {

						vector<Eigen::VectorXd>vdiag(mthread, Eigen::VectorXd::Zero(sub_var));


						auto calc_diag = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mdiag) {
							for (auto i = begin; i != end; i++) {
								if (i->y > 0) {


									auto idx = i->idx;
									mdiag[idx[1]] += 1;
									mdiag[idx[2]] += 1;
									mdiag[idx[0]] += 1;


									mtriv.push_back(idx);
								}
							}
						};
						vector<thread>vt;
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calc_diag, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vdiag[i])));

						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							diag += vdiag[i];
						}




					}

					//
#ifdef BDEG
					buildeigen(mtriv, ata, diag);
#else
					ata.setFromTriplets(mtriv.begin(), mtriv.end());
#endif

					//cout << ata.norm() << " " << diag.norm() << endl;
#ifdef DIAGMULT

					diag /= (scale_a * scale_a);
					diag *= sigma;
					
#else

					ata += diag.asDiagonal();
					ata /= (scale_a * scale_a);
					ata *= sigma;
#endif


				}


				time1.toc();
				//for (auto i : aj1) {
				//	int ci, cj, ck;
				//	ck = i % n;
				//	cj = i / n;
				//	ci = cj / n;
				//	cj = cj % n;

				//	int idx[3];
				//	idx[0] = ij(ci, cj, n);
				//	idx[1] = ij(ci, ck, n);
				//	idx[2] = ij(cj, ck, n);
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[1], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[0], idx[2], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[0], -sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[1], idx[2], sigma));
				//	test_vb.push_back(Eigen::Triplet<double>(idx[2], idx[1], sigma));
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_v[i].row()<<" "<<test_v[i].col()<<" "<<test_v[i].value() << endl;
				//}
				//for (int i = 0; i < 10; i++) {
				//	std::cout << test_vb[i].row() << " " << test_vb[i].col() << " " << test_vb[i].value() << endl;
				//}
				//return;





				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;
				//cout << "rhsnorm1:" << rhs.norm() << endl;

				//cout << "rhsnorm2:" << rhs.norm() << endl;
				prt_res << iter_1 << "  " << rhs.norm() << "  ";
				double rhsnm = rhs.norm();

				if (iter_1 > 0 && rhs.norm() < rlim) {
					cout << "rhsout:" << rhs.norm() << endl;
					break;
				}
				//std::cout << "rhs:" << rhs.norm() << endl;
				Eigen::VectorXd direction;
				//if (flaga) {
				//	direction = get_direction(sigma, rhs, &aj1, &diag);
				//}
				//else {
				//	direction= get_direction(sigma, rhs, &aj1);
				//}
				int cgnum;
				time2.tic();
				int kcounter = 0;
				{
					double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);

					Eigen::VectorXd wx = sub_x + sub_v / sigma;

					Eigen::VectorXd r = wx;


					//Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					//double ita = tol;
					//cg.setTolerance(ita / rhs.norm());
#ifdef DIAGMULT
					//ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
					diag += (a + cst) * Eigen::VectorXd::Ones(rhs.size());
#else
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
#endif


						Eigen::VectorXd aj2 = (((sub_x + sub_v / sigma).cwiseAbs() - sub_w / sigma).array() <= 0.0).cast<double>().matrix();

						


#ifdef DIAGMULT
						diag += sigma * aj2;
						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						Eigen::VectorXd newdiag = diag;



#else
						ata += sigma * aj2.asDiagonal();
						Eigen::VectorXd newdiag = ata.diagonal();
#endif




						newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);


#ifdef DIAGMULT
						kcounter = psqmr_ori_diag(ata,  diag, newdiag,  rhs, direction, sigma, rlim);
#else
						kcounter = psqmr_ori(ata, newdiag, rhs, direction, rlim);
#endif






					
					



					//std::cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

					//Eigen::PardisoLDLT<SparseMatrix<double> > cg;



					//cgnum = cg.iterations();
					//std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
					//prt_res << cg.iterations() << "\t";
					//std::cout << "pro_res:" << (A * direction+rhs).norm() << endl;
					//std::cout << "CG:" << cg.iterations() << std::endl;
					//std::cout << "di:" << direction.norm() << endl;

				}
				time2.toc();
				//cout << "di:" << setprecision(10) << direction.norm() << endl;
				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;
				time3.tic();
				string tst(":");
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;

					double galp = 0;
					double g0 = direction.dot(rhs);
					if (g0 >= 0) {
						std::cout << "Need a descent direction\n";
						exit(-1);
					}
					//if (mthread == 0) {
					//	axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
					//	(*axbu)();
					//}
					//else {
					//	m_calc_phi(mthread, lx0, n, *axbu, 0, &direction);
					//}


					auto square = [](double x) {return x * x; };



					Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(sub_x.size());




					{


						auto calcrhs = [&](concurrency::concurrent_vector<cons_full>::iterator begin, concurrency::concurrent_vector<cons_full>::iterator end, Eigen::VectorXd& mrhs) {
							for (auto i = begin; i != end; i++) {
								auto idx = i->idx;
								i->y += (direction[idx[0]] - direction[idx[1]] - direction[idx[2]]) / scale_a;
								if (i->y > 0) {
									mrhs(idx[1]) -= i->y / scale_a;
									mrhs(idx[2]) -= i->y / scale_a;
									mrhs(idx[0]) += i->y / scale_a;
								}
							}
						};

						vector<thread>vt;
						vector<Eigen::VectorXd>vrhs(mthread, Eigen::VectorXd::Zero(sub_var));
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {
							vt.push_back(thread(calcrhs, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								std::ref(vrhs[i]))
							);
						}
						for (int i = 0; i < mthread; i++) {
							vt[i].join();
							nrhs += vrhs[i];
						}
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);



					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs;
					nrhs += sigma * prox_p_star(sub_x + direction + sub_v / sigma, sigma, sub_w) + cst * (sub_x + direction - x0);
					rhs = nrhs;
					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";

						rhsflag = true;
						lxflag = false;
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					Eigen::VectorXd lxnorm = Eigen::VectorXd::Zero(searchnum);
					//omp_set_num_threads(mthread);
					double lx0 = 0;
					if (lxflag == false) {
#pragma omp parallel for reduction(+:lx0)
						for (int i = 0; i < generator.data.size(); i++) {
							auto idx = generator.data[i].idx;
							float y = generator.data[i].y - (direction[idx[0]] - direction[idx[1]] - direction[idx[2]]) / scale_a;

							if (y > 0) {
								lx0 += square(y);
							}
						}
						lx0 = sigma / 2 * lx0 + sigma / 2 * prox_p_star(sub_x +  sub_v / sigma, sigma, sub_w).squaredNorm();

						lx0 += func_p(prox_p(sub_x + sub_v / sigma, sigma, sub_w), sub_w) + cst / 2 * (sub_x - x0).squaredNorm();

					}
					else {





						//	cout << "lxbias:" << lx0 - old_lx << endl;

						lx0 = old_lx;
					}










					//lx0 = lxnorm[0];




					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();
					alpha[0] = 1;

					for (int i = 1; i < searchnum; i++) {

						double lx = 0;

#pragma omp parallel for reduction(+:lx)
						for (int ii = 0; ii < generator.data.size(); ii++) {
							auto idx = generator.data[ii].idx;
							generator.data[ii].y -= (alpha[i - 1] - alpha[i]) * (direction[idx[0]] - direction[idx[1]] - direction[idx[2]]) / scale_a;

							if (generator.data[ii].y > 0) {
								lx += square(generator.data[ii].y);
							}

						}


						double alp = alpha[i];



						lx = sigma / 2 * lx;
						lx += sigma / 2 * prox_p_star(sub_x + alp*direction+sub_v / sigma, sigma, sub_w).squaredNorm();
						lx += func_p(prox_p(sub_x +alp*direction+ sub_v / sigma, sigma, sub_w),sub_w) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();
						lxflag = true;
						old_lx = lx;
						if (lx <= down) {

							if (i == 1) {
								rhsflag = true;

							}
							else {
								rhsflag = false;
							}
							return alp;
						}

					}


					rhsflag = false;
					return alpha[searchnum - 1];





				}();
				time3.toc();
				//double alpha = 1;



				sub_x += dalpha * direction;
				//cout << "x:" << setprecision(10) << x.norm() << endl;
				//cout << "x:" << sub_x.norm() << endl;
				prt_res << dalpha << std::endl;
				std::cout << setprecision(10) << iter_1 << "\t" << rhsnm << "\t" << kcounter << "\t" << dalpha << endl;
			}


			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				consflag = true;

#pragma omp parallel for reduction(+:kkt_y,kkt_u,y_norm)
				for (int ii = 0; ii < generator.data.size(); ii++) {
					auto idx = generator.data[ii].idx;
					double y = generator.data[ii].y;
					float u = generator.data[ii].u;
					if (y < 0) {
						kkt_u += (u / sigma) * (u / sigma);
					}
					else {
						kkt_u += (y - u / sigma) * (y - u / sigma);
					}
					double ny = y < 0 ? y : 0;
					double nu = sigma * (y - ny);


					if (y - u / sigma > 0.01)consflag = false;

					//double y = axbu < 0 ? axbu : 0;
					y_norm += ny * ny;
					//	cout << "wtf\n";
						//*kkt_u += (axb - y) * (axb - y);


					double y_plus_u = (ny + nu > 0 ? ny : -nu);
					kkt_y += y_plus_u * y_plus_u;

					generator.data[ii].u = nu;

				}


				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();

				total_kkt_u = kkt_u;
				total_y_norm = y_norm;
				//std::cout << "u:" << u.size() << endl;
				Eigen::VectorXd z = prox_p(sub_x + sub_v / sigma, sigma, sub_w);
				sub_v = sigma * prox_p_star(sub_x + sub_v / sigma, sigma, sub_w);

				y_norm = sqrt(y_norm);
				kkt_u = sqrt(kkt_u) / (1 + y_norm);
				kkt_y = sqrt(kkt_y) / (1 + y_norm);

				std::cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
				Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(sub_x.size());


				Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_x.size());



				//m_calc_atu(mthread, &atu, n, au);

				{
					int sz = generator.data.size();
					int pnum = mthread;


					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(sub_var);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i]), scale_a));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu + sub_v).norm() / (1 + scale_a * scale_a * sub_v.norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_x - z).norm() / scale_a / (1 + sub_x.norm() / scale_a);
				//kkt_z = ((z - prox_p(z + sub_v, 1.0, sub_w)).norm() / scale_a) / (1 + z.norm() / scale_a);
				Eigen::VectorXd tmp_xv = sub_x + sub_v;
				Eigen::VectorXd tmp_xv2= (tmp_xv.cwiseAbs()-sub_w).cwiseMax(0.0).cwiseProduct(tmp_xv.cwiseSign());

				kkt_z= ((z - prox_p(sub_x + sub_v, 1.0, sub_w)).norm() / scale_a) / (1 + z.norm() / scale_a);
				//for (int i = 0; i < var_num; i++)
				//    t_res += abs(x[i] * w[i]);
				//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
				double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
				std::cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
				std::cout << "kkt:" << kkt_final << endl;
				runhist[string("kkt")] = kkt_final;
				rp = std::fmax(kkt_u, kkt_v);
				rd = kkt_x;
				kkt = kkt_final;

				min_kkt = min(kkt, min_kkt);
				tol_cg = min(min_kkt, 5.0e-3);
				rlim = max(1e-4, min(kkt, 1.0) * 0.5);
				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;


				dobj = 0;
#pragma omp parallel for reduction(+:dobj)
				for (int i = 0; i < generator.data.size(); i++) {
					int* idx = generator.data[i].idx;
					dobj += generator.data[i].u * (sub_b[idx[0]] - sub_b[idx[1]] - sub_b[idx[2]]);
				}
				dobj /= scale_a * scale_a;




			}
			time4.toc();


			double pobj = func_p(sub_x, sub_w) / scale_a;
			cout << "pobj:" << pobj << ",dobj:" << dobj << endl;
			runhist[string("R_g")] = fabs((pobj - dobj) / (1 + dobj));
			if (kkt < tol && consflag&&fabs((pobj-dobj)/(1+dobj))<1e-4) {

				sigini = max(sigini/ 10, 1.0e-4);

				cst_const = sigini*1e-3;

				//int pnum = mthread;


				//vector<thread>tv;

				//std::vector<almfunc_rhs>vf;
				//for (int i = 0; i < mthread; i++) {
				//	vf.push_back(almfunc_rhs(0, NULL, mrhs + i));
				//}
				//int sz = generator.data.size();
				//for (int i = 0; i < pnum; i++) {
				//	tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
				//		n, std::ref(sub_x), std::ref(sub_b), sigma));

				//}
				//for (int i = 0; i < pnum; i++)
				//	tv[i].join();
				//for (int i = 0; i < pnum; i++)
				//	nrhs += mrhs[i];
				//delete[]mrhs;
				break;
			}

			/*
			if (iter_1 == sub_iter) {
				sigma = 1;
				primal_win = 0;
				dual_win = 0;
			}
			else {
				if (5 * rp < rd)
					primal_win += 1;
				else if (rd < rp)
					dual_win += 1;


				if (iter_all % 3 == 2) {

					if (primal_win > std::fmax(1, 1.2 * dual_win)) {
						primal_win = 0;
						sigma = min(sigma / 1.2, 1.0e7);
					}
					else if ((dual_win > std::fmax(1, 1.2 * primal_win)) && (kkt > kkt0 * 0.7)) {

						dual_win = 0;
						sigma = max(sigma * 1.5, 1.0e0);
					}

				}
				else if (kkt > kkt0 * 0.7) {
					sigma = std::fmin(sigma * 1.5, 1.0e7);
				}
			}
			*/
			if (iter_1 == sub_iter && iter_all > 5) {


				sigma = max(sigma / 10, 1.0e-4);
			}

			else if (iter_1 != sub_iter) {
				if (kkt > kkt0 * 0.7)
					sigma = min(sigma * 1.5, 1.0e7);

			}
			kkt0 = kkt;





			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
			kkt0 = kkt;

		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
		itervec.push_back(iter_all);
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]];
				//b[i] = sub_b[colidx[i]];
				//w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}
	}


	//for (auto iter = u.begin(); iter != u.end(); iter++) {
	//	int64_t longindex = iter->first;
	//	int i, j, k;
	//	k = longindex % n;
	//	j = longindex / n;
	//	i = j / n;
	//	j = j % n;
	//	int idx[3];
	//	idx[0] = ij(i, j, n);
	//	idx[1] = ij(i, k, n);
	//	idx[2] = ij(j, k, n);
	//	dobj += iter->second * (b[idx[0]]-b[idx[1]]-b[idx[2]]);


	//}


	std::cout << "func_value:" << std::setprecision(8) << func_p(x)/scale_a << "\t" << dobj << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	cout << "diana:" << time5.time() << endl;
	x /= scale_a;
	b /= scale_a;
	runhist[string("time_of_cg")] = generator.ava.time();
	//x = x.cwiseQuotient(w);
	return func_p(x);




}
