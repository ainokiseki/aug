


#include"metric.h"
#include"cons_generate.h"

#include<thread>
#include<mutex>
#include<math.h>
#include<numeric>
#include<string>
#include<sstream>
#include<concurrent_unordered_map.h>



#include<omp.h>




class norm_ori :public metric
{
public:

	norm_ori(string graph_name);
	//norm_ori(string graph_name, string data);
	//norm_ori(metric met) :metric(met) {

	//}
	~norm_ori();


	Eigen::VectorXd prox_p_star(Eigen::VectorXd x, double sigma, Eigen::VectorXd w);

	Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma, Eigen::VectorXd w);

	double func_p(Eigen::VectorXd x, Eigen::VectorXd w);

	inline Eigen::VectorXd prox_p(Eigen::VectorXd x, double sigma) ;
	inline Eigen::VectorXd prox_p_star (Eigen::VectorXd x, double sigma);
	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1);

	Eigen::VectorXd get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1, Eigen::VectorXd*diag);

	Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag);

	Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag, int ld);

	double get_step_length(Eigen::VectorXd& direction, Eigen::VectorXd& sub_x, Eigen::VectorXd& sub_v, Eigen::VectorXd& sub_w, double sigma, double cst);



	Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd& x, Eigen::VectorXd& w, Eigen::VectorXd& v, Eigen::VectorXd& rhs, int ld, double cst, int& cgnum);

	double run_axbu(double cst_const, int con_const, int sub_iter);

	double run_small(double cst_const, int con_const, int sub_iter);

	double run_ori(double cst_const, int sub_iter);

	//Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, int ld, double cst);

	//Eigen::VectorXd get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag, int ld, double cst);

	double run(double cst_const, int con_const, int sub_iter);

	double run_run(double cst_const, int con_const, int sub_iter);

	//void run(double cst_const=1e-3, int con_const=10);

	inline double func_p(Eigen::VectorXd x);

	double get_kkt(Eigen::VectorXd& sub_x, Eigen::VectorXd& sub_v, Eigen::VectorXd& sub_w, double sigma, double& rp, double& rd, std::vector<int64_t>* cons);



	

	double get_step_length(Eigen::VectorXd& direction, double sigma);
	



	Eigen::VectorXd get_direction_hand(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1, Eigen::VectorXd* diag);
	//void run();
	void run_dense();
	double time_record = 0;
	double cst = 0;
	int iter_all = 0;
	double sigma = 1;
	int mthread = 32;
	ofstream prt_res;
	Eigen::VectorXd rhs,x0;


	std::vector<std::vector< std::pair<int64_t, double> > > u;


};


class norm_2 :public metric
{
public:

	norm_2(string graph_name) :metric(graph_name) {
		mthread = std::thread::hardware_concurrency();
		iter_all = 0;
	}

	double time_record = 0;
	double cst = 0;
	int iter_all = 0;
	double sigma = 1;

	ofstream prt_res;
	Eigen::VectorXd rhs, x0;
	~norm_2() = default;

	double func_p(Eigen::VectorXd x) {
		return x.norm();
	}
	Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma) {
		Eigen::VectorXd tmp = x.cwiseAbs();
		double gap = tmp.norm() - 1 / sigma;
		if (gap <= 0)
			return x;
		return x / x.norm() * (1 / sigma);
	}
	double run_axbu(double cst_const, int con_const, int sub_iter);
};

class norm_inf :public metric
{
public:

	norm_inf(string graph_name) :metric(graph_name){
		mthread = std::thread::hardware_concurrency();
		iter_all = 0;
	}
	norm_inf(int x,string graph_name) :metric(x,graph_name) {
		mthread = std::thread::hardware_concurrency();
		iter_all = 0;
	}

	double time_record = 0;
	double cst = 0;
	int iter_all = 0;
	double sigma = 1;

	ofstream prt_res;
	Eigen::VectorXd rhs, x0;
	~norm_inf()=default;

	double func_p(Eigen::VectorXd x) {
		return x.cwiseAbs().maxCoeff();
	}
	Eigen::VectorXd pi_c(Eigen::VectorXd x, double sigma) {
		Eigen::VectorXd tmp = x.cwiseAbs();
		double gap = tmp.sum() - 1 / sigma;

		if (gap <= 0)
			return x;

		sort(tmp.data(), tmp.data() + tmp.size());
		double act;
		int len = tmp.size();
		for (int i = 0; i < len; i++) {

			if ((len - i) * tmp[i] >= gap) {
				act = gap / (len - i);
				break;
			}
			else {
				gap -= tmp[i];
			}
		}
		return x.cwiseSign().cwiseProduct(x - act * x.cwiseSign()).cwiseMax(0).cwiseProduct(x.cwiseSign());
	
	}

	double run_run(double cst_const, int con_const, int sub_iter);
	double run_run_run(double cst_const, int con_const, int sub_iter);
	double run_wx(double cst_const, int con_const, int sub_iter);
	double run_axbu_sparse(double cst_const, int con_const, int sub_iter);
	double run_axbu(double cst_const, int con_const, int sub_iter);
	double run_wx_test(double cst_const, int con_const, int sub_iter);
};






class cons_generator
{
public:
	cons_generator(int n, int mthread,bool cgflag=false) :n(n), mthread(mthread),cgflag(cgflag) {
		int var_num = n * (n - 1) / 2;
		colidx.resize(var_num, -1);
		colidx_num = 0;

	}
	~cons_generator() = default;
	timer ava;
	const bool cgflag;
	bool consvio=false;
	double consvio_limit = 1e-2;
	concurrency::concurrent_vector<cons_full>data;
	int64_t initconsnum = 0;
	double k1, k2;
	int n;
	int mthread;
	std::set<int64_t>viscons;
	std::vector<int> colidx;
	void alm_generate_para(Eigen::VectorXd& x, Eigen::VectorXd& b, double limit) {
		ReaderWriterLatch rw;
		mutex mu;
		auto gf = [&](int ia, int ib, double* kt1, double* kt2) {
			double kkt1 = 0, kkt2 = 0;
			for (int i = ia; i < ib; i++) {
				int64_t index;
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						for (int tst = 0; tst < 3; tst++) {

							swap(idx[tst], idx[0]);
							double test = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
							if (test > 0) {
								kkt1 += test * test;
							}
							else {
								kkt2 += test * test;
							}


							if (cgflag || x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > limit) {
								index = (static_cast<int64_t>(i) * n * n + j * n + k) * 4 + tst;
								rw.RLock();
								auto flg = viscons.find(index) == viscons.end();
								rw.RUnlock();
								if (flg) {

									for (int tp = 0; tp < 3; tp++) {
										if (colidx[idx[tp]] == -1) {
											mu.lock();
											if (colidx[idx[tp]] == -1) {


												colidx[idx[tp]] = colidx_num;
												colidx_num++;
											}
											mu.unlock();

										}
									}

									data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]],index));


									rw.WLock();
									viscons.insert(index);
									rw.WUnlock();
								}
							}
						}
					}
				}
			}
			*kt1 = kkt1;
			*kt2 = kkt2;

		};
		vector<double>kt1(mthread, 0), kt2(mthread, 0);
		vector<thread>vt;
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(gf, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), &kt1[i], &kt2[i]));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];
		}

	}
	void alm_generate_tri_para(Eigen::VectorXd& x, Eigen::VectorXd& b, double limit) {
		ReaderWriterLatch rw;
		mutex mu;
		auto gf = [&](int ia, int ib, double* kt1, double* kt2) {
			double kkt1 = 0, kkt2 = 0;
			for (int i = ia; i < ib; i++) {
				int64_t index;
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						index = static_cast<int64_t>(i) * n * n + j * n + k;
						double test = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
						if (test > 0) {
							kkt1 += test * test;
						}
						else {
							kkt2 += test * test;
						}
						test = x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]];
						if (test > 0) {
							kkt1 += test * test;
						}
						else {
							kkt2 += test * test;
						}
						test = x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]];
						if (test > 0) {
							kkt1 += test * test;
						}
						else {
							kkt2 += test * test;
						}


						if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > limit ||
							x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > limit ||
							x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > limit) {

							rw.RLock();
							auto flg = viscons.find(index) == viscons.end();
							rw.RUnlock();
							if (flg) {
								mu.lock();
								for (int tp = 0; tp < 3; tp++) {
									if (colidx[idx[tp]] == -1) {

										colidx[idx[tp]] = colidx_num;
										colidx_num++;

									}
								}
								mu.unlock();
								data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]]));
								data.push_back(cons_full(colidx[idx[1]], colidx[idx[0]], colidx[idx[2]]));
								data.push_back(cons_full(colidx[idx[2]], colidx[idx[1]], colidx[idx[0]]));


								rw.WLock();
								viscons.insert(index);
								rw.WUnlock();
							}
						}
					}
				}
			}
			*kt1 = kkt1;
			*kt2 = kkt2;

		};
		vector<double>kt1(mthread, 0), kt2(mthread, 0);
		vector<thread>vt;
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(gf, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), &kt1[i], &kt2[i]));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];
		}

	}
	void alm_generate_loop(Eigen::VectorXd& x, Eigen::VectorXd& b, double limit) {
		concurrency::concurrent_priority_queue<std::pair<double, int64_t> >consheap;
		mutex mu1,mu2;
		//std::unordered_set<int64_t>vviscons;
		//for (auto i : viscons) {
		//	vviscons.insert(i);
		//}
		timer jiaran;
		jiaran.tic();
		auto sqr = [](double x) {return x * x; };
		atomic<int64_t>ct(0);
		double maxv = 0;
		auto para_cons_generate = [&](int istart, int iend, double& kkt1, double& kkt2) {
			double tpkkt1 = 0, tpkkt2 = 0;

			for (int i = istart; i < iend; i++) {
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						for (int tst = 0; tst < 3; tst++) {
							swap(idx[0], idx[tst]);
							int64_t index = (static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k) * 4 + tst;

							double test = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];

							if ( cgflag||test > 0 ) {
#ifdef SHOWVIO
								mu1.lock();
								maxv = max(maxv, test);
								mu1.unlock();

#endif
								if (viscons.find(index) != viscons.end())continue;
								tpkkt1 += test * test;
								if (test < limit)continue;
								ct++;
								if (consheap.size() < initconsnum) {
									consheap.push(std::pair<double, int64_t>(-test, index));
								}
								else {
									std::pair<double, int64_t>tp;
									while (!consheap.try_pop(tp)) {};
									if (test > -tp.first) {
										consheap.push(std::pair<double, int64_t>(-test, index));
									}
									else {
										consheap.push(tp);
									}
								}


							}
							else if(test < 0) {
								tpkkt2 += test * test;
							}



						}
					}
				}
				kkt1 = tpkkt1; kkt2 = tpkkt2;
			}
		};
		std::vector<double>kt1(mthread, 0.0), kt2(mthread, 0.0), vax(mthread, 0.0);
		vector<thread>vt;
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(para_cons_generate, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), std::ref(kt1[i]), std::ref(kt2[i])));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];

		}
#ifdef SHOWVIO
		cout << "maxvio is:" << maxv << endl;
#endif
		//cout << "counteris:" << ct<<"ct2:"<<ct2 << endl;
		auto square = [](double x) {return x * x; };

		for (auto i_ : viscons) {
			int tst = i_ % 4;
			i_ /= 4;

			int i, j, k;
			k = i_ % n;
			j = i_ / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = ij(i, j, n);
			idx[1] = ij(i, k, n);
			idx[2] = ij(j, k, n);
			std::swap(idx[0], idx[tst]);
			double vio = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
			if (vio < 0) {
				k2 -= vio * vio;
				//total_ax -= square(x[idx[0]] - x[idx[1]] - x[idx[2]]);
			}
		}

		//cout << "kt2:" << ttk2 << endl;
		consvio = true;
		//para_cons_generate(0, n);
		std::pair<double, int64_t>tpp;
		while (consheap.size() > initconsnum) {
			consheap.try_pop(tpp);
		}
		bool exitflag = consheap.size()+data.size()>static_cast<int64_t>(n)*(n-1)/2;
		exitflag &= ct > data.size();
		//exitflag = false;
		cout << "cs:" << consheap.size() << "\t" << data.size() << "\t" << ct << endl;
		concurrency::concurrent_priority_queue<std::pair<double, int> >exitheap;
		if (exitflag) {
#pragma omp parallel for
			for (int i = 0; i < data.size(); i++) {
				int* idx = data[i].idx;
				double test = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
				if (test < 0) {
					exitheap.push(std::pair<double, int>(-test, i));
				}

			}
		}
		int exitnum = exitheap.size() / 2;
		cout << "exitnum:" << exitnum << endl;
		auto consgen = [&]() {
			std::pair<double, int64_t>tp{ 0.0,-1 };
			while (true) {
				while (!consheap.try_pop(tp)) {
					if (consheap.empty())return;
				}

				if (std::abs(tp.first) > consvio_limit)consvio = false;
				//if (tp.second == -1)cout << "wtf\n";
				int64_t index = tp.second;

				int tst = tp.second % 4;
				tp.second /= 4;

				int i, j, k;
				k = tp.second % n;
				j = tp.second / n;
				i = j / n;
				j = j % n;
				std::pair<double, int>tpnum{ 0.0,-1 };

				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				if (tst == 1)swap(idx[0], idx[1]);
				else if (tst == 2)swap(idx[0], idx[2]);

				if (exitflag && exitheap.size() > exitnum) {
					while (!exitheap.try_pop(tpnum)) {
						if (exitheap.empty())break;
					}
					if (tpnum.second == -1)goto L;
					mu1.lock();
					viscons.erase(data[tpnum.second].index);
					viscons.insert(index);
					mu1.unlock();
					for (int tmp = 0; tmp < 3; tmp++) {
						if (colidx[idx[tmp]] == -1) {
							mu2.lock();
							if (colidx[idx[tmp]] == -1) {
								colidx[idx[tmp]] = colidx_num;
								colidx_num++;
							}
							mu2.unlock();
						}
						int* idxori = data[tpnum.second].idx;
						idxori[0] = colidx[idx[0]];
						idxori[1] = colidx[idx[1]];
						idxori[2] = colidx[idx[2]];
						data[tpnum.second].y = 0;
						data[tpnum.second].u = 0;
						data[tpnum.second].index = index;
					}

					continue;

				}
				L:



				for (int tmp = 0; tmp < 3; tmp++) {
					if (colidx[idx[tmp]] == -1) {
						mu2.lock();
						if (colidx[idx[tmp]] == -1) {
							colidx[idx[tmp]] = colidx_num;
							colidx_num++;
						}
						mu2.unlock();
					}
				}

				//cout << idx[1] << " " << idx[2] << " " << idx[0] << endl;
				mu1.lock();
				viscons.insert(index);
				mu1.unlock();
				data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]], index));


			}

		};
		vt.clear();
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(consgen));
		}
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
		}


	}



	std::atomic<int>colidx_num;
	bool initflag = true;
	void generate(Eigen::VectorXd& x, Eigen::VectorXd& b, double limit) {
		ava.tic();

		if (initflag) {

			alm_generate_para(x, b, limit);
			
			//initconsnum = 70000000;
			initflag = false;
		}
		else {

			if (cgflag) {
				consvio = true;
			}
			else {



				alm_generate_loop(x, b, limit);
			}
		}

		initconsnum = data.size();



		ava.toc();
		std::cout << "cg time:" << ava.time() << std::endl;
	}
};

class cons_generator_sparse
{
public:
	cons_generator_sparse(int n, int mthread) :n(n), mthread(mthread) {
		int var_num = n * (n - 1) / 2;
		colidx_num = 0;

	}
	~cons_generator_sparse() = default;
	timer ava;
	bool consvio = false;
	double consvio_limit = 1e-2;
	concurrency::concurrent_vector<cons_full>data;
	int64_t initconsnum = 0;
	double k1, k2;
	int n;
	int mthread;
	std::set<int64_t>viscons;
	//std::unordered_map<int,int> colidx;
	concurrency::concurrent_unordered_map<int, int>colidx;
	template<typename T1, typename T2>
	void alm_generate_para(T1& x, T2& b, double limit) {
		ReaderWriterLatch rw,rw2;
		mutex mu;
		auto gf = [&](int ia, int ib, double* kt1, double* kt2) {
			double kkt1 = 0, kkt2 = 0;
			for (int i = ia; i < ib; i++) {
				int64_t index;
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						for (int tst = 0; tst < 3; tst++) {

							swap(idx[tst], idx[0]);



							if (b[idx[0]] - b[idx[1]] - b[idx[2]] > 0) {
								index = (static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k) * 4 + tst;
								rw.RLock();
								auto flg = viscons.find(index) == viscons.end();
								rw.RUnlock();
								if (flg) {

									for (int tp = 0; tp < 3; tp++) {
										
										auto cfind = colidx.find(idx[tp]);
										
										if (cfind == colidx.end()) {
											rw2.WLock();

											if (colidx.find(idx[tp]) == colidx.end()) {

												colidx[idx[tp]] = colidx_num;
												colidx_num++;
											}
											rw2.WUnlock();

										}
									}

									data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]],index));


									rw.WLock();
									viscons.insert(index);
									rw.WUnlock();
								}
							}
						}
					}
				}
			}
			*kt1 = kkt1;
			*kt2 = kkt2;

		};

		vector<double>kt1(mthread, 0), kt2(mthread, 0);
		vector<thread>vt;
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(gf, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), &kt1[i], &kt2[i]));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];
		}

	}
	template<typename T1, typename T2>
	void alm_generate_para_new(T1& x, T2& b, double limit,vector<pair<int,int>>&wmap) {
		ReaderWriterLatch rw, rw2;
		mutex mu;
		auto gf = [&](int ia, int ib, double* kt1, double* kt2) {
			double kkt1 = 0, kkt2 = 0;
			for (int it = ia; it < ib; it++) {
				int ii, jj;
				ii = wmap[it].first;
				jj = wmap[it].second;
				int64_t index;

					for (int kk = 0; kk < n; kk++) {
						if (kk == ii || kk == jj)continue;
						int i, j, k;
						i = min(min(ii, jj), kk);
						k = max(max(ii, jj), kk);
						j = ii + kk + jj - i - k;
						int idx[3];
						idx[0] = ij(i, j, static_cast<int64_t>(n));
						idx[1] = ij(i, k, static_cast<int64_t>(n));
						idx[2] = ij(j, k, static_cast<int64_t>(n));
						if (b[idx[0]] == 1 && b[idx[1]] == 1 && b[idx[2]] == 1)continue;
							for (int tst = 0; tst < 3; tst++) {

								swap(idx[tst], idx[0]);



								if (b[idx[0]] - b[idx[1]] - b[idx[2]] > 0) {
									index = (static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k) * 4 + tst;
									rw.RLock();
									auto flg = viscons.find(index) == viscons.end();
									rw.RUnlock();
									if (flg) {

										for (int tp = 0; tp < 3; tp++) {

											auto cfind = colidx.find(idx[tp]);

											if (cfind == colidx.end()) {
												rw2.WLock();

												if (colidx.find(idx[tp]) == colidx.end()) {
													colidx[idx[tp]] = colidx_num;
													colidx_num++;
												}
												rw2.WUnlock();

											}
										}

										

										bool flg = false;
										rw.WLock();
										if (viscons.find(index) == viscons.end()) {
											viscons.insert(index);
											flg = true;
										}
										rw.WUnlock();
										if(flg)data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]], index));

									}
								}
							}
					}
				}
			*kt1 = kkt1;
			*kt2 = kkt2;

		};

		vector<double>kt1(mthread, 0), kt2(mthread, 0);
		vector<thread>vt;
		int sz = wmap.size();
		cout << "wsize:" << sz << endl;
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(gf, i* sz/ mthread,(i+1)*sz/mthread  , &kt1[i], &kt2[i]));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];
		}

	}
	template<typename T1, typename T2>
	void alm_generate_loop(T1& x, T2& b, double limit) {
		concurrency::concurrent_priority_queue<std::pair<double, int64_t> >consheap;
		mutex mu1, mu2;
		ReaderWriterLatch rw, rw2;
		//std::unordered_set<int64_t>vviscons;
		//for (auto i : viscons) {
		//	vviscons.insert(i);
		//}
		timer jiaran;
		jiaran.tic();
		auto sqr = [](double x) {return x * x; };
		atomic<int64_t>ct(0);
		auto para_cons_generate = [&](int istart, int iend, double& kkt1, double& kkt2) {
			double tpkkt1 = 0, tpkkt2 = 0;

			for (int i = istart; i < iend; i++) {
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						for (int tst = 0; tst < 3; tst++) {
							swap(idx[0], idx[tst]);
							int64_t index = (static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k) * 4 + tst;

							double test = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];

							if (test > 0) {

								
								//ct++;
								if (viscons.find(index) != viscons.end())continue;
								tpkkt1 += test * test;
								if (consheap.size() < initconsnum) {
									consheap.push(std::pair<double, int64_t>(-test, index));
								}
								else {
									std::pair<double, int64_t>tp;
									while (!consheap.try_pop(tp)) {};
									if (test > -tp.first) {
										consheap.push(std::pair<double, int64_t>(-test, index));
									}
									else {
										consheap.push(tp);
									}
								}


							}
							else {
								tpkkt2 += test * test;
							}



						}
					}
				}
				kkt1 = tpkkt1; kkt2 = tpkkt2;
			}
		};
		std::vector<double>kt1(mthread, 0.0), kt2(mthread, 0.0), vax(mthread, 0.0);
		vector<thread>vt;
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(para_cons_generate, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), std::ref(kt1[i]), std::ref(kt2[i])));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];

		}
		//cout << "counteris:" << ct << endl;
		auto square = [](double x) {return x * x; };

		for (auto i_ : viscons) {
			int tst = i_ % 4;
			i_ /= 4;

			int i, j, k;
			k = i_ % n;
			i_ = i_ / n;
			i = i_ / n;
			j = i_ % n;
			int idx[3];
			idx[0] = ij(i, j, n);
			idx[1] = ij(i, k, n);
			idx[2] = ij(j, k, n);
			std::swap(idx[0], idx[tst]);
			double vio = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
			if (vio < 0) {
				k2 -= vio * vio;
				//total_ax -= square(x[idx[0]] - x[idx[1]] - x[idx[2]]);
			}
		}

		//cout << "kt2:" << ttk2 << endl;
		consvio = true;
		//para_cons_generate(0, n);
		std::pair<double, int64_t>tpp;
		while (consheap.size() > initconsnum) {
			consheap.try_pop(tpp);
		}
		//exitflag = true;


		auto consgen = [&]() {
			std::pair<double, int64_t>tp{ 0.0,-1 };
			while (true) {
				while (!consheap.try_pop(tp)) {
					if (consheap.empty())return;
				}
#ifndef MDEBUG
				if (std::abs(tp.first) < limit)continue;
#endif
				if (std::abs(tp.first) > consvio_limit)consvio = false;
				//if (tp.second == -1)cout << "wtf\n";
				int64_t index = tp.second;

				int tst = tp.second % 4;
				tp.second /= 4;

				int i, j, k;
				k = tp.second % n;
				tp.second /=  n;
				i = tp.second / n;
				j = tp.second % n;
				std::pair<double, int>tpnum{ 0.0,-1 };

				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				if (tst == 1)swap(idx[0], idx[1]);
				else if (tst == 2)swap(idx[0], idx[2]);



				for (int tmp = 0; tmp < 3; tmp++) {
					
					auto cfind = colidx.find(idx[tmp]);
					

					if (cfind == colidx.end()) {
						mu2.lock();
						if (colidx.find(idx[tmp]) == colidx.end()) {
							colidx[idx[tmp]] = colidx_num;
							colidx_num++;
						}
						mu2.unlock();
					}
				}

				//cout << idx[1] << " " << idx[2] << " " << idx[0] << endl;
				mu1.lock();
				viscons.insert(index);
				mu1.unlock();
				data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]], index));


			}

		};
		vt.clear();
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(consgen));
		}
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
		}


	}
	
	template<typename T1, typename T2>
	void alm_generate_loop(T1& x, T2& b, double limit, vector<pair<int, int>>& wmap) {
		concurrency::concurrent_priority_queue<std::pair<double, int64_t> >consheap;
		mutex mu1, mu2;
		concurrency::concurrent_unordered_set<int64_t>subindex;
		ReaderWriterLatch rw, rw2;
		//std::unordered_set<int64_t>vviscons;
		//for (auto i : viscons) {
		//	vviscons.insert(i);
		//}
		timer jiaran;
		jiaran.tic();
		auto sqr = [](double x) {return x * x; };
		atomic<int64_t>ct(0);
		auto para_cons_generate = [&](int ia, int ib, double& kkt1, double& kkt2) {
			double tpkkt1 = 0, tpkkt2 = 0;

			for (int it = ia; it < ib; it++) {
				int ii, jj;
				ii = wmap[it].first;
				jj = wmap[it].second;
				int64_t index;
				
					for (int kk = 0; kk < n; kk++) {
						if (kk == ii || kk == jj)continue;
						int i, j, k;
						i = min(min(ii, jj), kk);
						k = max(max(ii, jj), kk);
						j = ii + kk + jj - i - k;

						int idx[3];
						idx[0] = ij(i, j, static_cast<int64_t>(n));
						idx[1] = ij(i, k, static_cast<int64_t>(n));
						idx[2] = ij(j, k, static_cast<int64_t>(n));
						for (int tst = 0; tst < 3; tst++) {
							swap(idx[0], idx[tst]);
							int64_t index = (static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k) * 4 + tst;

							double test = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];

							if (test > 0) {

								//ct++;
								if (viscons.find(index) != viscons.end())continue;
								bool flg;
								rw.RLock();
								flg = (subindex.find(index) == subindex.end());
								rw.RUnlock();
								if (flg) {
									flg = false;
									rw.WLock();
									if (subindex.find(index) == subindex.end()) {
										subindex.insert(index);
										flg = true;
									}
									rw.WUnlock();
								}
								if (!flg)continue;

								tpkkt1 += test * test;
								if (consheap.size() < initconsnum) {
									consheap.push(std::pair<double, int64_t>(-test, index));
								}
								else {
									std::pair<double, int64_t>tp;
									while (!consheap.try_pop(tp)) {};
									if (test > -tp.first) {
										consheap.push(std::pair<double, int64_t>(-test, index));
									}
									else {
										consheap.push(tp);
									}
								}
							}
							else {
								tpkkt2 += test * test;
							}
						}
					}
				
				kkt1 = tpkkt1; kkt2 = tpkkt2;
			}
		};
		std::vector<double>kt1(mthread, 0.0), kt2(mthread, 0.0), vax(mthread, 0.0);
		vector<thread>vt;
		int64_t sz = wmap.size();
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(para_cons_generate, sz*i / mthread,  sz*(i+1) / mthread, std::ref(kt1[i]), std::ref(kt2[i])));
		}
		k1 = k2 = 0;
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
			k1 += kt1[i];
			k2 += kt2[i];

		}
		//cout << "counteris:" << ct << endl;
		auto square = [](double x) {return x * x; };

		for (auto i_ : viscons) {
			int tst = i_ % 4;
			i_ /= 4;

			int i, j, k;
			k = i_ % n;
			i_ = i_ / n;
			i = i_ / n;
			j = i_ % n;
			int idx[3];
			idx[0] = ij(i, j, n);
			idx[1] = ij(i, k, n);
			idx[2] = ij(j, k, n);
			std::swap(idx[0], idx[tst]);
			double vio = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
			if (vio < 0) {
				k2 -= vio * vio;
				//total_ax -= square(x[idx[0]] - x[idx[1]] - x[idx[2]]);
			}
		}

		//cout << "kt2:" << ttk2 << endl;
		consvio = true;
		//para_cons_generate(0, n);
		std::pair<double, int64_t>tpp;
		while (consheap.size() > initconsnum) {
			consheap.try_pop(tpp);
		}
		//exitflag = true;


		auto consgen = [&]() {
			std::pair<double, int64_t>tp{ 0.0,-1 };
			while (true) {
				while (!consheap.try_pop(tp)) {
					if (consheap.empty())return;
				}
#ifndef MDEBUG
				if (std::abs(tp.first) < limit)continue;
#endif
				if (std::abs(tp.first) > consvio_limit)consvio = false;
				//if (tp.second == -1)cout << "wtf\n";
				int64_t index = tp.second;

				int tst = tp.second % 4;
				tp.second /= 4;

				int i, j, k;
				k = tp.second % n;
				tp.second /= n;
				i = tp.second / n;
				j = tp.second % n;
				std::pair<double, int>tpnum{ 0.0,-1 };

				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				if (tst == 1)swap(idx[0], idx[1]);
				else if (tst == 2)swap(idx[0], idx[2]);



				for (int tmp = 0; tmp < 3; tmp++) {

					auto cfind = colidx.find(idx[tmp]);


					if (cfind == colidx.end()) {
						mu2.lock();
						if (colidx.find(idx[tmp]) == colidx.end()) {
							colidx[idx[tmp]] = colidx_num;
							colidx_num++;
						}
						mu2.unlock();
					}
				}

				//cout << idx[1] << " " << idx[2] << " " << idx[0] << endl;
				mu1.lock();
				viscons.insert(index);
				mu1.unlock();
				data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]], index));


			}

		};
		vt.clear();
		for (int i = 0; i < mthread; i++) {
			vt.push_back(thread(consgen));
		}
		for (int i = 0; i < mthread; i++) {
			vt[i].join();
		}


	}

	template<typename T1>
	void update_wmap(T1& x, vector<pair<int, int>>& wmap) {
		for (auto& i : x.data) {
			auto rs = rij(i.first, n);
			wmap.emplace_back(rs);
		}
		sort(wmap.begin(), wmap.end());
		wmap.resize(unique(wmap.begin(), wmap.end())-wmap.begin());
	}

	std::atomic<int>colidx_num;
	bool initflag = true;

	template<typename T1,typename T2>
	void generate(T1& x,T2& b, double limit, vector<pair<int, int>>& wmap) {
		ava.tic();

		if (initflag) {

			alm_generate_para_new(x, b, limit,wmap);

			//initconsnum = 70000000;
			initflag = false;
		}
		else {


			update_wmap(x, wmap);
			cout << "wmap updated, size:" << wmap.size() << endl;
			alm_generate_loop(x, b, limit,wmap);

		}

		initconsnum = data.size();



		ava.toc();
		std::cout << "cg time:" << ava.time() << std::endl;
	}
};
