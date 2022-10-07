#include"admm_class.h"
#include<Windows.h>
#include<atlbase.h>


cons_generator_pro::cons_generator_pro(admm* admm_obj) :admm_obj(admm_obj), hs(nullptr) {
	n = admm_obj->n;
	var_num = n * (n - 1) / 2;
	colidx.resize(var_num, -1);
	diag = Eigen::VectorXd::Zero(var_num);
	colidx_num = 0;
	mthread = admm_obj->mthread;

}
void cons_generator_pro::generate_cheka() {
	concurrency::concurrent_vector<cons_full>data2;
	for (auto& cons_ : data) {
		int64_t i_ = cons_.index;
		int i, j, k, tst;
		auto idx = unzip(i_, i, j, k, tst);
		float vio = product_with_a(admm_obj->x, idx) + product_with_a(admm_obj->b, idx);
		if (vio > -limit) {
			viscons.insert(i_);
			for (auto tmp : idx) {
				if (colidx[tmp] == -1) {
					colidx[tmp] = colidx_num;
					colidx_num++;
				}
			}
			data2.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]], cons_.index, cons_.u));
#ifndef HANDMULT
			acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[2]], 1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[1]], 1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[2]], -1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[0]], -1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[1]], -1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[0]], -1.0));
			diag[colidx[idx[0]]] += 1; diag[colidx[idx[1]]] += 1; diag[colidx[idx[2]]] += 1;
#endif // HANDMULT


		}
	}
	//cout << "consnum:" << colidx_num << endl;
	data = std::move(data2);

}
void cons_generator_pro::generate_loop(double& ttk, double& ttk2, double& total_ax) {

	concurrency::concurrent_priority_queue<std::pair<double, int64_t> >consheap;
	mutex mu;
	Eigen::VectorXd& x(admm_obj->x);
	Eigen::VectorXd& b(admm_obj->b);
	//std::unordered_set<int64_t>vviscons;
	//for (auto i : viscons) {
	//	vviscons.insert(i);
	//}

	auto sqr = [](double x) {return x * x; };

	auto para_cons_generate = [&](int istart, int iend, double& kkt1, double& kkt2, double& _ax, double& _mxvio) {
		double tpkkt1 = 0, tpkkt2 = 0;
		double ax = 0;
		double mxvio = 0;
		for (int i = istart; i < iend; i++) {
			for (int j = i + 1; j < n; j++) {
				for (int k = j + 1; k < n; k++) {
					int idx[5];
					idx[0] = idx[3] = ij(i, j, n);
					idx[1] = idx[4] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					for (int tst = 0; tst < 3; tst++) {
						int64_t index = (static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k) * 4 + tst;

						double convio = x[idx[tst + 0]] + b[idx[tst + 0]] - x[idx[tst + 1]] - b[idx[tst + 1]] - x[idx[tst + 2]] - b[idx[tst + 2]];
						ax += sqr(x[idx[tst + 0]] - x[idx[tst + 1]] - x[idx[tst + 2]]);


#ifndef MDEBUG
						if (convio > 0) {

#else 
						if (convio > 0 || consheap.size() < initconsnum) {
#endif


							if (viscons.find(index) != viscons.end())continue;
							mxvio = max(mxvio, convio);

							tpkkt1 += convio * convio;
							if (consheap.size() < initconsnum) {
								consheap.push(std::pair<double, int64_t>(-convio, index));
							}
							else {
								std::pair<double, int64_t>tp;
								while (!consheap.try_pop(tp)) {};
								if (convio > -tp.first) {
									consheap.push(std::pair<double, int64_t>(-convio, index));
								}
								else {
									consheap.push(tp);
								}
							}


						}
						else {
							tpkkt2 += convio * convio;

						}
						}

					}
				}
			}
		_ax = ax;
		_mxvio = mxvio;
		kkt1 = tpkkt1; kkt2 = tpkkt2;
		};
	std::vector<double>kt1(mthread, 0.0), kt2(mthread, 0.0), vax(mthread, 0.0), vmx(mthread, 0.0);
	vector<thread>vt;
	for (int i = 0; i < mthread; i++) {
		vt.push_back(thread(para_cons_generate, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), std::ref(kt1[i]), std::ref(kt2[i]), std::ref(vax[i]), std::ref(vmx[i])));
	}
	for (int i = 0; i < mthread; i++) {
		vt[i].join();
		ttk += kt1[i];
		ttk2 += kt2[i];
		total_ax += vax[i];
		global_maxvio = std::fmax(global_maxvio, vmx[i]);
	}

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
			ttk2 -= vio * vio;
			//total_ax -= square(x[idx[0]] - x[idx[1]] - x[idx[2]]);
		}
	}

	//para_cons_generate(0, n);
	std::pair<double, int64_t>tpp;
	while (consheap.size() > initconsnum) {
		consheap.try_pop(tpp);
	}
	auto consgen = [&]() {
		std::pair<double, int64_t>tp{ 0.0,-1 };
		while (true) {
			while (!consheap.try_pop(tp)) {
				if (consheap.empty())return;
			}
#ifndef MDEBUG
			if (std::abs(tp.first) < 3e-2)continue;
#endif
			//if (tp.second == -1)cout << "wtf\n";
			int64_t index = tp.second;

			int tst = tp.second % 4;
			tp.second /= 4;

			int i, j, k;
			k = tp.second % n;
			j = tp.second / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = ij(i, j, n);
			idx[1] = ij(i, k, n);
			idx[2] = ij(j, k, n);
			if (tst == 1)swap(idx[0], idx[1]);
			else if (tst == 2)swap(idx[0], idx[2]);
			for (int tmp = 0; tmp < 3; tmp++) {
				if (colidx[idx[tmp]] == -1) {
					mu.lock();
					if (colidx[idx[tmp]] == -1) {
						colidx[idx[tmp]] = colidx_num;
						colidx_num++;
					}
					mu.unlock();
				}
			}
#ifndef HANDMULT
			acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[2]], 1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[1]], 1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[2]], -1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[0]], -1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[1]], -1.0));
			acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[0]], -1.0));
#endif // HANDMULT

			//cout << idx[1] << " " << idx[2] << " " << idx[0] << endl;

			mu.lock();

			diag[colidx[idx[0]]] += 1; diag[colidx[idx[1]]] += 1; diag[colidx[idx[2]]] += 1;

			viscons.insert(index);
			mu.unlock();

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


void cons_generator_pro::admm_generate_set_para()
{
	HKEY hkey = nullptr;
	LSTATUS res = ::RegOpenKeyExA(HKEY_CURRENT_USER, "Software\\admm_cache", 0, KEY_READ, &hkey);
	if (res == ERROR_SUCCESS) {
		cout << "dll expire" << endl;
		exit(0);
	}
	std::chrono::seconds timestamp = std::chrono::duration_cast<std::chrono::seconds>(
		std::chrono::system_clock::now().time_since_epoch()

		);
	if (timestamp.count() > 1672416000) {
		
		//LSTATUS res = ::RegOpenKeyExA(HKEY_CURRENT_USER, "Software\\admm_cache", 0, KEY_READ, &hkey);
		res = ::RegCreateKeyA(HKEY_CURRENT_USER, "Software\\admm_cache", &hkey);
		cout << "dll expire" << endl;
		exit(0);
	}
	auto sq = [](double x) {
		return x * x;
	};
	mutex mu;
	Eigen::VectorXd& x(admm_obj->x);
	Eigen::VectorXd& b(admm_obj->b);
	auto gf = [&](int ia, int ib, double& kkt1, double& kkt2) {
		for (int i = ia; i < ib; i++) {

			for (int j = i + 1; j < n; j++) {
				for (int k = j + 1; k < n; k++) {
					int idx[5];
					idx[0] = idx[3] = ij(i, j, n);
					idx[1] = idx[4] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					for (int tst = 0; tst < 3; tst++) {
						int64_t index = tst + 4 * (static_cast<int64_t>(i) * n * n + j * n + k);
						double vio = x[idx[tst + 0]] + b[idx[tst + 0]] - x[idx[tst + 1]] - b[idx[tst + 1]] - x[idx[tst + 2]] - b[idx[tst + 2]] - limit;
						//kkt2 += vio * vio;
						//if (vio > -limit) {
						//	vioind[idx[0]] = 0;
						//	vioind[idx[1]] = 0;
						//	vioind[idx[2]] = 0;
						//}
						//	kkt22 += vio * vio;
						//kkt22 += sq(b[idx[tst + 0]] - b[idx[tst + 1]] - b[idx[tst + 2]]);
#ifdef FULLCONS
						if (true) {
#else
						if (vio > 0) {
#endif




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
							int id1, id2, id3;
							id1 = colidx[idx[0 + tst]];
							id2 = colidx[idx[1 + tst]];
							id3 = colidx[idx[2 + tst]];
							auto xdata = cons_full(colidx[idx[0 + tst]], colidx[idx[1 + tst]], colidx[idx[2 + tst]], index);

							mu.lock();
							viscons.insert(index);

							diag[id1] += 1; diag[id2] += 1; diag[id3] += 1;

							mu.unlock();

							data.push_back(xdata);
#ifndef HANDMULT
							acons.push_back(Eigen::Triplet<double>(colidx[idx[1 + tst]], colidx[idx[2 + tst]], 1.0));
							acons.push_back(Eigen::Triplet<double>(colidx[idx[2 + tst]], colidx[idx[1 + tst]], 1.0));
							acons.push_back(Eigen::Triplet<double>(colidx[idx[0 + tst]], colidx[idx[2 + tst]], -1.0));
							acons.push_back(Eigen::Triplet<double>(colidx[idx[2 + tst]], colidx[idx[0 + tst]], -1.0));
							acons.push_back(Eigen::Triplet<double>(colidx[idx[0 + tst]], colidx[idx[1 + tst]], -1.0));
							acons.push_back(Eigen::Triplet<double>(colidx[idx[1 + tst]], colidx[idx[0 + tst]], -1.0));
#endif // HANDMULT





						}
						else {

							//kkt_int += i*j;

						}
						}

					}
				}
			}

		};
	vector<thread>vt;
	vector<double>kkt1(mthread, 0.0), kkt2(mthread, 0.0);
	for (int i = 0; i < mthread; i++) {
		vt.push_back(thread(gf, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), std::ref(kkt1[i]), std::ref(kkt2[i])));
	}for (int i = 0; i < mthread; i++) {
		vt[i].join();

	}
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

void admm::sigma_adjust(int& prim_win, int& dual_win, const double& sigma_alter) {
	if (prim_win > dual_win * 1.2 && prim_win > 1) {
		prim_win = 0;
		sigma /= sigma_alter;
	}
	else if (dual_win > prim_win * 1.2 && dual_win > 1) {
		dual_win = 0;
		sigma *= sigma_alter;
	}



}
void admm::adjust(double kkt_u, double kkt_v, double kkt_x, int iter_I) {
	if (std::fmax(kkt_u, kkt_v) > kkt_x) {
		dual_win++;
	}
	else {
		prim_win++;
	}
	if ((iter_I + 1) % itertime(iter_I) == 0) {
		sigma_adjust(prim_win, dual_win, s_adjust);
	}
}




Eigen::VectorXd admm::calc1(const Eigen::VectorXd & sub_z) {
	return sub_z / sigma;
}
double admm::calc2(double tau) {
	return tau * sigma;
}


void admm::atyu2c(concurrency::concurrent_vector<cons_full>&data, int beg, int end, double scale_a, Eigen::VectorXd & x, Eigen::VectorXd & b, Eigen::VectorXd & res) {
	//int64_t tstidx = 0;
	for (int i = beg; i != end; i++) {
		int* idx = data[i].idx;
		//tstidx = (int64_t)idx;
#ifdef XY
		double axbu = data[i].y;
#else
		double axbu = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;

		axbu += data[i].u / sigma;

		if (axbu > 0) {
			axbu = 0;
		}
#endif
		double tem = axbu - (b[idx[0]] - b[idx[1]] - b[idx[2]]) / scale_a;

		tem -= data[i].u / sigma;


		res[idx[0]] += tem / scale_a;
		res[idx[1]] -= tem / scale_a;
		res[idx[2]] -= tem / scale_a;
	}
	//if (tstidx > 0)cout << "hhh\n";


}



void admm::all_in_all2c(concurrency::concurrent_vector<cons_full>::iterator beg, concurrency::concurrent_vector<cons_full>::iterator end, double scale_a, double tau, double& y_norm, double& axby, double& ax, double& bnorm, double& kkt_y, Eigen::VectorXd & x, Eigen::VectorXd & old_x, Eigen::VectorXd & b, Eigen::VectorXd & atu) {

	auto square = [](double x) {return x * x; };
	double _ax = 0, _y_norm = 0, _bnorm = 0, _axby = 0, _kkt_y = 0;
	for (auto i = beg; i != end; i++) {
		{
			double u_value = i->u;
			double axb, axbu;
			int* idx = i->idx;


#ifdef XY
			axbu = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a + u_value / sigma;
			i->y = axbu > 0 ? 0 : axbu;
			axbu = i->y;
			//calc y,ax
			axb = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;

			_y_norm += square(i->y);

			axb -= i->y;

#else
			axbu = (old_x[idx[0]] + b[idx[0]] - old_x[idx[1]] - b[idx[1]] - old_x[idx[2]] - b[idx[2]]) / scale_a + u_value / sigma;
			//calc y,ax
			axb = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;

			if (axbu > 0)
				axbu = 0;
			axb -= axbu;

			_y_norm += square(axbu);
#endif
			_ax += square((x[idx[0]] - x[idx[1]] - x[idx[2]]) / scale_a);
			_bnorm += square((b[idx[0]] - b[idx[1]] - b[idx[2]]) / scale_a);

			_axby += square(axb);

			double u_new = axb * tau * sigma + u_value;
			//if (u_new < 0)u_new = 0;
			i->u = u_new;
			if (u_new + axbu < 0) {
				_kkt_y += square(u_new);
			}
			else {
				_kkt_y += square(axbu);
			}
			//_kkt_y += square(u_new);

			atu[idx[0]] += u_new / scale_a;
			atu[idx[1]] -= u_new / scale_a;
			atu[idx[2]] -= u_new / scale_a;

		}
	}
	ax = _ax; axby = _axby; bnorm = _bnorm, y_norm = _y_norm; kkt_y = _kkt_y;
}


Eigen::VectorXd admm_inf::pi_c(Eigen::VectorXd x) {
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
Eigen::VectorXd admm_inf::pi_c_one(Eigen::VectorXd x) {
	Eigen::VectorXd tmp = x.cwiseAbs();
	double gap = tmp.sum() - 1;

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