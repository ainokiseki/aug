#include"dk.h"
#include<iomanip>
#include<omp.h>
#include<queue>
#include<condition_variable>
#include<atomic>


jk::jk(string graph_name, int type, double g) : metric(graph_name, type)
{
	gam = g;

	y = Eigen::VectorXd::Zero(var_num);
	m = -gam * Eigen::VectorXd::Ones(var_num);
	p = Eigen::VectorXd::Zero(n * n);

	mthread = 30;

}
jk::~jk()
{

}
void jk::multi_triloop2(int threadnum, vector<queue<pid> >& correction, vector<queue<pid> >& next_correction) {

	queue<looplist>q;

	int block = 10;
	int blocknum = (n - 1) / block + 1;
	int sz = correction.size();
	mutex mut;
	condition_variable cond;
	atomic_int ant(0);
	vector<thread>v;
	for (int i = 0; i < threadnum; i++) {
		v.push_back(thread(&jk::worker, this, ref(q), ref(cond), ref(ant), ref(mut), ref(correction), ref(next_correction)));
	}

	//dk_mutex dkloop(n, &y, &m, &inv_w, &mat_b,correction,next_correction);
	vector<std::pair<vector<std::pair<int,int>>,vector<bool>> >ts;
	ts.push_back(std::pair<vector<std::pair<int, int>>, vector<bool>> (vector<std::pair<int,int>>(0),vector<bool>(blocknum,false)));


	for (int i = 0; i < blocknum; i++) {
		for (int j = blocknum - 1; j >= i;j--) {

			bool flg = false;
			for (auto & co : ts) {
				if (co.second[i] || co.second[j])
					continue;
				co.second[i] = co.second[j] = true;
				co.first.push_back(std::pair<int, int>(i, j));
				flg = true;
					break;
			}

			if (!flg) {
				ts.push_back(std::pair<vector<std::pair<int, int>>, vector<bool>>(vector<std::pair<int, int>>(0), vector<bool>(blocknum, false)));
				ts[ts.size()-1].first.push_back(std::pair<int, int>(i, j));
				ts[ts.size() - 1].second[i]= ts[ts.size() - 1].second[j]=true;
			}

		}
	}
	



	for (auto co : ts) {
		int tnum = 0;
		ant = 0;
		for (auto target : co.first) {
			//cout << target.first << " " << target.second << endl;
			int kend= std::fmin((target.second + 1) * block, n);
			std::unique_lock<std::mutex> lk(mut);
			q.push(looplist(target.first * block, (target.first + 1) * block, target.second * block, kend,tnum%sz ));
			lk.unlock();
			cond.notify_one();
			tnum++;
		}
		while (ant != tnum) {}
		//cout << "==========\n";
	}
	//return;
	//cout << "wtf\n";
	//return;


	//cout << "wtf\n";
	{

		std::unique_lock<std::mutex> lk(mut);
		q.push(looplist(-1, -1, -1, -1, -1));
		lk.unlock();
		cond.notify_one();
		
	}


	for (int i = 0; i < threadnum; i++) {
		v[i].join();
	}
}
void jk::multi_triloop(int threadnum, vector<queue<pid> >& correction, vector<queue<pid> >& next_correction) {

	queue<looplist>q;

	int block = 10;
	int blocknum = (n - 1) / block + 1;
	int sz = correction.size();
	mutex mut;
	condition_variable cond;
	atomic_int ant(0);
	vector<thread>v;
	for (int i = 0; i < threadnum; i++) {
		v.push_back(thread(&jk::worker, this, ref(q), ref(cond), ref(ant),  ref(mut), ref(correction), ref(next_correction)));
	}

	//dk_mutex dkloop(n, &y, &m, &inv_w, &mat_b,correction,next_correction);
	for (int i = 0; i < blocknum; i++) {
		ant = 0;
		int istart, iend, kstart, kend;
		istart = i;
		while (istart < (i + blocknum + 1) / 2) {

			kstart = blocknum + i - 1 - istart;
			kend = std::fmin((kstart + 1) * block, n);
			std::unique_lock<std::mutex> lk(mut);
			q.push(looplist(istart * block, (istart + 1) * block, kstart * block, kend, istart % sz));
			//cout << istart * block << "\t" << (istart + 1) * block << "\t" << kstart  * block << "\t" << kend << "\t" << istart % sz << endl;
			lk.unlock();
			cond.notify_one();
			istart++;
		}

		//cout << "ist:"<<istart <<"i:"<<i<< endl;
		while (ant != istart - i && istart != i) {}

	}
	//cout << "wtf\n";
	for (int i = 1; i < blocknum; i++) {
		ant = 0;
		int istart, iend, kstart, kend;
		istart = i;
		while (istart < (i + blocknum + 1) / 2) {

			kstart = blocknum + i - 1 - istart;
			kend = std::fmin((kstart - i + 1) * block, n);
			std::unique_lock<std::mutex> lk(mut);
			q.push(looplist((istart - i) * block, (istart - i + 1) * block, (kstart - i) * block, kend, istart % sz));
			//cout << (istart - i) * block << "\t" << (istart - i + 1) * block << "\t" << (kstart - i) * block << "\t" << kend << "\t" << istart % sz << endl;
			lk.unlock();
			cond.notify_one();

			istart++;
		}
		//cout << "ist:" << istart << "i:" << i << endl;
		while (ant != istart - i && istart != i) {}

	}


	//cout << "wtf\n";
	{

		std::unique_lock<std::mutex> lk(mut);
		q.push(looplist(-1, -1, -1, -1, -1));
		cond.notify_one();
		lk.unlock();
	}
	for (int i = 0; i < threadnum; i++) {
		v[i].join();
	}
}
void jk::worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, mutex& mut, vector<queue<pid> >& vc, vector<queue<pid> >& vnc)
{
	while (true) {



		std::unique_lock<std::mutex> lk(mut);
		cond.wait(lk, [&] {return !q.empty(); });
		auto x = q.front();
		if (x.cnum == -1) {
			lk.unlock();
			cond.notify_one();
			return;
		}
		q.pop();
		lk.unlock();



		int istart, iend;
		istart = x.istart;
		iend = x.iend;
		int kstart, kend;
		kstart = x.kstart;
		kend = x.kend;

		queue<pid>& correction = vc[x.cnum];
		queue<pid>& next_correction = vnc[x.cnum];



		int idx[3];
		int uindex = 0;
		for (int i = istart; i < iend; i++) {


			for (int k = kstart; k < kend; k++) {
				for (int j = i + 1; j < k; j++) {


					//int64_t ijk = static_cast<int64_t>(i) * n * n + j * n + k;
					//int64_t ijk2 = static_cast<int64_t>(i) * n * n + k * n + j;
					//int64_t ijk3 = static_cast<int64_t>(j) * n * n + k * n + i;

					//next_correction.push(pid(ijk, 0));
					//next_correction.push(pid(ijk2, 0));
					//next_correction.push(pid(ijk3, 0));
					idx[0] = ij(i, j, n);


					double wij, wik, wjk;
					double bij, bik, bjk;
					wij = inv_w[idx[0]];
					bij = b[idx[0]];
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					wik = inv_w[idx[1]];
					bik = b[idx[1]];
					wjk = inv_w[idx[2]];
					bjk = b[idx[2]];

					double fm = wij + wik + wjk;
					double mu;
					int64_t ijk;
					ijk = static_cast<int64_t>(i)*n*n+j*n+k;
					double theta;






					if (!correction.empty() && correction.front().first == ijk) {
						theta = correction.front().second;
						correction.pop();

						y[idx[0]] += theta * wij / fm;
						y[idx[1]] -= theta * wik / fm;
						y[idx[2]] -= theta * wjk / fm;
					}
					mu = y[idx[0]] + bij - (y[idx[1]] + bik) - (y[idx[2]] + bjk);
					if (mu > 0) {
						y[idx[0]] -= mu * wij / fm;
						y[idx[1]] += mu * wik / fm;
						y[idx[2]] += mu * wjk / fm;
						next_correction.push(std::pair<int64_t, double>(ijk, mu));

					}






					ijk = static_cast<int64_t>(i) * n * n + k * n + j;



					if (!correction.empty() && correction.front().first == ijk) {
						theta = correction.front().second;
						correction.pop();
						y[idx[1]] += theta * wik / fm;
						y[idx[0]] -= theta * wij / fm;
						y[idx[2]] -= theta * wjk / fm;
					}
					mu = y[idx[1]] + bik - (y[idx[0]] + bij) - (y[idx[2]] + bjk);
					if (mu > 0) {
						y[idx[1]] -= mu * wik / fm;
						y[idx[0]] += mu * wij / fm;
						y[idx[2]] += mu * wjk / fm;
						next_correction.push(std::pair<int64_t, double>(ijk, mu));


					}

					ijk = static_cast<int64_t>(j) * n * n + k * n + i;



					if (!correction.empty() && correction.front().first == ijk) {
						theta = correction.front().second;
						correction.pop();
						y[idx[2]] += theta * wjk / fm;
						y[idx[1]] -= theta * wik / fm;
						y[idx[0]] -= theta * wij / fm;
					}
					mu = y[idx[2]] + bjk - (y[idx[1]] + bik) - (y[idx[0]] + bij);
					if (mu > 0) {
						y[idx[2]] -= mu * wjk / fm;
						y[idx[1]] += mu * wik / fm;
						y[idx[0]] += mu * wij / fm;
						next_correction.push(std::pair<int64_t, double>(ijk, mu));


					}



				}




			}
		}

		ant++;

	}

}


void jk::triloop(int threadnum, std::queue<pid>& correction, std::queue<pid>& next_correction, mutex* mut)
{

	int istart, iend;

	istart = n - (int)(cbrt((threadnum + 1) * 1.0 / mthread) * n);
	iend = n - (int)(cbrt(threadnum * 1.0 / mthread) * 1.0 * n);
	int idx[3];
	int uindex = 0;
	for (int i = istart; i < iend; i++) {

		int64_t ijk1 = static_cast<int64_t>(i) * n * n;
		for (int j = i + 1; j < n; j++) {
			idx[0] = ij(i, j, n);
			int64_t ijk2, ijk3, ijk4;
			ijk2 = ijk1 + j * n;
			ijk3 = static_cast<int64_t>(j) * n * n + i;
			ijk4 = ijk1 + j;
			double wij, wik, wjk;
			double bij, bik, bjk;


			wij = inv_w[idx[0]];
			bij = b[idx[0]];
			for (int k = j + 1; k < n; k++) {
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				wik = inv_w[idx[1]];
				bik = b[idx[1]];
				wjk = inv_w[idx[2]];
				bjk = b[idx[2]];

				double fm = wij + wik + wjk;
				double mu;
				int64_t ijk;
				ijk = ijk2 + k;
				double theta;



				mut[i * n + j].lock();
				mut[i * n + k].lock();
				mut[j * n + k].lock();



				if (!correction.empty() && correction.front().first == ijk) {
					theta = correction.front().second;
					correction.pop();

					y[idx[0]] += theta * wij / fm;
					y[idx[1]] -= theta * wik / fm;
					y[idx[2]] -= theta * wjk / fm;
				}
				mu = y[idx[0]] + bij - (y[idx[1]] + bik) - (y[idx[2]] + bjk);
				if (mu > 0) {
					y[idx[0]] -= mu * wij / fm;
					y[idx[1]] += mu * wik / fm;
					y[idx[2]] += mu * wjk / fm;
					next_correction.push(std::pair<int64_t, double>(ijk, mu));

				}






				ijk = ijk4 + k * n;



				if (!correction.empty() && correction.front().first == ijk) {
					theta = correction.front().second;
					correction.pop();
					y[idx[1]] += theta * wik / fm;
					y[idx[0]] -= theta * wij / fm;
					y[idx[2]] -= theta * wjk / fm;
				}
				mu = y[idx[1]] + bik - (y[idx[0]] + bij) - (y[idx[2]] + bjk);
				if (mu > 0) {
					y[idx[1]] -= mu * wik / fm;
					y[idx[0]] += mu * wij / fm;
					y[idx[2]] += mu * wjk / fm;
					next_correction.push(std::pair<int64_t, double>(ijk, mu));


				}

				ijk = ijk3 + k * n;



				if (!correction.empty() && correction.front().first == ijk) {
					theta = correction.front().second;
					correction.pop();
					y[idx[2]] += theta * wjk / fm;
					y[idx[1]] -= theta * wik / fm;
					y[idx[0]] -= theta * wij / fm;
				}
				mu = y[idx[2]] + bjk - (y[idx[1]] + bik) - (y[idx[0]] + bij);
				if (mu > 0) {
					y[idx[2]] -= mu * wjk / fm;
					y[idx[1]] += mu * wik / fm;
					y[idx[0]] += mu * wij / fm;
					next_correction.push(std::pair<int64_t, double>(ijk, mu));


				}



				mut[j * n + k].unlock();

				mut[i * n + k].unlock();
				mut[i * n + j].unlock();
			}




		}
	}



}

void jk::biloop()
{


	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			int idx = ij(i, j, n);
			int idx1, idx2;
			idx1 = i * n + j;
			idx2 = j * n + i;
			double mu;

			y[idx] -= p[idx1];
			m[idx] -= p[idx1];
			mu = -y[idx] - m[idx];
			if (mu > 0) {
				y[idx] += mu / 2;
				m[idx] += mu / 2;
				p[idx1] = mu / 2;
			}
			else {
				p[idx1] = 0;
			}
			y[idx] += p[idx2];
			m[idx] -= p[idx2];
			mu = y[idx] - m[idx];
			if (mu > 0) {
				y[idx] -= mu / 2;
				m[idx] += mu / 2;
				p[idx2] = mu / 2;
			}
			else {
				p[idx2] = 0;
			}



		}
	}




}

void jk::multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction)
{



	vector<thread>v;
	mutex* mut = new mutex[n * n];

	for (int i = 0; i < mthread; i++) {
		v.push_back(thread(&jk::triloop, this, i, std::ref(correction[i]), std::ref(next_correction[i]), mut));
	}
	for (int i = 0; i < mthread; i++) {
		v[i].join();
	}
	delete[]mut;
}
void jk::multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction, int te) {
	int threadnum = mthread;
	int block = 10;
	int blocknum = (n - 1) / block + 1;

	//dk_mutex dkloop(n, &y, &m, &inv_w, &mat_b,correction,next_correction);
	vector<thread>v;
	for (int t = 0; t < blocknum; t++) {
		//Æ½¾ù·ÖÅäblock
		int tnum = std::min(threadnum, (blocknum - t + 1) / 2);
		//cout << "tnum:" << tnum << endl;

		v.clear();
		for (int th = 0; th < tnum; th++) {
			int istart, iend, kstart;
			iend = (blocknum + 1 - t) / 2 - int(cbrt(double(tnum - 1 - th) / tnum) * ((blocknum - t + 1) / 2));
			istart = (blocknum + 1 - t) / 2 - int(cbrt(double(tnum - th) / tnum) * ((blocknum - t + 1) / 2));
			kstart = std::min(blocknum - t - 1 - (istart), blocknum - t - 1 - (iend - 1));
			if (th == tnum - 1) {
				iend = (blocknum + 1 - t) / 2;
				kstart = (blocknum - t) / 2;
			}

			if (iend <= istart)continue;
			cout << istart << " " << iend << " " << kstart << endl;
			if (istart == kstart)cout << "--------------------\n";
		}

		if (t > 0) {


			v.clear();
			for (int th = 0; th < tnum; th++) {
				int istart, iend, kstart;
				iend = (blocknum + 1 - t) / 2 - int(cbrt(double(tnum - 1 - th) / tnum) * ((blocknum - t + 1) / 2));
				istart = (blocknum + 1 - t) / 2 - int(cbrt(double(tnum - th) / tnum) * ((blocknum - t + 1) / 2));
				kstart = std::min(blocknum - t - 1 - (istart), blocknum - t - 1 - (iend - 1));
				if (th == tnum - 1) {
					iend = (blocknum + 1 - t) / 2;
					kstart = (blocknum - t) / 2;
				}

				if (iend <= istart)continue;
				cout << istart + t << " " << iend + t << " " << kstart + t << endl;
				//v.push_back(std::thread(&jk::triloop2,this, istart + t, iend + t, kstart + t, block, std::ref(correction[th]), std::ref(next_correction[th])));
			}

		}

	}
}
bool jk::mconstrain_check(double tol)
{
	atomic_int ant(0);

	int pnum = mthread;
	vector<thread>v;
	for (int i = 0; i < pnum; i++) {
		
		v.push_back(thread(&jk::constrain_check,this, tol,ref(ant),n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
	}

	for (int i = 0; i < pnum; i++) {
		v[i].join();
	}
	if (ant > 0)return false;
	return true;
}
void jk::constrain_check(double tol,atomic_int& ant,int istart,int iend) {
	for (int i = istart; i < iend; i++) {

		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				if (ant > 0)return;

				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				if (y[idx[0]] + b[idx[0]] - (y[idx[1]] + b[idx[1]]) - (y[idx[2]] + b[idx[2]]) > tol ||
					y[idx[1]] + b[idx[1]] - (y[idx[0]] + b[idx[0]]) - (y[idx[2]] + b[idx[2]]) > tol ||
					y[idx[2]] + b[idx[2]] - (y[idx[1]] + b[idx[1]]) - (y[idx[0]] + b[idx[0]]) > tol
					) {
					ant++;
					return;
				}
			}
		}


	}



}



double jk::run(double tol)
{
	inv_w = Eigen::VectorXd::Ones(var_num).cwiseQuotient(w);
	double pobj = 0.0;
	double lastobj = 0.0;



	std::vector<std::queue<pid> >correction(500);
	std::vector<std::queue<pid> >next_correction(500);
	auto& c1 = correction;
	auto& c2 = next_correction;

	kkt_1 = kkt_2 = 0.0;
	int iter_I;
	for (iter_I = 0; iter_I < 10000; iter_I++) {
		kkt_2 = kkt_1 = 0.0;
		//triloop(c1,c2);





		multi_triloop2(mthread, c1, c2);




		//return 0;
		//triloop(0, c1[0], c2[0], NULL);
		//return 0;
		//std::map<int64_t, int>tst;
		//for (int i = 0; i < n; i++) {
		//	for (int j = i + 1; j < n; j++) {
		//		for (int k = j + 1; k < n; k++) {
		//			tst[static_cast<int64_t>(i) * n * n + j * n + k] = 0;
		//			//tst[static_cast<int64_t>(i) * n * n + k* n + j] = 0;
		//			//tst[static_cast<int64_t>(j) * n * n + k * n + i] = 0;
		//		}
		//	}
		//}


		//for (int co = 0; co < c2.size(); co++) {
		//	while (!c2[co].empty()) {

		//		int64_t longindex = c2[co].front().first;
		//		
		//		c2[co].pop();
		//		tst[longindex]++;
		//		
		//	}
		//}
		//for (auto iter : tst) {
		//	if (iter.second != 1) {

		//		int64_t longindex = iter.first;
		//		int i, j, k;
		//		k = longindex % n;
		//		j = longindex / n;
		//		i = j / n;
		//		j = j % n;
		//		cout << "i:" << i << " j:" << j << " k:" << k << " num:" << iter.second << endl;

		//	}
		//}
		////
		//return 0.0;


		//return 0.0;
		biloop();

		int tst = 0;
		int ttt = 0;
		for (auto i : c1) {
			tst += i.size();
		}
		for (auto i : c2) {
			ttt += i.size();
		}
		cout << tst << "," << ttt << endl;



		std::swap(c1, c2);
		pobj = 0.0;
		int index = 0;


		//cout << "====" << endl;


		if ((iter_I + 1) % 10 == 0) {

			cout << "=============================\n";
			pobj = m.cwiseProduct(w).sum();
			double wnorm = m.cwiseProduct(m).cwiseProduct(w).sum() + y.cwiseProduct(y).cwiseProduct(w).sum();

			lastobj = pobj;

			double dobj = 0;
#pragma omp parallel for reduction(+:dobj)

			for (int co = 0; co < c1.size(); co++) {
				while (!c1[co].empty()) {

					int64_t longindex = c1[co].front().first;
					double mu = c1[co].front().second;
					c2[co].push(c1[co].front());
					c1[co].pop();
					int i, j, k;
					k = longindex % n;
					j = longindex / n;
					i = j / n;
					j = j % n;
					int idx[3];
					idx[0] = ij(i, j, n);
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					double fm = inv_w[idx[0]] + inv_w[idx[1]] + inv_w[idx[2]];
					dobj += mu * (b[idx[0]] - b[idx[1]] - b[idx[2]]) / fm;
				}
			}
			std::swap(c1, c2);
			dobj /= gam;
			dobj -= wnorm / (2 * gam);
			pobj += wnorm / (2 * gam);
			std::cout << iter_I << " " << std::setprecision(8) << pobj << "\tdobj:" << dobj << "\twnorm:" << wnorm << std::endl;
			if (std::fabs((pobj - dobj) / dobj) < 1e-2&&mconstrain_check(1e-2))break;
		}

	}
	double trueobj = 0.0;
	auto square = [](double x) {return x * x; };
	int index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++, index++) {
			trueobj += square(y[ij(i, j, n)]) * w[index];
			trueobj += square(m[ij(i, j, n)]) * w[index];

		}
	}

	int ct = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			x[ct] = y[ij(i, j, n)];
			ct++;
		}
	}
	std::cout << "final res:" << trueobj / (2 * gam) + pobj << "k1:" << kkt_1 << "k2:" << kkt_2 << std::endl;
	std::cout << "gam:" << gam << "\tobj:" << lastobj << endl;




	return 0.0;
}

void jk::init()
{
}

double jk::run(double tol)
{
	inv_w = Eigen::VectorXd::Ones(var_num).cwiseQuotient(w);
	double pobj = 0.0;
	double lastobj = 0.0;



	std::vector<std::queue<pid> >correction(500);
	std::vector<std::queue<pid> >next_correction(500);
	auto& c1 = correction;
	auto& c2 = next_correction;

	kkt_1 = kkt_2 = 0.0;
	int iter_I;
	for (iter_I = 0; iter_I < 10000; iter_I++) {
		kkt_2 = kkt_1 = 0.0;
		//triloop(c1,c2);





		multi_triloop2(mthread, c1, c2);




		//return 0;
		//triloop(0, c1[0], c2[0], NULL);
		//return 0;
		//std::map<int64_t, int>tst;
		//for (int i = 0; i < n; i++) {
		//	for (int j = i + 1; j < n; j++) {
		//		for (int k = j + 1; k < n; k++) {
		//			tst[static_cast<int64_t>(i) * n * n + j * n + k] = 0;
		//			//tst[static_cast<int64_t>(i) * n * n + k* n + j] = 0;
		//			//tst[static_cast<int64_t>(j) * n * n + k * n + i] = 0;
		//		}
		//	}
		//}


		//for (int co = 0; co < c2.size(); co++) {
		//	while (!c2[co].empty()) {

		//		int64_t longindex = c2[co].front().first;
		//		
		//		c2[co].pop();
		//		tst[longindex]++;
		//		
		//	}
		//}
		//for (auto iter : tst) {
		//	if (iter.second != 1) {

		//		int64_t longindex = iter.first;
		//		int i, j, k;
		//		k = longindex % n;
		//		j = longindex / n;
		//		i = j / n;
		//		j = j % n;
		//		cout << "i:" << i << " j:" << j << " k:" << k << " num:" << iter.second << endl;

		//	}
		//}
		////
		//return 0.0;


		//return 0.0;
		biloop();

		int tst = 0;
		int ttt = 0;
		for (auto i : c1) {
			tst += i.size();
		}
		for (auto i : c2) {
			ttt += i.size();
		}
		cout << tst << "," << ttt << endl;



		std::swap(c1, c2);
		pobj = 0.0;
		int index = 0;


		//cout << "====" << endl;


		if ((iter_I + 1) % 10 == 0) {

			cout << "=============================\n";
			pobj = m.cwiseProduct(w).sum();
			double wnorm = m.cwiseProduct(m).cwiseProduct(w).sum() + y.cwiseProduct(y).cwiseProduct(w).sum();

			lastobj = pobj;

			double dobj = 0;
#pragma omp parallel for reduction(+:dobj)

			for (int co = 0; co < c1.size(); co++) {
				while (!c1[co].empty()) {

					int64_t longindex = c1[co].front().first;
					double mu = c1[co].front().second;
					c2[co].push(c1[co].front());
					c1[co].pop();
					int i, j, k;
					k = longindex % n;
					j = longindex / n;
					i = j / n;
					j = j % n;
					int idx[3];
					idx[0] = ij(i, j, n);
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					double fm = inv_w[idx[0]] + inv_w[idx[1]] + inv_w[idx[2]];
					dobj += mu * (b[idx[0]] - b[idx[1]] - b[idx[2]]) / fm;
				}
			}
			std::swap(c1, c2);
			dobj /= gam;
			dobj -= wnorm / (2 * gam);
			pobj += wnorm / (2 * gam);
			std::cout << iter_I << " " << std::setprecision(8) << pobj << "\tdobj:" << dobj << "\twnorm:" << wnorm << std::endl;
			if (std::fabs((pobj - dobj) / dobj) < 1e-2 && mconstrain_check(1e-2))break;
		}

	}
	double trueobj = 0.0;
	auto square = [](double x) {return x * x; };
	int index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++, index++) {
			trueobj += square(y[ij(i, j, n)]) * w[index];
			trueobj += square(m[ij(i, j, n)]) * w[index];

		}
	}

	int ct = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			x[ct] = y[ij(i, j, n)];
			ct++;
		}
	}
	std::cout << "final res:" << trueobj / (2 * gam) + pobj << "k1:" << kkt_1 << "k2:" << kkt_2 << std::endl;
	std::cout << "gam:" << gam << "\tobj:" << lastobj << endl;




	return 0.0;
}