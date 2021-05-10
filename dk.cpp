#include"dk.h"
#include<iomanip>
dk::dk(string graph_name) :metric(graph_name)
{
	gam =1;
}
dk::dk(string graph_name,double g) : metric(graph_name)
{
	gam = g;
}
dk::~dk()
{
}
inline long long ij(long long x, long long y, long long n) {
	long long res = n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;

	return res;
}
void dk::triloop(std::map<int64_t,double>*correction, std::map<int64_t, double>* next_correction)
{

	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				double fm = inv_w(i, j) + inv_w(i, k) + inv_w(j, k);
				double mu;
				int64_t ijk;
				ijk = static_cast<int64_t>(i) * n * n + j * n + k;
				if (correction->find(ijk) != correction->end()) {
					y(i, j) += (*correction)[ijk] * inv_w(i, j) / fm;
					y(i, k) -= (*correction)[ijk] * inv_w(i, k) / fm;
					y(j, k) -= (*correction)[ijk] * inv_w(j, k) / fm;
				}
				mu = y(i, j) + mat_b(i, j) - (y(i, k) + mat_b(i, k)) - (y(j, k) + mat_b(j, k));
				if (mu > 0) {
					y(i, j) -= mu * inv_w(i, j) / fm;
					y(i, k) += mu * inv_w(i, k) / fm;
					y(j, k) += mu * inv_w(j, k) / fm;
					next_correction->insert(std::pair<int64_t, double>(ijk, mu));
				}
				ijk = static_cast<int64_t>(i) * n * n + k * n + j;
				if (correction->find(ijk) != correction->end()) {
					y(i, k) += (*correction)[ijk] * inv_w(i, k) / fm;
					y(i, j) -= (*correction)[ijk] * inv_w(i, j) / fm;
					y(j, k) -= (*correction)[ijk] * inv_w(j, k) / fm;
				}
				mu = y(i, k) + mat_b(i, k) - (y(i, j) + mat_b(i, j)) - (y(j, k) + mat_b(j, k));
				if (mu > 0) {
					y(i, k) -= mu * inv_w(i, k) / fm;
					y(i, j) += mu * inv_w(i, j) / fm;
					y(j, k) += mu * inv_w(j, k) / fm;
					next_correction->insert(std::pair<int64_t, double>(ijk, mu));

				}
				ijk = static_cast<int64_t>(j) * n * n + i * n + k;
				if (correction->find(ijk) != correction->end()) {
					y(j, k) += (*correction)[ijk] * inv_w(j, k) / fm;
					y(i, k) -= (*correction)[ijk] * inv_w(i, k) / fm;
					y(i, j) -= (*correction)[ijk] * inv_w(i, j) / fm;
				}
				mu = y(j, k) + mat_b(j, k) - (y(i, k) + mat_b(i, k)) - (y(i, j) + mat_b(i, j));
				if (mu > 0) {
					y(j, k) -= mu * inv_w(j, k) / fm;
					y(i, k) += mu * inv_w(i, k) / fm;
					y(i, j) += mu * inv_w(i, j) / fm;
					next_correction->insert(std::pair<int64_t, double>(ijk, mu));

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
				y(i, j) += mu / 2;
				m(i, j) += mu / 2;
				p(i, j) = mu / 2;
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
void dk::multi_triloop(std::map<int64_t, double>* correction, std::map<int64_t, double>* next_correction) {
	int threadnum = 12;
	int block = 10;
	int blocknum = (n-1) / block+1;
	dk_triloop dkloop(n, &y, &m, &inv_w, &mat_b,correction,next_correction);
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
			if (th == tnum-1) {
				iend = (blocknum + 1 - t) / 2;
				kstart = (blocknum - t) / 2;
			}

			if (iend <= istart)continue;
			v.push_back(std::thread(std::ref(dkloop), istart, iend, kstart, block));
		}
		for (int i = 0; i <v.size(); i++) {
			v[i].join();
		}
		if (t > 0) {


			v.clear();
			for (int th = 0; th < tnum; th++) {
				int istart, iend, kstart;
				iend = (blocknum + 1 - t) / 2 - int(cbrt(double(tnum - 1 - th) / tnum) * ((blocknum - t + 1) / 2));
				istart = (blocknum + 1 - t) / 2 - int(cbrt(double(tnum - th) / tnum) * ((blocknum - t + 1) / 2));
				kstart = std::min(blocknum - t - 1 - (istart), blocknum - t - 1 - (iend - 1));
				if (th == tnum-1) {
					iend = (blocknum + 1 - t) / 2;
					kstart = (blocknum - t) / 2;
				}

				if (iend <= istart)continue;
				v.push_back(std::thread(std::ref(dkloop), istart + t, iend + t, kstart + t, block));
			}
			for (int i = 0; i < v.size(); i++) {
				v[i].join();
			}
		}

	}
}

void get_kkt(std::map<int64_t, double>* c1, std::map<int64_t, double>* c2, double* kkt_1,double* kkt_2, dk* th ,int start = 0,int end=-1) {
	if (end < 0)end = th->n;
	for (int ii = start; ii < end; ii++) {
		for (int jj = ii + 1; jj < th->n; jj++) {
			for (int kk = jj + 1; kk < th->n; kk++) {
				int64_t ijk = static_cast<int64_t>(ii) * th->n * th->n + jj * th->n + kk;
				double ty;
				double fm = th->inv_w(ii, jj) + th->inv_w(ii, kk) + th->inv_w(jj, kk);
				double mu = th->y(ii, jj) + th->mat_b(ii, jj) - (th->y(ii, kk) + th->mat_b(ii, kk)) - (th->y(jj, kk) + th->mat_b(jj, kk));

				if (c1->find(ijk) == c1->end()) {
					ty = 0;
				}
				else {
					ty = (*c1)[ijk] / (fm * th->gam);
				}
				if (mu + ty > 0) {
					(*kkt_1) += mu * mu;
				}
				(*kkt_2) += mu * mu;

				ijk = static_cast<int64_t>(ii) * th->n * th->n + kk * th->n + jj;
				mu = th->y(ii, kk) + th->mat_b(ii, kk) - (th->y(ii, jj) + th->mat_b(ii, jj)) - (th->y(jj, kk) + th->mat_b(jj, kk));
				if (c1->find(ijk) == c1->end()) {
					ty = 0;
				}
				else {
					ty = (*c1)[ijk] / (fm * th->gam);
				}
				if (mu + ty > 0) {
					(*kkt_1) += mu * mu;
				}
				(*kkt_2) += mu * mu;

				ijk = static_cast<int64_t>(jj) * th->n * th->n + ii * th->n + kk;
				mu = th->y(jj, kk) + th->mat_b(jj, kk) - (th->y(ii, kk) + th->mat_b(ii, kk)) - (th->y(ii, jj) + th->mat_b(ii, jj));
				if (c1->find(ijk) == c1->end()) {
					ty = 0;
				}
				else {
					ty = (*c1)[ijk] / (fm * th->gam);
				}
				if (mu + ty > 0) {
					(*kkt_1) += mu * mu;
				}
				(*kkt_2) += mu * mu;
			}
		}
	}

}


int dk::run_test(double threshod)
{
	init();
	auto i_start = std::chrono::high_resolution_clock::now();
	double pobj = 0.0;
	double lastobj = 0.0;
	std::map<int64_t, double>correction;
	std::map<int64_t, double>next_correction;
	std::map<int64_t, double>* c1;
	std::map<int64_t, double>* c2;
	c1 = &correction;
	c2 = &next_correction;
	kkt_1 = kkt_2 = 0.0;
	int k;
	for (k = 0; k < 1000; k++) {
		kkt_2 = kkt_1 = 0.0;
		//triloop(c1,c2);
		multi_triloop(c1, c2);
		biloop();
		c1->clear();
		std::swap(c1, c2);
		pobj = 0.0;
		int index = 0;

		pobj = Eigen::MatrixXd(m.cwiseProduct(inv_w.transpose()).triangularView<Eigen::StrictlyUpper>()).sum();


		kkt_1 = kkt_2 = 0.0;
		double kkt_3 = 0;
		if (k % 20 == 0) {
			double kkt_1 = 0;
			double kkt_2 = 0;
			vector<thread>v;
			int pnum = 8;


			double* mphi = new double[pnum];
			double* mphi2 = new double[pnum];
			for (int i = 0; i < pnum; i++)
				mphi[i] =mphi2[i]= 0;
			vector<thread>tv;
			for (int i = 0; i < pnum; i++) {
				int a, b;
				a = n - (int)(cbrt((i + 1) * 1.0 / pnum) * n);
				b = n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n);
				tv.push_back(std::thread(get_kkt,c1,c2,mphi+i,mphi2+i,this, a, b));

			}
			for (int i = 0; i < pnum; i++)
				tv[i].join();
			for (int i = 0; i < pnum; i++)
				kkt_1 += mphi[i];
			for (int i = 0; i < pnum; i++)
				kkt_2 += mphi2[i];
			delete[]mphi;
			delete[]mphi2;


			//for (int ii = 0; ii < n; ii++) {
			//	for (int jj = ii + 1; jj < n; jj++) {
			//		for (int kk = jj + 1; kk < n; kk++) {
			//			int64_t ijk = static_cast<int64_t>(ii) * n * n + jj * n + kk;
			//			double ty;
			//			double fm = inv_w(ii, jj) + inv_w(ii, kk) + inv_w(jj, kk);
			//			double mu = y(ii, jj) + mat_b(ii, jj) - (y(ii, kk) + mat_b(ii, kk)) - (y(jj, kk) + mat_b(jj, kk));
			//			kkt_2 += mu * mu;
			//			if (c1->find(ijk) == c1->end()) {
			//				ty = 0;
			//			}
			//			else {
			//				ty = (*c1)[ijk] / (fm * gam);
			//			}

			//			kkt_3 += ty * mu;

			//			ijk = static_cast<int64_t>(ii) * n * n + kk * n + jj;
			//			mu = y(ii, kk) + mat_b(ii, kk) - (y(ii, jj) + mat_b(ii, jj)) - (y(jj, kk) + mat_b(jj, kk));
			//			if (c1->find(ijk) == c1->end()) {
			//				ty = 0;
			//			}
			//			else {
			//				ty = (*c1)[ijk] / (fm * gam);
			//			}
			//			if (mu > 0) {
			//				kkt_1 += mu * mu;
			//			}
			//			else {
			//				kkt_1 += ty * ty;
			//			}
			//			kkt_3 += ty * mu;
			//			kkt_2 += mu * mu;
			//			ijk = static_cast<int64_t>(jj) * n * n + ii * n + kk;
			//			mu = y(jj, kk) + mat_b(jj, kk) - (y(ii, kk) + mat_b(ii, kk)) - (y(ii, jj) + mat_b(ii, jj));
			//			if (c1->find(ijk) == c1->end()) {
			//				ty = 0;
			//			}
			//			else {
			//				ty = (*c1)[ijk] / (fm * gam);
			//			}
			//			if (mu > 0) {
			//				kkt_1 += mu * mu;
			//			}
			//			else {
			//				kkt_1 += ty * ty;
			//			}
			//			kkt_2 += mu * mu;
			//			kkt_3 += ty * mu;
			//		}
			//	}
			//}
			/*
			Eigen::VectorXd ay = Eigen::VectorXd::Zero(var_num);
			for (int i = 0; i < n; i++) {
				for (int j = i + 1; j < n; j++) {
					for (int kk = j + 1; kk < n; kk++) {
						int idx[3];
						double fm = inv_w(i, j) + inv_w(i, kk) + inv_w(j, kk);
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, kk, n);
						idx[2] = ij(j, kk, n);
						int64_t ijk= static_cast<int64_t>(i) * n * n + j * n + kk;

						if (c1->find(ijk) != c1->end()) {
							double tem = (*c1)[ijk]/fm;
							ay[idx[0]] += tem;
							ay[idx[1]] -= tem;
							ay[idx[2]] -= tem;
						}
						 ijk = static_cast<int64_t>(i) * n * n + kk * n + j;

						if (c1->find(ijk) != c1->end()) {
							double tem = (*c1)[ijk] / fm;
							ay[idx[1]] += tem;
							ay[idx[0]] -= tem;
							ay[idx[2]] -= tem;
						}
						 ijk = static_cast<int64_t>(kk) * n * n + i * n + j;

						if (c1->find(ijk) != c1->end()) {
							double tem = (*c1)[ijk] / fm;
							ay[idx[2]] += tem;
							ay[idx[1]] -= tem;
							ay[idx[0]] -= tem;
						}
					}
				}

			}
			std::cout << "spkkt:" << ay.norm() << endl;
			*/

			lastobj = pobj;
			kkt_3 = sqrt(kkt_1) / (1 + sqrt(kkt_2));
			std::cout << k << " " << std::setprecision(8) << pobj << "\tkkt_3:" << kkt_3 << std::endl;
			if (std::fabs(kkt_3) < threshod)break;
		}

	}
	double trueobj = 0.0;
	auto square = [](double x) {return x * x; };
	int index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++, index++) {
			trueobj += square(y(i, j)) * w[index];
			trueobj += square(m(i, j)) * w[index];

		}
	}
	std::cout << "final res:" << trueobj / (2 * gam) + pobj << "k1:" << kkt_1 << "k2:" << kkt_2 << std::endl;
	std::cout << "gam:" << gam << "\tobj:"<< lastobj << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	cout << "time is:" << tm << endl;
	return k;
}
void dk::run()
{
	init();
	double pobj = 0.0;
	double lastobj = 0.0;
	std::map<int64_t, double>correction;
	std::map<int64_t, double>next_correction;
	std::map<int64_t, double>* c1;
	std::map<int64_t, double>* c2;
	c1 = &correction;
	c2 = &next_correction;
	kkt_1 = kkt_2 = 0.0;
	for (int k = 0; k < 1000; k++) {
		kkt_2 = kkt_1=0.0;
		//triloop(c1,c2);
		multi_triloop(c1, c2);
		biloop();
		c1->clear();
		std::swap(c1, c2);
		pobj = 0.0;
		int index = 0;

		pobj =Eigen::MatrixXd( m.cwiseProduct(inv_w.transpose()).triangularView<Eigen::StrictlyUpper>()).sum();


		kkt_1 = kkt_2 = 0.0;
		double kkt_3 = 0;
		if (k % 1 == 0) {


			for (int ii = 0; ii < n; ii++) {
				for (int jj = ii + 1; jj < n; jj++) {
					for (int kk = jj + 1; kk < n; kk++) {
						int64_t ijk = static_cast<int64_t>(ii) * n * n + jj * n + kk;
						double ty;
						double fm = inv_w(ii, jj) + inv_w(ii, kk) + inv_w(jj, kk);
						double mu = y(ii, jj) + mat_b(ii, jj) - (y(ii, kk) + mat_b(ii, kk)) - (y(jj, kk) + mat_b(jj, kk));
						kkt_2 += mu * mu;
						if (c1->find(ijk) == c1->end()) {
							ty = 0;
						}else {
							ty = (*c1)[ijk]/(fm*gam);
						}
						if (mu+ty > 0) {
							kkt_1 += mu * mu;
						}
						else {
							kkt_1 += ty * ty;
						}
						kkt_3 += ty* mu;

						ijk = static_cast<int64_t>(ii) * n * n + kk * n + jj;
						mu = y(ii, kk) + mat_b(ii, kk) - (y(ii, jj) + mat_b(ii, jj)) - (y(jj, kk) + mat_b(jj, kk));
						if (c1->find(ijk) == c1->end()) {
							ty = 0;
						}else {
							ty = (*c1)[ijk] / (fm * gam);
						}
						if (mu + ty > 0) {
							kkt_1 += mu * mu;
						}else {
							kkt_1 += ty * ty;
						}
						kkt_3 += ty * mu;
						kkt_2 += mu * mu;
						ijk = static_cast<int64_t>(jj) * n * n + ii * n + kk;
						mu = y(jj, kk) + mat_b(jj, kk) - (y(ii, kk) + mat_b(ii, kk)) - (y(ii, jj) + mat_b(ii, jj));
						if (c1->find(ijk) == c1->end()) {
							ty = 0;
						}else {
							ty = (*c1)[ijk] / (fm * gam);
						}
						if (mu + ty > 0) {
							kkt_1 += mu * mu;
						}else {
							kkt_1 += ty * ty;
						}
						kkt_2 += mu * mu;
						kkt_3 += ty * mu;
					}
				}
			}
			/*
			Eigen::VectorXd ay = Eigen::VectorXd::Zero(var_num);
			for (int i = 0; i < n; i++) {
				for (int j = i + 1; j < n; j++) {
					for (int kk = j + 1; kk < n; kk++) {
						int idx[3];
						double fm = inv_w(i, j) + inv_w(i, kk) + inv_w(j, kk);
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, kk, n);
						idx[2] = ij(j, kk, n);
						int64_t ijk= static_cast<int64_t>(i) * n * n + j * n + kk;

						if (c1->find(ijk) != c1->end()) {							
							double tem = (*c1)[ijk]/fm;
							ay[idx[0]] += tem;
							ay[idx[1]] -= tem;
							ay[idx[2]] -= tem;
						}
						 ijk = static_cast<int64_t>(i) * n * n + kk * n + j;

						if (c1->find(ijk) != c1->end()) {
							double tem = (*c1)[ijk] / fm;
							ay[idx[1]] += tem;
							ay[idx[0]] -= tem;
							ay[idx[2]] -= tem;
						}
						 ijk = static_cast<int64_t>(kk) * n * n + i * n + j;

						if (c1->find(ijk) != c1->end()) {
							double tem = (*c1)[ijk] / fm;
							ay[idx[2]] += tem;
							ay[idx[1]] -= tem;
							ay[idx[0]] -= tem;
						}
					}
				}

			}
			std::cout << "spkkt:" << ay.norm() << endl;
			*/
			kkt_1 = sqrt(kkt_1) / (1 + sqrt(kkt_2));
			std::cout << k<<" "<<std::setprecision(8)<<pobj << "\t" << lastobj - pobj << "\tkkt1:" << kkt_1<<"\tkkt_3:"<<kkt_3<< std::endl;
			lastobj = pobj;
			if (std::fabs(kkt_1) < 1e-5)break;
		}

	}
	double trueobj = 0.0;
	auto square = [](double x) {return x * x; };
	int index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++,index++) {
			trueobj += square(y(i, j)) * w[index];
			trueobj += square(m(i, j)) * w[index];

		}
	}
	std::cout << "final res:" << trueobj / (2 * gam) + pobj << "k1:"<<kkt_1<<"k2:"<<kkt_2<<std::endl;

}

void dk::init()
{
	metric::init();
	inv_w = Eigen::MatrixXd::Zero(n, n);
	mat_b = Eigen::MatrixXd::Zero(n, n);
	int index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++, index++) {
			inv_w(i, j) = 1 / w[index];
			inv_w(j, i) = w[index];
		}

	}
	index = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++, index++) {
			mat_b(i, j) = b[index];
		}

	}
	y = Eigen::MatrixXd::Zero(n, n);
	m = -gam * Eigen::MatrixXd::Ones(n, n);
	p = Eigen::MatrixXd::Zero(n, n);

}


dk_triloop::dk_triloop(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::map<int64_t, double>* c1, std::map<int64_t, double>* c2) :
	n(n), y(y), m(m), inv_w(inv_w), mat_b(mat_b), correction(c1), next_correction(c2)
{

}

dk_triloop::~dk_triloop()
{
}

void dk_triloop::operator()(int istart, int iend, int kstart, int blocksize)
{

	int kend = kstart + iend - istart;
	for (int iblock = istart; iblock < iend; iblock++) {


		for (int i = iblock * blocksize; i < (iblock + 1) * blocksize && i < n; i++) {
			int true_k_end = kstart + iend - istart;
			int kblock = kend - 1 - (iblock - istart);


			for (int k = std::max(kblock * blocksize, i + 2); k < n && k < (kblock + 1) * blocksize; k++) {
				std::vector<std::pair<int64_t, double>>tv;
				for (int j = i + 1; j < k; j++) {
					double fm = (*inv_w)(i, j) + (*inv_w)(i, k) + (*inv_w)(j, k);
					double mu;
					int64_t ijk;
					ijk = static_cast<int64_t>(i) * n * n + j * n + k;
					bool flg;
					
					flg = correction->find(ijk) != correction->end();
					
					if (flg) {
						(*y)(i, j) += (*correction)[ijk] * (*inv_w)(i, j) / fm;
						(*y)(i, k) -= (*correction)[ijk] * (*inv_w)(i, k) / fm;
						(*y)(j, k) -= (*correction)[ijk] * (*inv_w)(j, k) / fm;
					}
					mu = (*y)(i, j) + (*mat_b)(i, j) - ((*y)(i, k) + (*mat_b)(i, k)) - ((*y)(j, k) + (*mat_b)(j, k));
					if (mu > 0) {
						(*y)(i, j) -= mu * (*inv_w)(i, j) / fm;
						(*y)(i, k) += mu * (*inv_w)(i, k) / fm;
						(*y)(j, k) += mu * (*inv_w)(j, k) / fm;
						tv.push_back(std::pair<int64_t, double>(ijk, mu));

					}
					ijk = static_cast<int64_t>(i) * n * n + k * n + j;
					
					flg = correction->find(ijk) != correction->end();
					
					if (flg) {
						(*y)(i, k) += (*correction)[ijk] * (*inv_w)(i, k) / fm;
						(*y)(i, j) -= (*correction)[ijk] * (*inv_w)(i, j) / fm;
						(*y)(j, k) -= (*correction)[ijk] * (*inv_w)(j, k) / fm;
					}
					mu = (*y)(i, k) + (*mat_b)(i, k) - ((*y)(i, j) + (*mat_b)(i, j)) - ((*y)(j, k) + (*mat_b)(j, k));
					if (mu > 0) {
						(*y)(i, k) -= mu * (*inv_w)(i, k) / fm;
						(*y)(i, j) += mu * (*inv_w)(i, j) / fm;
						(*y)(j, k) += mu * (*inv_w)(j, k) / fm;
						tv.push_back(std::pair<int64_t, double>(ijk, mu));


					}
					ijk = static_cast<int64_t>(j) * n * n + i * n + k;
					
					flg = correction->find(ijk) != correction->end();
					
					if (flg) {
						(*y)(j, k) += (*correction)[ijk] * (*inv_w)(j, k) / fm;
						(*y)(i, k) -= (*correction)[ijk] * (*inv_w)(i, k) / fm;
						(*y)(i, j) -= (*correction)[ijk] * (*inv_w)(i, j) / fm;
					}
					mu = (*y)(j, k) + (*mat_b)(j, k) - ((*y)(i, k) + (*mat_b)(i, k)) - ((*y)(i, j) + (*mat_b)(i, j));
					if (mu > 0) {
						(*y)(j, k) -= mu * (*inv_w)(j, k) / fm;
						(*y)(i, k) += mu * (*inv_w)(i, k) / fm;
						(*y)(i, j) += mu * (*inv_w)(i, j) / fm;
						tv.push_back(std::pair<int64_t, double>(ijk, mu));


					}
				}
				mut.lock();
				for (int index = 0; index < tv.size(); index++) {
					next_correction->insert(tv[index]);
				}
				mut.unlock();
			}
		}
	}
}
