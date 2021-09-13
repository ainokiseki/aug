#include"admm_class.h"
#include"mytime.h"


admm_square::admm_square(string graph_name, int type) :admm(graph_name, type)
{



}

admm_square::~admm_square()
{

}
double admm_square::func_p(Eigen::VectorXd x, Eigen::VectorXd w) {
	return w.cwiseProduct(x).squaredNorm();
}


double admm_square::run(double tau, double sigma, double scale_b, double scale_c) {
	vector<vector<pair<int64_t, double>>>u1(mthread), u2(mthread);
	vector<vector<pair<int64_t, double>>>& u = u1;
	vector<vector<pair<int64_t, double>>>& u_new = u2;
	double scale_a;
	//scale_a = 3;
	//scale_b = b.norm();
	//b = b / b.norm();
	//w = w / 9;
	//scale_b = b.norm();
	//scale_c = 3;
	w = w.cwiseSqrt();
	std::vector<int> colidx(var_num, -1);
	int colidx_num = 0;
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

	timer madoka, homula;
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
				vt.push_back(thread(&admm::constrain_generate, this, ta, tb, ref(mut), ref(colidx), ref(colidx_num), ref(acons), ref(cons), ref(viscons)));
			}
			for (int i = 0; i < mthread; i++) {
				vt[i].join();
			}
		}
		cout << "cons:" << viscons.size() << endl;

		if (loopflag == 2 && old_size == viscons.size()) {


			cout << "timeis:" << madoka.time() << " homulais:" << homula.time() << endl;

			x = x / scale_a * scale_b;
			b = b / 3;
			return finalres;
		}
		else if (loopflag == 2 || old_size == viscons.size()) {
			loopflag = 1;
		}

		old_size = viscons.size();
		Eigen::SparseMatrix<double> hst(viscons.size() * 3, colidx_num);
		Eigen::SparseMatrix<double> hs(colidx_num, colidx_num);
		hst.setFromTriplets(acons.begin(), acons.end());
		//acons.clear();

		hs = hst.transpose() * hst / (scale_a * scale_a);

		//hs += Eigen::VectorXd::Ones(colidx_num).asDiagonal();
		Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
		//Eigen::PardisoLDLT<SparseMatrix<double>> cg;
		//cg.compute(hs);

		sub_x = Eigen::VectorXd::Zero(colidx_num);
		sub_b = Eigen::VectorXd::Zero(colidx_num);
		sub_w = Eigen::VectorXd::Zero(colidx_num);




		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				sub_x[colidx[i]] = x[i];
				sub_b[colidx[i]] = b[i];
				sub_w[colidx[i]] = w[i];
			}
		}


		int sub_var = colidx_num;
		bool sigmaflag = true;
		double sigma_old = 1;

		hs += 2 / sigma_old * sub_w.cwiseProduct(sub_w).asDiagonal();
		for (int iter_I = 0;; iter_I++)
		{
			if (sigmaflag) {
				sigmaflag = false;

				hs-= 2 / sigma_old * sub_w.cwiseProduct(sub_w).asDiagonal();
				hs += 2 / sigma * sub_w.cwiseProduct(sub_w).asDiagonal();
				cg.compute(hs);
				sigma_old = sigma;
			}


			madoka.tic();
			Eigen::VectorXd old_x = sub_x;
			rhs = Eigen::VectorXd::Zero(sub_var);
			//atyu(n, x, b, &u, rhs, sigma,scale_a);
			//cout << u.size() << endl;
			m_atyu(n, sub_x, sub_b, u, rhs, sigma, cons, colidx, scale_a);
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



			m_all_in_all(n, u, u_new, sub_x, sub_b, atu, old_x, kkt_u, y_norm, sigma, cons, colidx, tau, scale_a);
			std::swap(u, u_new);

			Eigen::VectorXd wtwx = 2 * sub_w.cwiseProduct(sub_w).cwiseProduct(sub_x);

			double kkt_x = (atu + wtwx).norm() / (1 + wtwx.norm());
			double true_ku = scale_b * sqrt(kkt_u) / (1 + scale_b * sqrt(y_norm));
			kkt_u = sqrt(kkt_u) / (1 + sqrt(y_norm));
			

			//cout << "allinall\n";
			if (kkt_u > kkt_x) {
				dual_win++;
			}
			else {
				prim_win++;
			}
			if ((iter_I + 1) % itertime(iter_I) == 0) {
				if (prim_win > dual_win * 1.2 && prim_win > 1) {
					prim_win = 0;
					sigma /= 3;
					sigmaflag = true;
				}
				else if (dual_win > prim_win * 1.2 && dual_win > 1) {
					dual_win = 0;
					sigma *= 3;
					sigmaflag = true;
				}
			}
			kkt_x = scale_c * (atu + wtwx).norm() / (1 + scale_c * wtwx.norm());
			//double kkt_z=scale_b/scale_a*(z - prox_p(z + sub_v, 1.0, sub_w)).norm() / (1 + scale_b / scale_a * z.norm());
			

			//cout << kkt_x << " " << kkt_v << " " << true_ku <<" !:"<<kkt_z<< endl;
			double kkt = std::fmax(kkt_x, true_ku);
			maxsize = fmax(maxsize, u.size());

			cout << iter_I << "\t" << kkt << "\t" << func_p(sub_x, sub_w) * scale_b * scale_c / 3 << "\t" << sigma << "\t" << u.size() << endl;
			if (kkt < 1e-3) {


				double pobj = func_p(sub_x, sub_w) * scale_b / 3;
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
				double d2 = atu.cwiseQuotient(sub_w).squaredNorm() / 12;
				dobj -= d2;
				cout << "pobj:" << pobj << " dobj:" << dobj <<" d2:"<<d2<< endl;
				if (fabs(pobj - dobj) / dobj < 1e-2) {

					if (loopflag == 0) {

						for (int i = 0; i < var_num; i++) {
							if (colidx[i] != -1) {
								x[i] = sub_x[colidx[i]];
								b[i] = sub_b[colidx[i]];
								w[i] = sub_w[colidx[i]];
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