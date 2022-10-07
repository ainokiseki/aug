
Eigen::VectorXd norm_ori::get_direction(Eigen::SparseMatrix<double>* ata, Eigen::VectorXd& x, Eigen::VectorXd& w, Eigen::VectorXd& v, Eigen::VectorXd& rhs, int ld, double cst, int& cgnum) {

	//std::cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;
	Eigen::VectorXd aj2 = (((x + v / sigma).cwiseAbs() - w / sigma).array() <= 0.0).cast<double>().matrix();
	double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);
	(*ata) += sigma * aj2.asDiagonal();
	(*ata) += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
	//std::cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

	//Eigen::PardisoLDLT<SparseMatrix<double> > cg;
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
	double ita = std::fmin(0.5, 0.1 * rhs.norm());
	cg.setTolerance(ita / rhs.norm());
	cg.compute(*ata);

	Eigen::VectorXd direction = -cg.solve(rhs);
	cgnum = 666;
	//cgnum = cg.iterations();
	//std::cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
	//prt_res << cg.iterations() << "\t";
	//std::cout << "pro_res:" << (A * direction+rhs).norm() << endl;
	//std::cout << "CG:" << cg.iterations() << std::endl;
	//std::cout << "di:" << direction.norm() << endl;
	return direction;
}



double norm_ori::run_run(double cst_const, int con_const, int sub_iter) {

	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();


	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 1e-4;




	double c_t = 0;
	//std::unordered_map<int, int>colidx;


	timer time1, time2, time3, time4;

	cons_generator generator(n, mthread);



	while (true) {
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {

		}
		else {
			generator.generate(x, b, 1e-4);
		}

		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		if (viscons.size() < con_const * var_num) {


			for (int i = 0; i < n; i++) {
				int64_t index;
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						index = static_cast<int64_t>(i) * n * n + j * n + k;


						if (viscons.find(index) == viscons.end()) {
							for (int tp = 0; tp < 3; tp++) {
								if (colidx[idx[tp]] == -1) {

									colidx[idx[tp]] = colidx_num;
									colidx_num++;

								}
							}
							generator.data.push_back(cons_full(colidx[idx[0]], colidx[idx[1]], colidx[idx[2]]));
							generator.data.push_back(cons_full(colidx[idx[1]], colidx[idx[0]], colidx[idx[2]]));
							generator.data.push_back(cons_full(colidx[idx[2]], colidx[idx[1]], colidx[idx[0]]));
							viscons.insert(index);
							if (viscons.size() >= con_const * var_num)break;
						}

					}
					if (viscons.size() >= con_const * var_num)break;
				}
				if (viscons.size() >= con_const * var_num)break;
			}
		}
		int vssize = generator.data.size();

		double totalkkt = sqrt(generator.k1) / (1 + sqrt(generator.k2));
		cout << "totalkkt:" << totalkkt << endl;
		if (totalkkt < 1e-4)
			break;


		std::cout << "size:" << vssize << endl;

		sigma = 1e-2;
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



		cout << sub_var << endl;
		for (iter_all = 0; iter_all < 100; iter_all++) {


			x0 = sub_x;
			//cst = cst_const / sigma;
			cst = sigma * 1e-3;
			//if (pwd != "skip" && iter_all > 0) {
			//	ofstream model(string("model_save/") + pwd + string("_") + graph_name + string("_") + to_string(iter_all) + string(".dat"), std::ios::out | std::ios::binary);
			//	auto i_end = std::chrono::high_resolution_clock::now();
			//	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
			//	double tm = double(eld.count()) * 1.0 * 1e-9;
			//	model << setprecision(20);
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

			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {


				//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
				//std::cout << "x:" << x.norm() << endl;
				rhs = Eigen::VectorXd::Zero(sub_var);
				vector<int64_t>aj1;

				time1.tic();
				Eigen::SparseMatrix<double>ata(sub_var, sub_var);
				{
					Eigen::VectorXd diag = Eigen::VectorXd::Zero(sub_var);
					std::vector<Eigen::Triplet<double> >test_v;


					//m_calc_rhs_and_aj1(mthread, &rhs, &aj1, n, *axbu);
					//m_calc_rhs_and_aj1_diag(mthread, &rhs, &aj1, &diag, n, *axbu, &test_v);
					{

						//vector<int64_t>* maj1 = new vector<int64_t>[pnum];
						int var_num = rhs.size();
						vector<thread>tv;





						std::vector<almfunc_rhs_aj1_diag>vf;
						mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs_aj1_diag(&rhs, &diag, &test_v, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {

							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								n, std::ref(sub_x), std::ref(sub_b), sigma));
						}
						for (int i = 0; i < mthread; i++)
							tv[i].join();
						//cout << cnum[0]<<" "<<cnum[1]<<" "<<cnum[2] << endl;
					}

					ata.setFromTriplets(test_v.begin(), test_v.end());

					ata += diag.asDiagonal();
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




				ata *= sigma;
				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;

				rhs = sigma * rhs + sigma * prox_p_star(sub_x + sub_v / sigma, sigma, sub_w) + cst * (sub_x - x0);

				prt_res << iter_1 << "  " << rhs.norm() << "  ";


				if (rhs.norm() < rlim) {
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
				direction = get_direction(&ata, sub_x, sub_w, sub_v, rhs, 0, cst, cgnum);
				time2.toc();

				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;

				time3.tic();
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;
					double lx0 = 0;
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


					sub_x += direction;

					{
						int pnum = mthread;
						int var_num = nrhs.size();

						vector<thread>tv;

						std::vector<almfunc_rhs>vf;
						std::mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs(0, NULL, &nrhs, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < pnum; i++) {
							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								n, std::ref(sub_x), std::ref(sub_b), sigma));

						}
						for (int i = 0; i < pnum; i++)
							tv[i].join();
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);
					sub_x -= direction;


					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs;
					nrhs += sigma * prox_p_star(sub_x + direction + sub_v / sigma, sigma, sub_w) + cst * (sub_x + direction - x0);

					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";
						prt_res << "empty   ";
						return 1.0;

					}


					double alpha[20] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < 20; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}
					double lxnorm[20];
					for (int i = 0; i < 20; i++)
						lxnorm[i] = 0.0;

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					double* mphi = new double[pnum * 20];
					for (int i = 0; i < pnum * 20; i++)
						mphi[i] = 0;
					vector<thread>tv;
					std::vector<almfunc_superphi>vf;
					for (int i = 0; i < mthread; i++) {
						vf.push_back(almfunc_superphi(alpha, &direction, mphi + 20 * i, 20));
					}
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
							n, std::ref(sub_x), std::ref(sub_b), sigma));


					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++) {
						for (int j = 0; j < 20; j++) {
							lxnorm[j] += mphi[j + i * 20];
						}
					}

					delete[]mphi;





					lx0 = lxnorm[0];
					lx0 = sigma / 2 * lx0 + sigma / 2 * square(prox_p_star(sub_x + sub_v / sigma, sigma, sub_w).norm());

					lx0 += func_p(prox_p(sub_x + sub_v / sigma, sigma, sub_w), sub_w) + cst / 2 * (sub_x - x0).squaredNorm();



					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, 20));
					//(*axbu)();


					for (int i = 1; i < 20; i++) {

						double alp = alpha[i];

						double lx = lxnorm[i];

						lx = sigma / 2 * lx;
						lx += sigma / 2 * square(prox_p_star(sub_x + alp * direction + sub_v / sigma, sigma, sub_w).norm());
						lx += func_p(prox_p(sub_x + alp * direction + sub_v / sigma, sigma, sub_w), sub_w) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();

						if (lx <= down) {

							return alp;
						}

					}
					return alpha[19];





				}();
				time3.toc();
				//double alpha = 1;

				//std::cout << "alpha=" << alpha << endl;
				//if (rhs.norm() < 1e-5 && alpha < 1e-4)
				//	break;
				//flaga = true;

				sub_x += dalpha * direction;
				prt_res << dalpha << std::endl;
				std::cout << iter_1 << "\t" << rhs.norm() << "\t" << (rhs + ata * direction).norm() << "\t" << direction.norm() << "\t" << dalpha << endl;
			}
			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				int pnum = mthread;
				vector<thread>tv;

				std::vector<almfunc_u>vf;
				for (int i = 0; i < mthread; i++) {
					vf.push_back(almfunc_u(sigma));
				}
				int sz = generator.data.size();
				for (int i = 0; i < pnum; i++) {


					tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
						n, std::ref(sub_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				tv.clear();


				vector<double>vy_norm(mthread, 0), vkkt_y(mthread, 0), vkkt_u(mthread, 0);
				std::vector<almfunc_y>vf2;
				for (int i = 0; i < mthread; i++) {
					vf2.push_back(almfunc_y(sigma, &(vy_norm[i]), &(vkkt_y[i]), &(vkkt_u[i])));
				}


				for (int i = 0; i < pnum; i++) {
					tv.push_back(thread(new_axbu, std::ref(vf2[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
						n, std::ref(sub_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();
				kkt_y = std::accumulate(vkkt_y.begin(), vkkt_y.end(), 0.0);
				y_norm = std::accumulate(vy_norm.begin(), vy_norm.end(), 0.0);
				kkt_u = std::accumulate(vkkt_u.begin(), vkkt_u.end(), 0.0);

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

					int var_num = atu.size();
					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(var_num);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
							std::ref(mrhs[i])));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu + sub_v).norm() / (1 + sub_v.norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_x - z).norm() / (1 + sub_x.norm());
				kkt_z = (z - prox_p(z + sub_v, 1.0, sub_w)).norm() / (1 + z.norm());
				//for (int i = 0; i < var_num; i++)
				//    t_res += abs(x[i] * w[i]);
				//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
				double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
				std::cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
				std::cout << "kkt:" << kkt_final << endl;

				rp = std::fmax(kkt_u, kkt_v);
				rd = kkt_x;
				kkt = kkt_final;

				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;







			}
			time4.toc();



			prt_res << "==================\n";
			if (kkt < 1e-4)
				break;
			//rlim = k * 1e-4;
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

			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			kkt0 = kkt;

			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_x) << endl;
		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_x, sub_w) << endl;

		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]];
				b[i] = sub_b[colidx[i]];
				w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}
	}

	double dobj = 0;
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


	std::cout << "func_value:" << std::setprecision(8) << func_p(x, w) << "\t" << dobj << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	return func_p(x, w);




}



double norm_inf::run_run(double cst_const, int con_const, int sub_iter) {

	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();


	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 1e-3;

	double c_t = 0;
	//std::unordered_map<int, int>colidx;

	timer time1, time2, time3, time4;

	cons_generator generator(n, mthread);



	while (true) {
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(x, b, 1e-3);
		}
		double totalkkt = sqrt(generator.k1) / (1 + sqrt(generator.k2));
		cout << "totalkkt:" << totalkkt << endl;
		if (totalkkt < 1e-3)
			break;


		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = 1;
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



		cout << sub_var << endl;
		for (iter_all = 0; iter_all < 100; iter_all++) {


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

			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {


				//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
				//std::cout << "x:" << x.norm() << endl;
				rhs = Eigen::VectorXd::Zero(sub_var);
				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double>ata(sub_var, sub_var);
				{
					diag = Eigen::VectorXd::Zero(sub_var);
					std::vector<Eigen::Triplet<double> >test_v;

					//mthread = 1;
					//m_calc_rhs_and_aj1(mthread, &rhs, &aj1, n, *axbu);
					//m_calc_rhs_and_aj1_diag(mthread, &rhs, &aj1, &diag, n, *axbu, &test_v);
					{

						//vector<int64_t>* maj1 = new vector<int64_t>[pnum];
						int var_num = rhs.size();
						vector<thread>tv;





						std::vector<almfunc_rhs_aj1_diag>vf;
						mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs_aj1_diag(&rhs, &diag, &test_v, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {

							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
								n, std::ref(sub_x), std::ref(sub_b), sigma));
						}
						for (int i = 0; i < mthread; i++)
							tv[i].join();
						//cout << cnum[0]<<" "<<cnum[1]<<" "<<cnum[2] << endl;
					}

					ata.setFromTriplets(test_v.begin(), test_v.end());

					ata += diag.asDiagonal();
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




				ata *= sigma;
				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;

				rhs = sigma * rhs + sigma * sub_w.cwiseProduct(pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma)) + cst * (sub_x - x0);

				prt_res << iter_1 << "  " << rhs.norm() << "  ";


				if (rhs.norm() < rlim) {
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

					Eigen::VectorXd wx = sub_w.cwiseProduct(sub_x) + sub_v / sigma;

					Eigen::VectorXd r = pi_c(wx, sigma);
					r = Eigen::VectorXd::Ones(r.size()) - r.cwiseEqual(0).cast<double>();

					Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					double ita = std::fmin(0.5, 0.000001 * rhs.norm());
					cg.setTolerance(ita / rhs.norm());
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();

					if (wx.cwiseAbs().sum() < 1 / sigma) {
						ata += sigma * sub_w.cwiseProduct(sub_w).asDiagonal();
						cg.compute(ata);

						direction = -cg.solve(rhs);
					}
					else {
						if (true) {

							ata += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r).asDiagonal();
							double rnum = r.sum();
							Eigen::VectorXd rnew = r.cwiseProduct(wx.cwiseSign()).cwiseProduct(sub_w) * sqrt(sigma / rnum);

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

#else

							//{
							//	int k = 0; Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size()); Eigen::VectorXd r = rhs; double rou = r.squaredNorm(); double delta = std::fmin(0.5, 0.000001 * rhs.norm());
							//	Eigen::VectorXd p = r;
							//	Eigen::VectorXd omega;
							//	double tau, rou_last,mu;
							//	while (sqrt(rou) > delta) {
							//		k++;
							//		if (k > 1) {
							//			tau = rou / rou_last; p = r + tau * p;
							//		}
							//		omega = ata *p-rnew.dot(p)*rnew;
							//		mu = rou / (omega.dot(p)); x += mu * p; r -= mu * omega; rou_last = rou; rou = r.squaredNorm();
							//	}

							//	direction = -x;
							//	kcounter = k;
							//	//cout <<"eigencg:"<<cg2.iterations() << "cg:" << k<<" "<<(direction2-direction).norm() << endl;

							//}

							diag = sigma * diag + (a + cst) * Eigen::VectorXd::Ones(diag.size());
							diag += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r);
							diag = Eigen::VectorXd::Ones(diag.size()).cwiseQuotient(diag);

							//Eigen::VectorXd newdiag = ata.diagonal();
							//newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);




							kcounter = psqmr(ata, diag, rnew, rhs, direction);
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

							ata += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r.cwiseProduct(wx.cwiseSign().cwiseAbs())).asDiagonal();
							double rnum = r.sum();
							double rnum2 = r.dot(r);
							cout << rnum << " " << rnum2 << endl;
							Eigen::VectorXd rnew = r.cwiseProduct(wx.cwiseSign()).cwiseProduct(sub_w) * sqrt(sigma / rnum);
							cg.compute(ata);

							direction = -cg.solve(rhs);
							Eigen::VectorXd ar = cg.solve(rnew);
							direction += (ar * rnew.dot(direction)) / (1 - rnew.dot(ar));
							cout << "woodbury:" << direction.norm() << endl;

							//ata = ata - rnew * rnew.transpose();
							ata = ata - sigma * sub_w.cwiseProduct(wx.cwiseSign()).asDiagonal() * ((r * r.transpose()) / (r.dot(r))) * (wx.cwiseSign().cwiseProduct(sub_w).asDiagonal());
							cg.compute(ata);
							direction = -cg.solve(rhs);
							cout << "direct:" << direction.norm() << endl;



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

				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;

				time3.tic();
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;
					double lx0 = 0;
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


					sub_x += direction;

					{
						int pnum = mthread;
						int var_num = nrhs.size();

						vector<thread>tv;

						std::vector<almfunc_rhs>vf;
						std::mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs(0, NULL, &nrhs, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < pnum; i++) {
							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								n, std::ref(sub_x), std::ref(sub_b), sigma));

						}
						for (int i = 0; i < pnum; i++)
							tv[i].join();
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);
					sub_x -= direction;


					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs;
					nrhs += sigma * sub_w.cwiseProduct(pi_c(sub_w.cwiseProduct(sub_x + direction) + sub_v / sigma, sigma)) + cst * (sub_x + direction - x0);

					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";
						prt_res << "empty   ";
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}
					double lxnorm[searchnum];
					for (int i = 0; i < searchnum; i++)
						lxnorm[i] = 0.0;

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					double* mphi = new double[pnum * searchnum];
					for (int i = 0; i < pnum * searchnum; i++)
						mphi[i] = 0;
					vector<thread>tv;
					std::vector<almfunc_superphi>vf;
					for (int i = 0; i < mthread; i++) {
						vf.push_back(almfunc_superphi(alpha, &direction, mphi + searchnum * i, searchnum));
					}
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							n, std::ref(sub_x), std::ref(sub_b), sigma));


					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++) {
						for (int j = 0; j < searchnum; j++) {
							lxnorm[j] += mphi[j + i * searchnum];
						}
					}

					delete[]mphi;





					lx0 = lxnorm[0];
					lx0 = sigma / 2 * lx0 + sigma / 2 * pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma).squaredNorm();

					lx0 += func_p(sub_w.cwiseProduct(sub_x) + sub_v / sigma - pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma)) + cst / 2 * (sub_x - x0).squaredNorm();



					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();


					for (int i = 1; i < searchnum; i++) {

						double alp = alpha[i];

						double lx = lxnorm[i];

						lx = sigma / 2 * lx;
						lx += sigma / 2 * pi_c(sub_w.cwiseProduct(sub_x + alp * direction) + sub_v / sigma, sigma).squaredNorm();
						lx += func_p(sub_w.cwiseProduct(sub_x + alp * direction) + sub_v / sigma - pi_c(sub_w.cwiseProduct(sub_x + alp * direction) + sub_v / sigma, sigma)) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();

						if (lx <= down) {

							return alp;
						}

					}
					return alpha[19];





				}();
				time3.toc();
				//double alpha = 1;

				//std::cout << "alpha=" << alpha << endl;
				//if (rhs.norm() < 1e-5 && alpha < 1e-4)
				//	break;
				//flaga = true;

				sub_x += dalpha * direction;
				prt_res << dalpha << std::endl;
				std::cout << iter_1 << "\t" << rhs.norm() << "\t" << kcounter << "\t" << dalpha << endl;
			}
			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				int pnum = mthread;
				vector<thread>tv;

				std::vector<almfunc_u>vf;
				for (int i = 0; i < mthread; i++) {
					vf.push_back(almfunc_u(sigma));
				}
				int sz = generator.data.size();
				for (int i = 0; i < pnum; i++) {


					tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
						n, std::ref(sub_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				tv.clear();


				vector<double>vy_norm(mthread, 0), vkkt_y(mthread, 0), vkkt_u(mthread, 0);
				std::vector<almfunc_y>vf2;
				for (int i = 0; i < mthread; i++) {
					vf2.push_back(almfunc_y(sigma, &(vy_norm[i]), &(vkkt_y[i]), &(vkkt_u[i])));
				}


				for (int i = 0; i < pnum; i++) {
					tv.push_back(thread(new_axbu, std::ref(vf2[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
						n, std::ref(sub_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();
				kkt_y = std::accumulate(vkkt_y.begin(), vkkt_y.end(), 0.0);
				y_norm = std::accumulate(vy_norm.begin(), vy_norm.end(), 0.0);
				kkt_u = std::accumulate(vkkt_u.begin(), vkkt_u.end(), 0.0);

				//std::cout << "u:" << u.size() << endl;
				Eigen::VectorXd z = sub_w.cwiseProduct(sub_x) + sub_v / sigma - pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma);
				sub_v = sigma * pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma);

				y_norm = sqrt(y_norm);
				kkt_u = sqrt(kkt_u) / (1 + y_norm);
				kkt_y = sqrt(kkt_y) / (1 + y_norm);

				std::cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
				Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(sub_x.size());


				Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_x.size());



				//m_calc_atu(mthread, &atu, n, au);

				{

					int var_num = atu.size();
					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(var_num);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i])));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu + sub_w.cwiseProduct(sub_v)).norm() / (1 + sub_w.cwiseProduct(sub_v).norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_w.cwiseProduct(sub_x) - z).norm() / (1 + sub_w.cwiseProduct(sub_x).norm());
				kkt_z = (pi_c(z + sub_v, 1.0) - sub_v).norm() / (1 + z.norm());
				//for (int i = 0; i < var_num; i++)
				//    t_res += abs(x[i] * w[i]);
				//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
				double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
				std::cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
				std::cout << "kkt:" << kkt_final << endl;

				rp = std::fmax(kkt_u, kkt_v);
				rd = kkt_x;
				kkt = kkt_final;

				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;







			}
			time4.toc();



			prt_res << "==================\n";
			if (kkt < 1e-3) {
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

			//rlim = k * 1e-4;
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

			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_w.cwiseProduct(sub_x)) << endl;
			kkt0 = kkt;

		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_w.cwiseProduct(sub_x)) << endl;

		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]];
				b[i] = sub_b[colidx[i]];
				w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}
	}

	double dobj = 0;
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


	std::cout << "func_value:" << std::setprecision(8) << func_p(w.cwiseProduct(x)) << "\t" << dobj << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	return func_p(w.cwiseProduct(x));




}
double norm_inf::run_run_run(double cst_const, int con_const, int sub_iter) {

	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();
	b = b * scale_a;

	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = tol;

	double c_t = 0;
	//std::unordered_map<int, int>colidx;

	timer time1, time2, time3, time4, time5;

	cons_generator generator(n, mthread);



	while (true) {
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(x, b, scale_a * tol);
		}
		double totalkkt = sqrt(generator.k1) / (1 + sqrt(generator.k2));
		cout << "totalkkt:" << totalkkt << endl;
		if (totalkkt < tol)
			break;


		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = 1;
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



		cout << sub_var << endl;
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

			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {


				//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
				//std::cout << "x:" << x.norm() << endl;
				rhs = Eigen::VectorXd::Zero(sub_var);
				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double>ata(sub_var, sub_var);
				{
					diag = Eigen::VectorXd::Zero(sub_var);
					std::vector<Eigen::Triplet<double> >test_v;

					//mthread = 1;
					//m_calc_rhs_and_aj1(mthread, &rhs, &aj1, n, *axbu);
					//m_calc_rhs_and_aj1_diag(mthread, &rhs, &aj1, &diag, n, *axbu, &test_v);
					{

						//vector<int64_t>* maj1 = new vector<int64_t>[pnum];
						int var_num = rhs.size();
						vector<thread>tv;





						std::vector<almfunc_rhs_aj1_diag>vf;
						mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs_aj1_diag(&rhs, &diag, &test_v, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {

							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
								n, std::ref(sub_x), std::ref(sub_b), sigma));
						}
						for (int i = 0; i < mthread; i++)
							tv[i].join();
						//cout << cnum[0]<<" "<<cnum[1]<<" "<<cnum[2] << endl;
					}
					time5.tic();
					ata.setFromTriplets(test_v.begin(), test_v.end());

					ata += diag.asDiagonal();

					ata /= (scale_a * scale_a);
					time5.toc();
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




				ata *= sigma;
				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;

				rhs = sigma * rhs + sigma * sub_w.cwiseProduct(pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma)) + cst * (sub_x - x0);

				prt_res << iter_1 << "  " << rhs.norm() << "  ";


				if (rhs.norm() < rlim) {
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

					Eigen::VectorXd wx = sub_w.cwiseProduct(sub_x) + sub_v / sigma;

					Eigen::VectorXd r = pi_c(wx, sigma);
					r = Eigen::VectorXd::Ones(r.size()) - r.cwiseEqual(0).cast<double>();

					Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					double ita = std::fmin(0.5, 0.000001 * rhs.norm());
					cg.setTolerance(ita / rhs.norm());
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();

					if (wx.cwiseAbs().sum() < 1 / sigma) {
						ata += sigma * sub_w.cwiseProduct(sub_w).asDiagonal();
						cg.compute(ata);

						direction = -cg.solve(rhs);
					}
					else {
						if (true) {

							ata += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r).asDiagonal();
							double rnum = r.sum();
							Eigen::VectorXd rnew = r.cwiseProduct(wx.cwiseSign()).cwiseProduct(sub_w) * sqrt(sigma / rnum);

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

#else

							//{
							//	int k = 0; Eigen::VectorXd x = Eigen::VectorXd::Zero(rhs.size()); Eigen::VectorXd r = rhs; double rou = r.squaredNorm(); double delta = std::fmin(0.5, 0.000001 * rhs.norm());
							//	Eigen::VectorXd p = r;
							//	Eigen::VectorXd omega;
							//	double tau, rou_last,mu;
							//	while (sqrt(rou) > delta) {
							//		k++;
							//		if (k > 1) {
							//			tau = rou / rou_last; p = r + tau * p;
							//		}
							//		omega = ata *p-rnew.dot(p)*rnew;
							//		mu = rou / (omega.dot(p)); x += mu * p; r -= mu * omega; rou_last = rou; rou = r.squaredNorm();
							//	}

							//	direction = -x;
							//	kcounter = k;
							//	//cout <<"eigencg:"<<cg2.iterations() << "cg:" << k<<" "<<(direction2-direction).norm() << endl;

							//}

							//diag = sigma * diag + (a + cst) * Eigen::VectorXd::Ones(diag.size());
							//diag += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r);
							//diag = Eigen::VectorXd::Ones(diag.size()).cwiseQuotient(diag);

							Eigen::VectorXd newdiag = ata.diagonal();
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);




							kcounter = psqmr(ata, newdiag, rnew, rhs, direction, rlim);
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

							ata += sigma * sub_w.cwiseProduct(sub_w).cwiseProduct(r.cwiseProduct(wx.cwiseSign().cwiseAbs())).asDiagonal();
							double rnum = r.sum();
							double rnum2 = r.dot(r);
							cout << rnum << " " << rnum2 << endl;
							Eigen::VectorXd rnew = r.cwiseProduct(wx.cwiseSign()).cwiseProduct(sub_w) * sqrt(sigma / rnum);
							cg.compute(ata);

							direction = -cg.solve(rhs);
							Eigen::VectorXd ar = cg.solve(rnew);
							direction += (ar * rnew.dot(direction)) / (1 - rnew.dot(ar));
							cout << "woodbury:" << direction.norm() << endl;

							//ata = ata - rnew * rnew.transpose();
							ata = ata - sigma * sub_w.cwiseProduct(wx.cwiseSign()).asDiagonal() * ((r * r.transpose()) / (r.dot(r))) * (wx.cwiseSign().cwiseProduct(sub_w).asDiagonal());
							cg.compute(ata);
							direction = -cg.solve(rhs);
							cout << "direct:" << direction.norm() << endl;



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

				//Eigen::VectorXd direction = get_direction(rhs);
				//std::cout << "di:" << direction.norm() << endl;

				//for (int i = 0; i < var_num; i++)
				//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;

				time3.tic();
				double dalpha = [&]()->double
				{


					//cg end,result saved at direction
					//find step length

					const double c1 = 1e-4;
					double lx0 = 0;
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


					sub_x += direction;

					{
						int pnum = mthread;
						int var_num = nrhs.size();

						vector<thread>tv;

						std::vector<almfunc_rhs>vf;
						std::mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs(0, NULL, &nrhs, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < pnum; i++) {
							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								n, std::ref(sub_x), std::ref(sub_b), sigma));

						}
						for (int i = 0; i < pnum; i++)
							tv[i].join();
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);
					sub_x -= direction;


					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs;
					nrhs += sigma * sub_w.cwiseProduct(pi_c(sub_w.cwiseProduct(sub_x + direction) + sub_v / sigma, sigma)) + cst * (sub_x + direction - x0);

					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";
						prt_res << "empty   ";
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}
					double lxnorm[searchnum];
					for (int i = 0; i < searchnum; i++)
						lxnorm[i] = 0.0;

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					double* mphi = new double[pnum * searchnum];
					for (int i = 0; i < pnum * searchnum; i++)
						mphi[i] = 0;
					vector<thread>tv;
					std::vector<almfunc_superphi>vf;
					for (int i = 0; i < mthread; i++) {
						vf.push_back(almfunc_superphi(alpha, &direction, mphi + searchnum * i, searchnum));
					}
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							n, std::ref(sub_x), std::ref(sub_b), sigma));


					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++) {
						for (int j = 0; j < searchnum; j++) {
							lxnorm[j] += mphi[j + i * searchnum];
						}
					}

					delete[]mphi;





					lx0 = lxnorm[0];
					lx0 = sigma / 2 * lx0 + sigma / 2 * pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma).squaredNorm();

					lx0 += func_p(sub_w.cwiseProduct(sub_x) + sub_v / sigma - pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma)) + cst / 2 * (sub_x - x0).squaredNorm();



					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();


					for (int i = 1; i < searchnum; i++) {

						double alp = alpha[i];

						double lx = lxnorm[i];

						lx = sigma / 2 * lx;
						lx += sigma / 2 * pi_c(sub_w.cwiseProduct(sub_x + alp * direction) + sub_v / sigma, sigma).squaredNorm();
						lx += func_p(sub_w.cwiseProduct(sub_x + alp * direction) + sub_v / sigma - pi_c(sub_w.cwiseProduct(sub_x + alp * direction) + sub_v / sigma, sigma)) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();

						if (lx <= down) {

							return alp;
						}

					}
					return alpha[searchnum - 1];





				}();
				time3.toc();
				//double alpha = 1;

				//std::cout << "alpha=" << alpha << endl;
				//if (rhs.norm() < 1e-5 && alpha < 1e-4)
				//	break;
				//flaga = true;

				sub_x += dalpha * direction;
				prt_res << dalpha << std::endl;
				std::cout << iter_1 << "\t" << rhs.norm() << "\t" << kcounter << "\t" << dalpha << endl;
			}
			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				int pnum = mthread;
				vector<thread>tv;

				std::vector<almfunc_u>vf;
				for (int i = 0; i < mthread; i++) {
					vf.push_back(almfunc_u(sigma));
				}
				int sz = generator.data.size();
				for (int i = 0; i < pnum; i++) {


					tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
						n, std::ref(sub_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				tv.clear();


				vector<double>vy_norm(mthread, 0), vkkt_y(mthread, 0), vkkt_u(mthread, 0);
				std::vector<almfunc_y>vf2;
				for (int i = 0; i < mthread; i++) {
					vf2.push_back(almfunc_y(sigma, &(vy_norm[i]), &(vkkt_y[i]), &(vkkt_u[i])));
				}


				for (int i = 0; i < pnum; i++) {
					tv.push_back(thread(new_axbu, std::ref(vf2[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
						n, std::ref(sub_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();
				kkt_y = std::accumulate(vkkt_y.begin(), vkkt_y.end(), 0.0);
				y_norm = std::accumulate(vy_norm.begin(), vy_norm.end(), 0.0);
				kkt_u = std::accumulate(vkkt_u.begin(), vkkt_u.end(), 0.0);

				//std::cout << "u:" << u.size() << endl;
				Eigen::VectorXd z = sub_w.cwiseProduct(sub_x) + sub_v / sigma - pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma);
				sub_v = sigma * pi_c(sub_w.cwiseProduct(sub_x) + sub_v / sigma, sigma);

				y_norm = sqrt(y_norm);
				kkt_u = sqrt(kkt_u) / (1 + y_norm);
				kkt_y = sqrt(kkt_y) / (1 + y_norm);

				std::cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
				Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(sub_x.size());


				Eigen::VectorXd atu = Eigen::VectorXd::Zero(sub_x.size());



				//m_calc_atu(mthread, &atu, n, au);

				{

					int var_num = atu.size();
					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(var_num);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i])));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu + sub_w.cwiseProduct(sub_v)).norm() / (1 + sub_w.cwiseProduct(sub_v).norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_w.cwiseProduct(sub_x) - z).norm() / (1 + sub_w.cwiseProduct(sub_x).norm());
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

				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;




				rlim = fmin(0.1 * kkt, tol);


			}
			time4.toc();




			prt_res << "==================\n";
			if (kkt < tol) {
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

			//rlim = k * 1e-4;
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

			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_w.cwiseProduct(sub_x)) << endl;
			kkt0 = kkt;

		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_w.cwiseProduct(sub_x)) << endl;

		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]];
				b[i] = sub_b[colidx[i]];
				w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}
	}

	double dobj = 0;
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


	std::cout << "func_value:" << std::setprecision(8) << func_p(w.cwiseProduct(x)) << "\t" << dobj << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	std::cout << "total time is:" << tm << endl;
	std::cout << "cons generate time is:" << c_t << endl;
	cout << time1.time() << " " << time2.time() << " " << time3.time() << " " << time4.time() << endl;
	cout << "diana:" << time5.time() << endl;
	x /= scale_a;
	b /= scale_a;

	return func_p(w.cwiseProduct(x));




}

double norm_inf::run_wx(double cst_const, int con_const, int sub_iter) {

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
	//mthread = 1;
	cout << "m:" << mthread << endl;
	cons_generator generator(n, mthread);
	double total_kkt_u = 0, total_y_norm = 0;


	while (true) {
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {
			//generator.alm_generate(n, viscons, colidx, colidx_num, x, b, 1e-4);
		}
		else {
			generator.generate(x, b, scale_a * tol);
		}
		double totalkkt = sqrt(generator.k1 + total_kkt_u) / (1 + sqrt(generator.k2 + total_y_norm));
		cout << "totalkkt:" << totalkkt << endl;
		if (totalkkt < tol)
			break;


		std::vector<int>& colidx(generator.colidx);
		std::set<int64_t>& viscons(generator.viscons);
		int colidx_num = generator.colidx_num;


		int vssize = generator.data.size();



		std::cout << "size:" << vssize << endl;

		sigma = 1;
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
		Eigen::VectorXd ori_x;

		cout << sub_var << endl;
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

			for (iter_1 = 0; iter_1 < sub_iter; iter_1++) {
				ori_x = inv_w.cwiseProduct(sub_x);

				//std::coutstd::cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
				//std::cout << "x:" << x.norm() << endl;
				rhs = Eigen::VectorXd::Zero(sub_var);
				vector<int64_t>aj1;
				Eigen::VectorXd diag;
				time1.tic();
				Eigen::SparseMatrix<double>ata(sub_var, sub_var);
				{
					diag = Eigen::VectorXd::Zero(sub_var);
					std::vector<Eigen::Triplet<double> >test_v;

					//mthread = 1;
					//m_calc_rhs_and_aj1(mthread, &rhs, &aj1, n, *axbu);
					//m_calc_rhs_and_aj1_diag(mthread, &rhs, &aj1, &diag, n, *axbu, &test_v);
					{

						//vector<int64_t>* maj1 = new vector<int64_t>[pnum];
						int var_num = rhs.size();
						vector<thread>tv;





						std::vector<almfunc_rhs_aj1_diag>vf;
						mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs_aj1_diag(&rhs, &diag, &test_v, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < mthread; i++) {

							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
								n, std::ref(ori_x), std::ref(sub_b), sigma));
						}
						for (int i = 0; i < mthread; i++)
							tv[i].join();
						//cout << cnum[0]<<" "<<cnum[1]<<" "<<cnum[2] << endl;
					}
					time5.tic();
					ata.setFromTriplets(test_v.begin(), test_v.end());

					ata += diag.asDiagonal();

					ata /= (scale_a * scale_a);
					time5.toc();
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



				ata = inv_w.asDiagonal() * ata * inv_w.asDiagonal();
				ata *= sigma;
				//std::cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
				//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


				//std::cout << "r1:" << rhs.norm() << endl;
				//cout << "rhsnorm1:" << rhs.norm() << endl;
				rhs = sigma * rhs.cwiseProduct(inv_w) + sigma * pi_c(sub_x + sub_v / sigma, sigma) + cst * (sub_x - x0);
				//cout << "rhsnorm2:" << rhs.norm() << endl;
				prt_res << iter_1 << "  " << rhs.norm() << "  ";


				if (rhs.norm() < rlim) {
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

					Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
					double ita = tol;
					cg.setTolerance(ita / rhs.norm());
					ata += (a + cst) * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();

					if (wx.cwiseAbs().sum() < 1 / sigma) {
						//ata += sigma * sub_w.cwiseProduct(sub_w).asDiagonal();

						//ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();
						//cg.compute(ata);

						//direction = -cg.solve(rhs);
						//kcounter = cg.iterations();
						ata += sigma * Eigen::VectorXd::Ones(rhs.size()).asDiagonal();


						Eigen::VectorXd newdiag = ata.diagonal();
						newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);

						Eigen::VectorXd rnew = Eigen::VectorXd::Zero(rhs.size());
						kcounter = psqmr(ata, newdiag, rnew, rhs, direction, rlim);




					}
					else {
						if (true) {

							ata += sigma * r.asDiagonal();
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

							Eigen::VectorXd newdiag = ata.diagonal();
							newdiag = Eigen::VectorXd::Ones(newdiag.size()).cwiseQuotient(newdiag);




							kcounter = psqmr(ata, newdiag, rnew, rhs, direction, rlim);
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
					double lx0 = 0;
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


					ori_x += wdirection;

					{
						int pnum = mthread;
						int var_num = nrhs.size();

						vector<thread>tv;

						std::vector<almfunc_rhs>vf;
						std::mutex mu;
						for (int i = 0; i < mthread; i++) {
							vf.push_back(almfunc_rhs(0, NULL, &nrhs, &mu));
						}
						int sz = generator.data.size();
						for (int i = 0; i < pnum; i++) {
							tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + static_cast<int64_t>(i) * sz / mthread, generator.data.begin() + static_cast<int64_t>(i + 1) * sz / mthread,
								n, std::ref(ori_x), std::ref(sub_b), sigma));

						}
						for (int i = 0; i < pnum; i++)
							tv[i].join();
					}
					//m_calc_rhs(mthread, &nrhs, n, *axbu);
					ori_x -= wdirection;


					//std::cout << nrhs.norm() << endl;
					nrhs = sigma * nrhs.cwiseProduct(inv_w);
					nrhs += sigma * pi_c(sub_x + direction + sub_v / sigma, sigma) + cst * (sub_x + direction - x0);

					if (nrhs.dot(direction) < 0) {
						//std::cout << "green\n";
						prt_res << "empty   ";
						return 1.0;

					}


					double alpha[searchnum] = { 0.0 };
					alpha[0] = 0; alpha[1] = 1;
					for (int i = 2; i < searchnum; i++) {
						alpha[i] = alpha[i - 1] / 2;
					}
					double lxnorm[searchnum];
					for (int i = 0; i < searchnum; i++)
						lxnorm[i] = 0.0;

					//m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
					int sz = generator.data.size();
					int pnum = mthread;
					double* mphi = new double[pnum * searchnum];
					for (int i = 0; i < pnum * searchnum; i++)
						mphi[i] = 0;
					vector<thread>tv;
					std::vector<almfunc_superphi>vf;
					for (int i = 0; i < mthread; i++) {
						vf.push_back(almfunc_superphi(alpha, &wdirection, mphi + searchnum * i, searchnum));
					}
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							n, std::ref(ori_x), std::ref(sub_b), sigma));


					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++) {
						for (int j = 0; j < searchnum; j++) {
							lxnorm[j] += mphi[j + i * searchnum];
						}
					}

					delete[]mphi;





					lx0 = lxnorm[0];
					lx0 = sigma / 2 * lx0 + sigma / 2 * pi_c(sub_x + sub_v / sigma, sigma).squaredNorm();

					lx0 += func_p(sub_x + sub_v / sigma - pi_c(sub_x + sub_v / sigma, sigma)) + cst / 2 * (sub_x - x0).squaredNorm();



					double down = lx0 + c1 * rhs.dot(direction);
					//std::cout << "down:" << down << "\tlx0:" << lx0 << endl;

					//std::cout << "lx0:" << lx0 << "\tdown" << down << endl;
					prt_res << lx0 << "  ";


					//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, searchnum));
					//(*axbu)();


					for (int i = 1; i < searchnum; i++) {

						double alp = alpha[i];

						double lx = lxnorm[i];

						lx = sigma / 2 * lx;
						lx += sigma / 2 * pi_c(sub_x + alp * direction + sub_v / sigma, sigma).squaredNorm();
						lx += func_p(sub_x + alp * direction + sub_v / sigma - pi_c(sub_x + alp * direction + sub_v / sigma, sigma)) + cst / 2 * (sub_x + alp * direction - x0).squaredNorm();

						if (lx <= down) {

							return alp;
						}

					}
					return alpha[searchnum - 1];





				}();
				time3.toc();
				//double alpha = 1;



				sub_x += dalpha * direction;
				//cout << "x:" << setprecision(10) << x.norm() << endl;
				//cout << "x:" << sub_x.norm() << endl;
				prt_res << dalpha << std::endl;
				std::cout << setprecision(10) << iter_1 << "\t" << rhs.norm() << "\t" << kcounter << "\t" << dalpha << endl;
			}


			double rp, rd;
			time4.tic();

			{


				double y_norm = 0;
				double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
				kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
				int pnum = mthread;
				vector<thread>tv;

				//std::vector<almfunc_u>vf;
				//for (int i = 0; i < mthread; i++) {
				//	vf.push_back(almfunc_u(sigma));
				//}
				int sz = generator.data.size();
				//for (int i = 0; i < pnum; i++) {


				//	tv.push_back(thread(new_axbu, std::ref(vf[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
				//		n, std::ref(ori_x), std::ref(sub_b), sigma));

				//}
				//for (int i = 0; i < pnum; i++)
				//	tv[i].join();
				//tv.clear();


				vector<double>vy_norm(mthread, 0), vkkt_y(mthread, 0), vkkt_u(mthread, 0);
				std::vector<almfunc_uy>vf2;
				for (int i = 0; i < mthread; i++) {
					vf2.push_back(almfunc_uy(sigma, &(vy_norm[i]), &(vkkt_y[i]), &(vkkt_u[i])));
				}


				for (int i = 0; i < pnum; i++) {
					tv.push_back(thread(new_axbu, std::ref(vf2[i]), generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
						n, std::ref(ori_x), std::ref(sub_b), sigma));

				}
				for (int i = 0; i < pnum; i++)
					tv[i].join();
				//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
				//(*axbu)();
				kkt_y = std::accumulate(vkkt_y.begin(), vkkt_y.end(), 0.0);
				y_norm = std::accumulate(vy_norm.begin(), vy_norm.end(), 0.0);
				kkt_u = std::accumulate(vkkt_u.begin(), vkkt_u.end(), 0.0);
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

					int var_num = atu.size();
					Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

					for (int i = 0; i < pnum; i++)
						mrhs[i] = Eigen::VectorXd::Zero(var_num);
					vector<thread>tv;
					for (int i = 0; i < pnum; i++) {

						tv.push_back(thread(new_atu, generator.data.begin() + i * sz / mthread, generator.data.begin() + (i + 1) * sz / mthread,
							std::ref(mrhs[i])));
					}
					for (int i = 0; i < pnum; i++)
						tv[i].join();
					for (int i = 0; i < pnum; i++)
						atu += mrhs[i];
					delete[]mrhs;

				}
				//cout <<"sub_v:"<<sub_v.norm()<< "\tatu:" << atu.norm() << endl;
				//kkt_x = au_eigen(n, &u, v, tmp1,cons)();
				kkt_x = (atu.cwiseProduct(inv_w) + sub_v).norm() / (1 + sub_v.norm());
				//kkt_x = kkt_x / (1 + v.norm());

				//std::cout << "z_norm" << z.norm() << endl;
				kkt_v = (sub_x - z).norm() / (1 + sub_x.norm());
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

				//std::cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;







			}
			time4.toc();




			prt_res << "==================\n";
			if (kkt < tol) {
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

			std::cout << "primal:" << primal_win << "\t" << "dual:" << dual_win << " " << kkt0 << endl;
			std::cout << "prim_func_value:" << std::setprecision(8) << func_p(sub_x / scale_a) << endl;
			kkt0 = kkt;
			rlim = std::fmax(tol, 0.1 * kkt);
		}
		std::cout << "sub_func_value:" << std::setprecision(8) << func_p(sub_x / scale_a) << endl;

		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {


				x[i] = sub_x[colidx[i]] / sub_w[colidx[i]];
				//b[i] = sub_b[colidx[i]];
				//w[i] = sub_w[colidx[i]];
				v[i] = sub_v[colidx[i]];
			}
		}
	}

	double dobj = 0;
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


	std::cout << "func_value:" << std::setprecision(8) << func_p(w.cwiseProduct(x)) << "\t" << dobj << endl;
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
	return func_p(w.cwiseProduct(x));




}

