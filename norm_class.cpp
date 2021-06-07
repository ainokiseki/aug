
#include"norm_class.h"
#include<iomanip>

void m_calc_uy(int pnum, axbu_eigen axbu, std::map<int64_t,double>*u,double sigma,double*yn,double*ky,double*ku) {

	int n = axbu.n;
	double* y_norm = new double[pnum];
	double* kkt_y = new double[pnum];
	double* kkt_u = new double[pnum];
	for (int i = 0; i < pnum; i++) {
		y_norm[i] = kkt_y[i] = kkt_u[i] = 0;
	}
	std::mutex mu;
	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.set_function(new func_eigen_uy_mul(u,sigma,y_norm+i,kkt_y+i,kkt_u+i,&mu));
		int a, b;
		a = n - (int)(cbrt((i + 1) * 1.0 / pnum) * n);
		b = n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n);
		tv.push_back(thread(axbu, a, b));

	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++) {
		*yn += y_norm[i];
		*ky += kkt_y[i];
		*ku += kkt_u[i];
	}

	delete[]y_norm;
	delete[]kkt_y;
	delete[]kkt_u;
}




void m_calc_superphi(int pnum, axbu_eigen axbu, double* alpha, Eigen::VectorXd* direction,double*lxnorm) {
	int n = axbu.n;
	double* mphi = new double[pnum*20];
	for (int i = 0; i < pnum*20; i++)
		mphi[i] = 0;
	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.set_function(new func_eigen_superphi(alpha, direction, mphi+20*i,20));
		int a, b;
		a = n - (int)(cbrt((i + 1) * 1.0 / pnum) * n);
		b = n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n);
		tv.push_back(thread(axbu, a, b));

	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++) {
		for (int j = 0; j < 20; j++) {
			lxnorm[j] += mphi[j + i * 20];
		}
	}

	delete[]mphi;
}

void m_calc_rhs_and_aj1_diag(int pnum, Eigen::VectorXd* rhs, vector<int64_t>* aj1,Eigen::VectorXd*diag, int n, axbu_eigen axbu, vector<Eigen::Triplet<double> >*trip) {
	int var_num = n * (n - 1) / 2;
	Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];
	Eigen::VectorXd* mdiag = new Eigen::VectorXd[pnum];
	vector< vector<int64_t>*>maj1(pnum);
	vector< vector<Eigen::Triplet<double> >*>mtrip(pnum);
	for (int i = 0; i < pnum; i++) {
		maj1[i] = new vector<int64_t>;
		mtrip[i] = new vector<Eigen::Triplet<double> >;
	}
	//vector<int64_t>* maj1 = new vector<int64_t>[pnum];
	for (int i = 0; i < pnum; i++) {
		mdiag[i] = Eigen::VectorXd::Zero(var_num);

		mrhs[i] = Eigen::VectorXd::Zero(var_num);
	}

	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.set_function(new func_eigen_rhs_aj1_diag(mrhs + i, maj1[i],mdiag+i,mtrip[i]));
		tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++) {
		*rhs += mrhs[i];
		*diag += mdiag[i];
	}

	for (int i = 0; i < pnum; i++) {
		delete maj1[i];
		trip->insert(trip->end(), mtrip[i]->begin(), mtrip[i]->end());
		delete mtrip[i];
	}

	delete[]mrhs;

	delete[]mdiag;
}
void m_calc_rhs_and_aj1(int pnum, Eigen::VectorXd* rhs, vector<int64_t>* aj1, int n, axbu_eigen axbu) {
	int var_num = n * (n - 1) / 2;
	Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];
	vector<int64_t>* maj1 = new vector<int64_t>[pnum];
	for (int i = 0; i < pnum; i++)
		mrhs[i] = Eigen::VectorXd::Zero(var_num);
	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.set_function(new func_eigen_rhs_aj1(mrhs + i, maj1 + i));
		tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++)
		*rhs += mrhs[i];
	for (int i = 0; i < pnum; i++) {
		aj1->insert(aj1->end(), maj1[i].begin(), maj1[i].end());
		maj1[i].~vector();
	}

	delete[]mrhs;
	delete[]maj1;
}
void m_calc_rhs(int pnum, Eigen::VectorXd* rhs, int n, axbu_eigen axbu) {
	int var_num = n * (n - 1) / 2;
	Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

	for (int i = 0; i < pnum; i++)
		mrhs[i] = Eigen::VectorXd::Zero(var_num);
	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.set_function(new func_eigen_rhs_aj1(mrhs + i, NULL));
		tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++)
		* rhs += mrhs[i];
	delete[]mrhs;
}
void m_calc_phi(int pnum, double& phi, int n, axbu_eigen axbu, double alpha, Eigen::VectorXd* direction) {
	int var_num = n * (n - 1) / 2;
	double* mphi = new double[pnum];
	for (int i = 0; i < pnum; i++)
		mphi[i] = 0;
	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.set_function(new func_eigen_phi(alpha, direction, mphi + i, i));
		int a, b;
		a = n - (int)(cbrt((i + 1) * 1.0 / pnum) * n);
		b = n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n);
		tv.push_back(thread(axbu, a, b));

	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++)
		phi += mphi[i];

	delete[]mphi;
}

norm_ori::norm_ori(string graph_name) :metric(graph_name)
{
	iter_all = 0;
}
norm_ori::norm_ori(string graph_name,string data) : metric(graph_name)
{
	ifstream f(data,std::ios::in|std::ios::binary);
	string s;
	f >> s;
	f >> iter_all;
	f >> n >> var_num >> sigma;
	double tm;
	f >> tm;
	cout << "last time:"<<tm << endl;
	time_record = tm;
	for (int i = 0; i < var_num; i++) {
		f >> tm;
		x[i] = tm;
	}
	for (int i = 0; i < var_num; i++) {
		f >> tm;
		v[i] = tm;
	}
	int unum;
	f >> unum;
	int64_t tm2;
	for (int i = 0; i < unum; i++) {
		f >> tm2 >> tm;
		u.insert(pair<int64_t, double>(tm2, tm));
	}
	cout << "reload succeed\n" << "u:" << u.size() << endl;
}
norm_ori::~norm_ori()
{
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
double norm_ori::get_kkt(double sigma) {

	double y_norm = 0;
	double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
	kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
	//axbu->set_function(new func_eigen_u(&u, sigma));
	//(*axbu)();

	//axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
	//(*axbu)();
	//axbu->set_function(new func_eigen_uy(&u, sigma, &y_norm, &kkt_y, &kkt_u));
	//(*axbu)();
	m_calc_uy(mthread, *axbu, &u, sigma, &y_norm, &kkt_y, &kkt_u);

	cout << "u:" << u.size() << endl;
	Eigen::VectorXd z = prox_p(x+v/sigma,sigma);
	v = sigma * prox_p_star(x + v / sigma, sigma);

	y_norm = sqrt(y_norm);
	kkt_u = sqrt(kkt_u) / (1 + y_norm);
	kkt_y = sqrt(kkt_y) / (1 + y_norm);
	cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
	Eigen::VectorXd tmp1 = Eigen::VectorXd::Ones(x.size());

	kkt_x = au_eigen(n, &u, v, tmp1)();

	kkt_x = kkt_x / (1 + v.norm());

	//cout << "z_norm" << z.norm() << endl;
	kkt_v = (x - z).norm() / (1 + z.norm());
	kkt_z = (z - prox_p(z+v,1.0)).norm() / (1 + z.norm());
	//for (int i = 0; i < var_num; i++)
	//    t_res += abs(x[i] * w[i]);
	//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
	double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
	cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
	cout << "kkt:" << kkt_final << endl;


	return kkt_final;

	//cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;




}

double norm_ori::get_step_length(Eigen::VectorXd& direction, double sigma) {

	//cg end,result saved at direction
	//find step length
	double alpha = 1;
	const double c1 = 1e-4;
	double lx0 = 0;
	double galp = 0;
	double g0 = direction.dot(rhs);
	if (g0 >= 0) {
		std::cout << "Need a descent direction\n";
		exit(-1);
	}
	if (mthread == 0) {
		axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
		(*axbu)();
	}
	else {
		m_calc_phi(mthread, lx0, n, *axbu, 0, &direction);
	}


	auto square = [](double x) {return x * x; };


	lx0 = sigma / 2 * lx0 + sigma / 2 * square(prox_p_star(x + v / sigma, sigma).norm());

	lx0 += func_p(prox_p(x+v/sigma,sigma));
	


	double down = lx0 + c1 * alpha * rhs.dot(direction);
	//cout << "down:" << down << "\tlx0:" << lx0 << endl;

	cout << "lx0:" << lx0 << "\tdown" << down << endl;
	prt_res << lx0 << "  ";

	Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(var_num);
	if (mthread == 0) {
		axbu->set_function(new func_eigen_rhs(alpha, &direction, &nrhs));
		(*axbu)();
	}
	else {
		x += direction;
		m_calc_rhs(mthread, &nrhs, n, *axbu);
		x -= direction;
	}


	nrhs = sigma * nrhs;
	nrhs += sigma * prox_p_star(x+alpha*direction + v / sigma, sigma);

	if (nrhs.dot(direction) < 0) {
		cout << "green\n";
		return 1.0;
	}

	double lb, ub;
	lb = 0; ub = 1;
	double glb = -1; double gub = 1;
	for (int i = 0; i < 20; i++) {
		nrhs = Eigen::VectorXd::Zero(var_num);

		//cout << "start find step,alpha=" << alpha << endl;
		double lx = 0;
		if (mthread==0) {
			axbu->set_function(new func_eigen_phi(alpha, &direction, &lx));
			(*axbu)();
		}
		else {
			m_calc_phi(mthread, lx, n, *axbu, alpha, &direction);
		}

		lx = sigma / 2 * lx;
		lx += sigma / 2 * square(prox_p_star(x+alpha*direction + v / sigma, sigma).norm());
		lx += func_p(prox_p(x + alpha*direction +v / sigma, sigma));

		if (lx <= down) {
			break;
		}

		//if (mthread == 0) {
		//	axbu->set_function(new func_eigen_rhs(alpha, &direction, &nrhs));
		//	(*axbu)();
		//}
		//else {
		//	x += alpha*direction;
		//	m_calc_rhs(mthread, &nrhs, n, *axbu);
		//	x -= alpha*direction;
		//}
		//nrhs *= sigma;
		//nrhs += sigma * prox_p_star(x + alpha * direction + v / sigma, sigma);
		//
		//galp = nrhs.dot(direction);
		//if (galp*gub < 0) {
		//	lb = alpha; glb = galp;
		//}
		//else {
		//	ub = alpha; gub = galp;
		//}
		ub = alpha;
		alpha = 0.5 * (lb + ub);
	}
	return alpha;


}
double norm_ori::get_step_length(Eigen::VectorXd& direction, double sigma,int numer) {

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



	Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(var_num);
	if (mthread == 0) {
		axbu->set_function(new func_eigen_rhs(1.0, &direction, &nrhs));
		(*axbu)();
	}
	else {
		x += direction;
		m_calc_rhs(mthread, &nrhs, n, *axbu);
		x -= direction;
	}


	nrhs = sigma * nrhs;
	nrhs += sigma * prox_p_star(x + direction + v / sigma, sigma);

	if (nrhs.dot(direction) < 0) {
		cout << "green\n";
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
	m_calc_superphi(mthread, *axbu, alpha, &direction, lxnorm);
	lx0 = lxnorm[0];
	lx0 = sigma / 2 * lx0 + sigma / 2 * square(prox_p_star(x + v / sigma, sigma).norm());

	lx0 += func_p(prox_p(x + v / sigma, sigma));



	double down = lx0 + c1 * rhs.dot(direction);
	//cout << "down:" << down << "\tlx0:" << lx0 << endl;

	cout << "lx0:" << lx0 << "\tdown" << down << endl;
	prt_res << lx0 << "  ";


	//axbu->set_function(new func_eigen_superphi(alpha, &direction, lxnorm, 20));
	//(*axbu)();


	for (int i = 1; i < 20; i++) {
		double alp = alpha[i];

		double lx = lxnorm[i];

		lx = sigma / 2 * lx;
		lx += sigma / 2 * square(prox_p_star(x + alp * direction + v / sigma, sigma).norm());
		lx += func_p(prox_p(x + alp * direction + v / sigma, sigma));

		if (lx <= down) {
			return alp;
		}

	}
	return alpha[19];


}
Eigen::VectorXd norm_ori::get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1) {
	Eigen::VectorXd aj2 = (((x + v / sigma).cwiseAbs() - w / sigma).array() <= 0.0).cast<double>().matrix() ;

	//cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;

	double a = min(pow(rhs.norm(), 1.5) * 0.5,20.0);

	MatrixReplacement A(aj1, n, a, &aj2, sigma);

	//Eigen::BiCGSTAB<MatrixReplacement,  Eigen::IdentityPreconditioner> cg;
	Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
	double ita = std::fmin(0.5, rhs.norm());
	cg.setTolerance(ita);
	cg.compute(A);

	Eigen::VectorXd direction = cg.solve(-rhs);

	cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
	prt_res << cg.iterations() << "\t";
	//cout << "pro_res:" << (A * direction+rhs).norm() << endl;
	//cout << "CG:" << cg.iterations() << std::endl;
	//cout << "di:" << direction.norm() << endl;
	return direction;
}
Eigen::VectorXd norm_ori::get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1,Eigen::VectorXd*diag ) {
	Eigen::VectorXd aj2 = (((x + v / sigma).cwiseAbs() - w / sigma).array() <= 0.0).cast<double>().matrix();

	//cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;

	double a = min(pow(rhs.norm(), 1.5) * 0.5, 20.0);
	//cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;
	*diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient((*diag + sigma*aj2 + a * Eigen::VectorXd::Ones(var_num)).cwiseSqrt());
	//*diag = Eigen::VectorXd::Ones(var_num);
	MatrixReplacementd A(aj1, n, a, &aj2, sigma,diag);
	//cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

	//Eigen::BiCGSTAB<MatrixReplacement,  Eigen::IdentityPreconditioner> cg;
	Eigen::ConjugateGradient<MatrixReplacementd, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
	double ita = std::fmin(0.5, rhs.norm());
	cg.setTolerance(ita);
	cg.compute(A);
	
	Eigen::VectorXd direction = cg.solve(-diag->cwiseProduct(rhs));
	direction = diag->cwiseProduct(direction);
	cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
	prt_res << cg.iterations() << "\t";
	//cout << "pro_res:" << (A * direction+rhs).norm() << endl;
	//cout << "CG:" << cg.iterations() << std::endl;
	//cout << "di:" << direction.norm() << endl;
	return direction;
}
Eigen::VectorXd norm_ori::get_direction(Eigen::SparseMatrix<double>*ata, Eigen::VectorXd rhs, Eigen::VectorXd* diag) {

	//cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;
	Eigen::VectorXd aj2 = (((x + v / sigma).cwiseAbs() - w / sigma).array() <= 0.0).cast<double>().matrix();
	double a = min(pow(rhs.norm(), 1.5) * 0.5, 20.0);
	(*ata) += sigma*aj2.asDiagonal();
	(*ata) += a * Eigen::VectorXd::Ones(var_num).asDiagonal();
	//cout << (*diag)[0] << "\t" << (*diag)[1] << "\t" << (*diag)[2] << endl;

	Eigen::ConjugateGradient<SparseMatrix<double>, Eigen::Lower | Eigen::Upper, Eigen::DiagonalPreconditioner<double> > cg;
	double ita = std::fmin(0.5, rhs.norm());
	//cg.setTolerance(ita);
	cg.compute(*ata);

	Eigen::VectorXd direction = cg.solve(-rhs);

	cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
	prt_res << cg.iterations() << "\t";
	//cout << "pro_res:" << (A * direction+rhs).norm() << endl;
	//cout << "CG:" << cg.iterations() << std::endl;
	//cout << "di:" << direction.norm() << endl;
	return direction;
}

void norm_ori::run() {
	beta = b.norm();
	b = b / beta;
	string pwd;
	cout << "plz input pwd:";
	cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();

	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 1e-6;
	double s_a = 3;
	sigma = 1 / (s_a * s_a);
	for (; iter_all < 100; iter_all++, sigma *= 2) {


		if (pwd!="skip"&&iter_all > 0) {
			ofstream model(string("model_save/") +pwd+string("_")+graph_name+string("_")+to_string(iter_all)+string(".dat"),std::ios::out|std::ios::binary);
			auto i_end = std::chrono::high_resolution_clock::now();
			auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
			double tm = double(eld.count()) * 1.0 * 1e-9;
			model << setprecision(20);
			model << graph_name << endl;
			model << iter_all << endl;
			model << n << "\t" << var_num <<"\t"<<sigma<<"\t"<<tm+time_record<<endl;



			for (int i = 0; i < var_num; i++)
				model << x[i] << endl;
			for (int i = 0; i < var_num; i++)
				model << v[i] << endl;
			model << u.size() << endl;
			for (auto i = u.begin(); i != u.end();i++) {
				model << i->first << " " << i->second << endl;
			}
		}

		cout << "start main loop,sigma=" << sigma << endl;
		axbu = new axbu_eigen(n, &x, &b, &u, sigma);
		bool flaga = true;

		for (int iter_1 = 0; iter_1 < 1000; iter_1++) {


			cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
			cout << "x:" << x.norm() << endl;
			rhs = Eigen::VectorXd::Zero(var_num);
			vector<int64_t>aj1;
			
			Eigen::VectorXd diag = Eigen::VectorXd::Zero(var_num);
			std::vector<Eigen::Triplet<double> >test_v;
			if (!flaga) {
				if (mthread == 0) {
					axbu->set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
					(*axbu)();
				}
				else {
					//m_calc_rhs_and_aj1(mthread, &rhs, &aj1, n, *axbu);
					m_calc_rhs_and_aj1(mthread, &rhs, &aj1,n, *axbu);
				}
			}
			else {


				if (mthread == 0) {
					axbu->set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
					(*axbu)();
				}
				else {
					//m_calc_rhs_and_aj1(mthread, &rhs, &aj1, n, *axbu);
					m_calc_rhs_and_aj1_diag(mthread, &rhs, &aj1, &diag, n, *axbu,&test_v);
				}
			}

			
			std::vector<Eigen::Triplet<double> >test_vb;
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
			//	cout << test_v[i].row()<<" "<<test_v[i].col()<<" "<<test_v[i].value() << endl;
			//}
			//for (int i = 0; i < 10; i++) {
			//	cout << test_vb[i].row() << " " << test_vb[i].col() << " " << test_vb[i].value() << endl;
			//}
			//return;

			Eigen::SparseMatrix<double>ata(var_num, var_num);
			ata.setFromTriplets(test_v.begin(), test_v.end());

			ata += diag.asDiagonal();
			ata /= 9.0;
			ata *= sigma;
			//cout << diag[0] << "\t" << diag[1] << "\t" << diag[2] << endl;
			//diag = Eigen::VectorXd::Ones(var_num).cwiseQuotient(diag);


			//cout << "r1:" << rhs.norm() << endl;



			rhs = sigma * rhs + sigma * prox_p_star(x+v/sigma,sigma);

			prt_res << iter_1<<"  "<<rhs.norm() << "  ";


			if (rhs.norm() < rlim) {
				break;
			}
			cout << "rhs:" << rhs.norm() << endl;
			Eigen::VectorXd direction;
			//if (flaga) {
			//	direction = get_direction(sigma, rhs, &aj1, &diag);
			//}
			//else {
			//	direction= get_direction(sigma, rhs, &aj1);
			//}
			direction = get_direction(&ata, rhs, &diag);
			//Eigen::VectorXd direction = get_direction(rhs);
			//cout << "di:" << direction.norm() << endl;

			//for (int i = 0; i < var_num; i++)
			//	direction[i] = double(ceil(direction[i] * 1e12)) / 1e12;

			
			double alpha = get_step_length(direction, sigma,404);
			cout << "alpha=" << alpha << endl;
			if (rhs.norm() < 1e-5 && alpha < 1e-4)
				break;
				//flaga = true;
			x = x + alpha * direction;
			prt_res << alpha << std::endl;

		}

		double k = get_kkt(sigma);
		prt_res << "==================\n";
		if (k < 1e-5)
			break;
		//rlim = k * 1e-4;
		rlim = k * 1e-1;

		delete axbu;

	}
	cout << "func_value:" << std::setprecision(8)<< func_p(x)*beta << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	cout << "time is:" << tm << endl;
}
