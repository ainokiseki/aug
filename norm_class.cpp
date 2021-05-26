
#include"norm_class.h"
#include<iomanip>

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
void m_calc_rhs_mat(int pnum, Eigen::VectorXd* rhs,int n, axbu_aj1 axbu, Eigen::SparseMatrix<double>*mat) {
	int var_num = n * (n - 1) / 2;
	Eigen::VectorXd* mrhs = new Eigen::VectorXd[pnum];

	vector<Eigen::SparseMatrix<double>* >mmat(pnum);
	for (int i = 0; i < pnum; i++) {

		mmat[i] = new Eigen::SparseMatrix<double>(var_num, var_num);
	}
	//vector<int64_t>* maj1 = new vector<int64_t>[pnum];
	for (int i = 0; i < pnum; i++) {
		mrhs[i] = Eigen::VectorXd::Zero(var_num);
	}

	vector<thread>tv;
	for (int i = 0; i < pnum; i++) {
		axbu.mat = mmat[i];
		axbu.rhs = mrhs+i;
		tv.push_back(thread(axbu, n - (int)(cbrt((i + 1) * 1.0 / pnum) * n), n - (int)(cbrt(i * 1.0 / pnum) * 1.0 * n)));
	}
	for (int i = 0; i < pnum; i++)
		tv[i].join();
	for (int i = 0; i < pnum; i++) {
		*rhs += mrhs[i];
	}


	delete[]mrhs;
	for (int i = 0; i < pnum; i++) {
		(*mat) += (*mmat[i]);
		delete mmat[i];
	}


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
norm_basic::norm_basic(string graph_name) :metric(graph_name)
{
}

norm_basic::~norm_basic()
{
}



double norm_basic::get_kkt(double sigma) {
	double y_norm = 0;
	double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
	kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
	axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
	(*axbu)();
	cout << "u:" << u.size() << endl;
	Eigen::VectorXd z = w.cwiseProduct(x) + v / sigma - pi_c(w.cwiseProduct(x) + v / sigma, sigma);
	v = sigma * pi_c(w.cwiseProduct(x) + v / sigma, sigma);

	y_norm = sqrt(y_norm);
	kkt_u = sqrt(kkt_u) / (1 + y_norm);
	kkt_y = sqrt(kkt_y) / (1 + y_norm);
	cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
	kkt_x = au_eigen(n, &u, v, w)();
	kkt_x = kkt_x / (1 + v.norm());

	//cout << "z_norm" << z.norm() << endl;
	kkt_v = (w.cwiseProduct(x) - z).norm() / (1 + z.norm());
	kkt_z = (v - pi_c(z + v, 1)).norm() / (1 + z.norm());
	//for (int i = 0; i < var_num; i++)
	//    t_res += abs(x[i] * w[i]);
	//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
	double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
	cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
	cout << "kkt:" << kkt_final << endl;


	return kkt_final;

	//cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;




}

double norm_basic::get_step_length(Eigen::VectorXd& direction, double sigma) {

	//cg end,result saved at direction
	//find step length
	double alpha = 1;
	const double c1 = 1e-4;
	double lx0 = 0;
	axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
	(*axbu)();

	auto square = [](double x) {return x * x; };


	lx0 = sigma / 2 * lx0 + sigma / 2 * square(pi_c(w.cwiseProduct(x) + v / sigma, sigma).norm());

	lx0 += func_p(w.cwiseProduct(x) + v / sigma - pi_c(w.cwiseProduct(x) + v / sigma, sigma));

	double down = lx0 + c1 * alpha * rhs.dot(direction);
	//cout << "down:" << down << "\tlx0:" << lx0 << endl;

	cout << "lx0:" << lx0 << "\tdown" << down << endl;

	Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(var_num);
	axbu->set_function(new func_eigen_rhs(alpha, &direction, &nrhs));
	(*axbu)();
	//m_calc_rhs(8, nrhs, n, axbu, alpha, &direction);
	nrhs = sigma * nrhs;
	nrhs += sigma * w.cwiseProduct(pi_c(w.cwiseProduct(x + alpha * direction) + v / sigma, sigma));
	if (nrhs.dot(direction) < 0)
		return 1.0;

	for (int i = 0; i < 20; i++) {
		//cout << "start find step,alpha=" << alpha << endl;
		double lx = 0;
		axbu->set_function(new func_eigen_phi(alpha, &direction, &lx));
		(*axbu)();
		//m_calc_phi(8, lx, n, axbu, alpha, &direction);
		lx = sigma / 2 * lx;
		lx += sigma / 2 * square(pi_c(w.cwiseProduct(x + alpha * direction) + v / sigma, sigma).norm());
		lx += func_p(w.cwiseProduct(x + alpha * direction) + v / sigma - pi_c(w.cwiseProduct(x + alpha * direction) + v / sigma, sigma));
		if (lx <= down) {
			break;
		}
		alpha *= 0.5;
	}
	return alpha;


}

void norm_basic::run() {
	init();
	double sigma = 1;

	for (int iter_all = 0; iter_all < 1000; iter_all++, sigma *= 2) {


		cout << "start main loop,sigma=" << sigma << endl;
		axbu = new axbu_eigen(n, &x, &b, &u, sigma);
		for (int iter_1 = 0; iter_1 < 150; iter_1++) {



			cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
			rhs = Eigen::VectorXd::Zero(var_num);
			vector<int64_t>aj1;
			axbu->set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
			(*axbu)();

			rhs = sigma * rhs + sigma * w.cwiseProduct(pi_c(w.cwiseProduct(x) + v / sigma, sigma));
			if (rhs.norm() < 1e-6) {

				break;
			}
			cout << "rhs:" << rhs.norm() << endl;
			Eigen::VectorXd direction = get_direction(sigma, rhs, &aj1);
			double alpha = get_step_length(direction, sigma);
			cout << "alpha=" << alpha << endl;
			x = x + alpha * direction;

		}
		double k = get_kkt(sigma);

		if (k < 1e-5)
			break;

		delete axbu;

	}

	cout << "func_value:"<<func_p(w.cwiseProduct(x)) << endl;



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
	axbu->set_function(new func_eigen_u(&u, sigma));
	(*axbu)();

	axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
	(*axbu)();

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
	cg.setTolerance(ita);
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
	string pwd;
	cout << "plz input pwd:";
	cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();

	prt_res = ofstream("123.log");
	prt_res.setf(std::ios::scientific);
	//double sigma = 1;
	double rlim = 1e-6;

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
			axbu_aj1 axbumat(n,&x,&b,&u,sigma,NULL,NULL);
			Eigen::SparseMatrix<double>ata(var_num, var_num);
			m_calc_rhs_mat(mthread, &rhs, n, axbumat, &ata);
			/*
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

			*/
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


			//ata.setFromTriplets(test_v.begin(), test_v.end());

			//ata += diag.asDiagonal();
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


		delete axbu;

	}
	cout << "func_value:" << std::setprecision(8)<< func_p(x) << endl;
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	cout << "time is:" << tm << endl;
}
norm_1::norm_1(string graph_name) :norm_basic(graph_name)
{
}

norm_1::~norm_1()
{
}

Eigen::VectorXd norm_1::get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1) {
	Eigen::VectorXd aj2 = ((w.cwiseProduct(x) + v / sigma).cwiseAbs().array() < 1 / sigma).matrix().cast<double>();
	aj2 = aj2.cwiseProduct(w.cwiseProduct(w));

	//cout << "rnorm:" << r_norm << "\trnew:" << r_new << endl;

	double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);

	MatrixReplacement A(aj1, var_num, a, &aj2, sigma);
	//Eigen::BiCGSTAB<MatrixReplacement,  Eigen::IdentityPreconditioner> cg;
	Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
	cg.compute(A);
	Eigen::VectorXd direction = cg.solve(-rhs);
	//cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;
	cout << "CG:" << cg.iterations() << std::endl;

	return direction;
}

double norm_1::func_p(Eigen::VectorXd x) {
	return x.cwiseAbs().sum();
}

inline Eigen::VectorXd norm_1::pi_c(Eigen::VectorXd x, double sigma) {
	return x.cwiseSign().cwiseProduct(x.cwiseAbs().cwiseMin(1 / sigma));
}



norm_2::norm_2(string graph_name) :norm_basic(graph_name)
{
}

norm_2::~norm_2()
{
}

Eigen::VectorXd norm_2::get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1) {

	double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);
	Eigen::VectorXd wx = w.cwiseProduct(x) + v / sigma;
	cout << "wx::" << wx.norm() << "\t" << 1 / sigma << endl;

	MatrixReplacement2 A(aj1, var_num, a, &wx, &w, sigma);

	Eigen::ConjugateGradient<MatrixReplacement2, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
	cg.setTolerance(1e-10);
	cg.compute(A);


	Eigen::VectorXd direction = cg.solve(-rhs);

	cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;


	return direction;
}

double norm_2::func_p(Eigen::VectorXd x) {
	return x.norm();
}

inline Eigen::VectorXd norm_2::pi_c(Eigen::VectorXd x, double sigma) {
	if (x.norm() <= 1 / sigma)
		return x;
	return x * (1 / sigma) / x.norm();
}
void norm_2::init() {
	norm_basic::init();
	for (int i = 0; i < var_num; i++)
		w[i] = sqrt(w[i]);
	cout << w.sum() << endl;
}



norm_inf::norm_inf(string graph_name) :norm_basic(graph_name)
{
}

norm_inf::~norm_inf()
{
}

Eigen::VectorXd norm_inf::get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1) {

	double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);
	Eigen::VectorXd wx = w.cwiseProduct(x) + v / sigma;
	//cout << "wx::" << wx.norm() << "\t" << 1 / sigma << endl;
	Eigen::VectorXd r = pi_c(wx,sigma);
	r = Eigen::VectorXd::Ones(var_num)-r.cwiseEqual(0).cast<double>();
	MatrixReplacement_inf A(aj1, var_num, a, &wx, &w,&r, sigma);
	
	Eigen::ConjugateGradient<MatrixReplacement_inf, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
	cg.setTolerance(1e-12);
	cg.compute(A);


	Eigen::VectorXd direction = cg.solve(-rhs);

	cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;


	return direction;
}

double norm_inf::func_p(Eigen::VectorXd x) {
	return x.cwiseAbs().maxCoeff();
}

inline Eigen::VectorXd norm_inf::pi_c(Eigen::VectorXd x, double sigma) {
	

	Eigen::VectorXd tmp = x.cwiseAbs();
	double gap = tmp.sum() - 1 / sigma;

	if (gap <= 0)
		return x;
	
	sort(tmp.data(), tmp.data() + tmp.size());
	double act;
	int len = tmp.size();
	for (int i = 0; i < len; i++) {

		if ((len - i ) * tmp[i] >= gap) {
			act = gap / (len - i );
			break;
		}
		else {
			gap -= tmp[i];
		}
	}
	return x.cwiseSign().cwiseProduct(x - act * x.cwiseSign()).cwiseMax(0).cwiseProduct(x.cwiseSign());
}



norm_2_if::norm_2_if(string graph_name) :graph_name(graph_name)
{
}

norm_2_if::~norm_2_if()
{
}



void norm_2_if::init() {
	ifstream graph(graph_name, ios::in);
	map<int, set<int> > g;
	int tmp_a, tmp_b;
	int input_ct = 0;
	while (graph >> tmp_a >> tmp_b) {
		g[tmp_a].insert(tmp_b);
		g[tmp_b].insert(tmp_a);
		input_ct++;
	}

	cout << "line:" << input_ct << endl;
	//point list
	vector<int> target;
	//int n= 50;
	for (auto it = g.begin(); it != g.end(); it++) {
		target.push_back(it->first);
	}
	n = target.size();
	cout << "n:" << n << endl;
	var_num = n * (n - 1) / 2;
	cons_num = n * (n - 1) * (n - 2) / 6;
	w = Eigen::VectorXd::Zero(var_num);
	b = Eigen::VectorXd::Zero(var_num);

	int tmp_it = 0;
	const double delta = 0.05;
	const double epsi = 0.01;
	//calculate w
	for (int i = 0; i < n; i++)
	{
		int t_i = target[i];
		for (int j = i + 1; j < n; j++, tmp_it++) {
			//calculate w
			set<int>tmp1;
			set<int>tmp2;
			int t_j = target[j];
			set_intersection(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp1, tmp1.begin()));
			set_union(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp2, tmp2.begin()));
			double tmp = 1.0 * tmp1.size() / tmp2.size();
			tmp = log((1 + tmp - delta) / (1 - tmp + delta));
			if (tmp > 0) {
				tmp += epsi;
			}
			else if (tmp < 0) {
				tmp -= epsi;
			}
			else {
				if (g[t_i].find(t_j) != g[t_i].end())
					tmp = epsi;
				else
					tmp = -epsi;
			}
			b[tmp_it] = (tmp < 0);
			w[tmp_it] = sqrt(abs(tmp));
			//calculate d
		}
	}

	cout << "w:" << w.sum() << endl;
	cout << "b:" << b.sum() << endl;
	x = Eigen::VectorXd::Ones(var_num);

	return;
}

double norm_2_if::get_kkt(double sigma) {
	double y_norm = 0;
	double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
	kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
	axbu->set_function(new func_eigen_y(&u, sigma, &y_norm, &kkt_y, &kkt_u));
	(*axbu)();
	cout << "u:" << u.size() << endl;

	y_norm = sqrt(y_norm);
	kkt_u = sqrt(kkt_u) / (1 + y_norm);
	kkt_y = sqrt(kkt_y) / (1 + y_norm);
	cout << "kkt_u:" << kkt_u << "\tkkt_y:" << kkt_y << endl;
	Eigen::VectorXd v = Eigen::VectorXd::Zero(var_num);
	kkt_x = au_eigen(n, &u, v, w)();
	
	kkt_x = kkt_x / (1 + x.norm());

	//cout << "z_norm" << z.norm() << endl;
	//for (int i = 0; i < var_num; i++)
	//    t_res += abs(x[i] * w[i]);
	//prt_res << "f(Ax-b)+p(x) \t" << t_res << endl;
	double kkt_final = max(max(kkt_u, kkt_y), max(max(kkt_v, kkt_z), kkt_x));
	cout << "kkt_x:" << kkt_x << "\tkkt_z:" << kkt_z << "\tkkt_v:" << kkt_v << endl;
	cout << "kkt:" << kkt_final << endl;


	return kkt_final;

	//cout << "v_norm:" << v.norm() << "y_norm" << y_norm <<"x_norm"<<x.norm()<< endl;




}

double norm_2_if::get_step_length(Eigen::VectorXd& direction, double sigma) {

	//cg end,result saved at direction
	//find step length
	double alpha = 1;
	const double c1 = 1e-4;
	double lx0 = 0;
	axbu->set_function(new func_eigen_phi(0, &direction, &lx0));
	(*axbu)();
	auto square = [](double x) {return x * x; };


	lx0 = sigma / 2 * lx0 ;

	lx0 += func_p(w.cwiseProduct(x));

	double down = lx0 + c1 * alpha * rhs.dot(direction);
	//cout << "down:" << down << "\tlx0:" << lx0 << endl;

	cout << "lx0:" << lx0 << "\tdown" << down << endl;

	Eigen::VectorXd nrhs = Eigen::VectorXd::Zero(var_num);
	axbu->set_function(new func_eigen_rhs(alpha, &direction, &nrhs));
	(*axbu)();
	//m_calc_rhs(8, nrhs, n, axbu, alpha, &direction);
	nrhs = sigma * nrhs;
	nrhs += w.cwiseProduct(w).cwiseProduct(x+alpha*direction)/(w.cwiseProduct(x+alpha*direction).norm());
	if (nrhs.dot(direction) < 0)
		return 1.0;

	for (int i = 0; i < 20; i++) {
		//cout << "start find step,alpha=" << alpha << endl;
		double lx = 0;
		axbu->set_function(new func_eigen_phi(alpha, &direction, &lx));
		(*axbu)();
		//m_calc_phi(8, lx, n, axbu, alpha, &direction);
		lx = sigma / 2 * lx;
		lx += func_p(w.cwiseProduct(x+alpha*direction));
		if (lx <= down) {
			break;
		}
		alpha *= 0.5;
	}
	return alpha;


}

void norm_2_if::run() {
	init();
	double sigma = 1;
	for (int iter_all = 0; iter_all < 1000; iter_all++, sigma *= 2) {
		cout << "start main loop,sigma=" << sigma << endl;
		axbu = new axbu_eigen(n, &x, &b, &u, sigma);
		for (int iter_1 = 0; iter_1 < 150; iter_1++) {
			cout << "start find minimum of phi(x),iter time:" << iter_1 << endl;
			rhs = Eigen::VectorXd::Zero(var_num);
			vector<int64_t>aj1;
			axbu->set_function(new func_eigen_rhs_aj1(&rhs, &aj1));
			(*axbu)();
			rhs = sigma * rhs + w.cwiseProduct(w).cwiseProduct(x)/(w.cwiseProduct(x).norm());
			if (rhs.norm() < 1e-6) {
				break;
			}
			cout << "rhs:" << rhs.norm() << endl;
			Eigen::VectorXd direction = get_direction(sigma, rhs, &aj1);
			double alpha = get_step_length(direction, sigma);
			cout << "alpha=" << alpha << endl;
			x = x + alpha * direction;
		}
		double k = get_kkt(sigma);

		if (k < 1e-5)
			break;

		delete axbu;

	}
	cout << "func_value:" << func_p(w.cwiseProduct(x)) << endl;



}

Eigen::VectorXd norm_2_if::get_direction(double sigma, Eigen::VectorXd rhs, vector<int64_t>* aj1) {

	double a = min(pow(rhs.norm(), 1.5) * 0.5, 0.5);
	Eigen::VectorXd wx = w.cwiseProduct(x);
	cout << "wx::" << wx.norm() << "\t" << 1 / sigma << endl;

	MatrixReplacement3 A(aj1, var_num, a, &wx, &w, 1e9);

	Eigen::ConjugateGradient<MatrixReplacement3, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
	cg.setTolerance(1e-12);
	cg.compute(A);


	Eigen::VectorXd direction = cg.solve(-rhs);

	cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;


	return direction;
}

double norm_2_if::func_p(Eigen::VectorXd x) {
	return x.norm();
}
