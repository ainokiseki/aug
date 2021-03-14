#include"admm_class.h"

admm_basic::admm_basic(string graph_name):metric(graph_name)
{
	
}

admm_basic::~admm_basic()
{
}

double admm_basic::get_kkt(double sigma) {

	double y_norm = 0;
	double kkt_y, kkt_u, kkt_x, kkt_z, kkt_v;
	kkt_x = kkt_y = kkt_u = kkt_v = kkt_z = 0;
	axbu->set_function(new func_eigen_u(&u, sigma, &y_norm, &kkt_y, &kkt_u));
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

void admm_basic::run() {

	init();
	double sigma = 1e-2;
	for (int iter_all = 0; iter_all < 100000; iter_all++,sigma=1e-2/(sqrt(1+iter_all))) {
		//cout << "start main loop,sigma=" << sigma << endl;
		axbu_pro* axbp = new axbu_pro(n, &x, &b, &u, sigma);
		axbu = new axbu_eigen(n, &x, &b, &u, sigma);
		Eigen::VectorXd rhs = (*axbp)();
		//cout << "rhs_sum:" << rhs.norm() << endl;
		Eigen::VectorXd z = w.cwiseProduct(x) + v / sigma - pi_c(w.cwiseProduct(x) + v / sigma, sigma);
		rhs += w.cwiseProduct(z - v / sigma);
		//cout << "rhs_sum:" << rhs.norm() << endl;
		MatrixReplacement_admm A(n, &w);

		Eigen::ConjugateGradient<MatrixReplacement_admm, Eigen::Lower | Eigen::Upper, Eigen::IdentityPreconditioner> cg;
		cg.setTolerance(1e-12);
		cg.compute(A);
		x = cg.solve(rhs);
		//cout << "xs:" << x.size() << endl;

		//cout << "CG:       #iterations: " << cg.iterations() << ", estimated error: " << cg.error() << std::endl;

		double k = get_kkt(sigma);

		if (k < 1e-5)
			break;

		delete axbu;
		delete axbp;

	}
	cout << "func_value:" << func_p(w.cwiseProduct(x)) << endl;






}

admm_1::admm_1(string graph_name):admm_basic(graph_name)
{
}

admm_1::~admm_1()
{
}

inline Eigen::VectorXd admm_1::pi_c(Eigen::VectorXd x, double sigma) {
	return x.cwiseSign().cwiseProduct(x.cwiseAbs().cwiseMin(1 / sigma));
}

double admm_1::func_p(Eigen::VectorXd x) {
	return x.cwiseAbs().sum();
}
