#include"metric.h"
#include"saveload.h"


extern inline int ij(int x, int y, int n);



void metric::matrix_generate(char type)
{

	int tmp_it = 0;
	if (type == 'r') {

		std::ofstream matrix("matrix.3", std::ios::out);
		matrix << n << " " << n << endl;
		for (int i = 0; i < n; i++)
		{

			for (int j = 0; j < n; j++, tmp_it++) {
				if (i == j) {
					matrix << 0 << " ";
				}
				else {
					matrix << b[ij(i, j, n)] << " ";
				}
				//calculate d
			}
			matrix << endl;
		}
	}
	else if (type == 'o') {
		std::ofstream matrix("matrixo.3", std::ios::out);
		for (int i = 0; i < n; i++)
		{

			for (int j = 0; j < n; j++, tmp_it++) {
				if (i == j) {
					matrix << 0 << " ";
				}
				else {
					matrix << fmax(x[ij(i,j,n)]+b[ij(i, j, n)],0) << " ";
				}
				//calculate d
			}
			matrix << endl;
		}
	}
}

void saveHashB(ofstream&fout,metric::hash_b&b) {
	save4(&b.default_value, fout);
	saveSet(b.data, fout);
}
void readHashB(ifstream& fin, metric::hash_b& b) {
	read4(&b.default_value, fin);
	readSet(b.data, fin);
}

void saveHashX(ofstream& fout, metric::hash_x& x) {
	save8(&x.default_value, fout);
	saveMap(x.data, fout);
}
void readHashX(ifstream& fin, metric::hash_x& x) {
	read8(&x.default_value, fin);
	readMap(x.data, fin);
}

void metric::SparseSave(ofstream& fout)
{
	save4(&n, fout);
	saveHashB(fout,sv_b);
	saveHashX(fout, sv_x);
	saveHashX(fout, sv_v);
	saveHashX(fout, sv_w);
}
metric::metric(string graph_name, ifstream& fin) :sv_x(sparse_w, 0), sv_v(sparse_w, 0), sv_w(sparse_w, 0), sv_b(sparse_b) {
	read4(&n, fin);
	cout << "n:" << n << endl;
	var_num = static_cast<int>(static_cast<int64_t>(n) * (n - 1) / 2);
	readHashB(fin,sv_b);
	readHashX(fin, sv_x);
	readHashX(fin, sv_v);
	readHashX(fin, sv_w);
	double wrs = 0;
	for (int i = 0; i < var_num; i++)wrs += sv_w[i];
	cout << "w:" << wrs << "\nb:" << sparse_b.size() << "\nwdef:" << sv_w.default_value << endl;

}


double calc_bnorm(int mthread, Eigen::VectorXd& b, int n) {

	auto calc = [&](int ia, int ib, double& res) {
		double nres = 0;
		for (int i = ia; i < ib; i++) {
			for (int j = i + 1; j < n; j++) {
				for (int k = j + 1; k < n; k++) {
					int idx[3];
					idx[0] = ij(i, j, n);
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);
					int tp = floor(b[idx[0]] + b[idx[1]] + b[idx[2]]);
					if (tp == 3) {
						nres += 3;
					}
					else if (tp == 2) {
						nres += 4;
					}
					else if (tp == 1) {
						nres += 3;
					}

				}
			}
		}
		res = nres;

	};
	vector<thread>vt;
	vector<double>res(mthread, 0.0);
	double finalres = 0;
	for (int i = 0; i < mthread; i++) {
		vt.push_back(thread(calc, n - (int)(cbrt((i + 1) * 1.0 / mthread) * n), n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n), std::ref(res[i])));
	}
	for (int i = 0; i < mthread; i++) {
		vt[i].join();
		finalres += res[i];
	}
	return finalres;
}