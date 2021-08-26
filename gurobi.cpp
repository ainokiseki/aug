#include"gurobi.h"
#include<chrono>
#include"mytime.h"
#include"tri.h"
#include<thread>
#include<mutex>
void cons_generate(int start, int end, int n, Eigen::VectorXd& x, Eigen::VectorXd& b, std::mutex* mu, std::unordered_map<int, int>& colidx, std::set<int64_t>& viscons)
{

	for (int i = start; i < end; i++) {
		int64_t index;
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				index = static_cast<int64_t>(i) * n * n + j * n + k;

				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 1e-4) {

					if (viscons.find(index) == viscons.end()) {

						for (int tp = 0; tp < 3; tp++) {
							mu[1].lock();
							if (colidx.find(idx[tp]) == colidx.end()) {

								int tmmp = colidx.size();
								colidx[idx[tp]] = tmmp;

							}
							mu[1].unlock();
						}

						mu[0].lock();

						viscons.insert(index);
						mu[0].unlock();
					}
				}
			}
		}
	}


}
gurobi::gurobi(string graphname) :metric(graphname)
{

}

gurobi::~gurobi()
{
}

void gurobi::run()
{

	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();



	double rlim = 1e-4;

	int mthread = 24;
	std::vector<int64_t>cons;
	std::set<int64_t>viscons;
	double c_t = 0;
	std::unordered_map<int, int>colidx;
	int old_size = viscons.size();
	timer time1, time2, time3, time4;
	while (true) {
		auto c_start = std::chrono::high_resolution_clock::now();

		if (false) {



			for (int i = 0; i < n; i++) {
				int64_t index;
				for (int j = i + 1; j < n; j++) {
					for (int k = j + 1; k < n; k++) {
						int idx[3];
						idx[0] = ij(i, j, n);
						idx[1] = ij(i, k, n);
						idx[2] = ij(j, k, n);
						index = static_cast<int64_t>(i) * n * n + j * n + k;

						if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
							x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
							x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 1e-4) {

							if (viscons.find(index) == viscons.end()) {
								for (int tp = 0; tp < 3; tp++) {
									if (colidx.find(idx[tp]) == colidx.end()) {
										int tmmp = colidx.size();
										colidx[idx[tp]] = tmmp;
									}
								}

								viscons.insert(index);

							}
						}
					}
				}
			}

		}
		else {
			std::vector<thread>vt;
			mutex mu[2];
			for (int i = 0; i < mthread; i++) {
				int ta, tb;
				ta = n - (int)(cbrt((i + 1) * 1.0 / mthread) * n);
				tb = n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n);
				vt.push_back(std::thread(cons_generate, ta, tb, n, std::ref(x), std::ref(b), mu, std::ref(colidx), std::ref(viscons)));
			}
			for (int i = 0; i < mthread; i++) {
				vt[i].join();
			}
		}

		if (old_size == viscons.size())
			break;
		old_size = viscons.size();
		std::cout << "size:" << cons.size() << "vis:" << viscons.size() << endl;

		if (viscons.size() < 0 * var_num / 3) {
			for (int i = 0; i < n && viscons.size() < 10 * var_num / 3; i++) {
				int64_t index;
				for (int j = i + 1; j < n && viscons.size() < 10 * var_num / 3; j++) {
					for (int k = j + 1; k < n && viscons.size() < 10 * var_num / 3; k++) {

						index = static_cast<int64_t>(i) * n * n + j * n + k;



						if (viscons.find(index) == viscons.end()) {
							int idx[3];
							idx[0] = ij(i, j, n);
							idx[1] = ij(i, k, n);
							idx[2] = ij(j, k, n);

							for (int tp = 0; tp < 3; tp++) {
								if (colidx.find(idx[tp]) == colidx.end()) {
									int tmmp = colidx.size();
									colidx[idx[tp]] = tmmp;
								}
							}

							viscons.insert(index);
						}


					}
				}
			}
		}




		std::cout << "size:" << cons.size() << "vis:" << viscons.size() << endl;

		std::cout << "total num of variables:" << var_num << " num of current variables:" << colidx.size() << endl;

		auto c_end = std::chrono::high_resolution_clock::now();
		auto c_eld = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start);
		c_t += double(c_eld.count()) * 1.0 * 1e-9;

		n1(viscons, colidx);

	}

}

void gurobi::n1(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx)
{

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar* y = model.addVars(colidx.size(), GRB_CONTINUOUS);
	GRBVar* m = model.addVars(colidx.size(), GRB_CONTINUOUS);
	try
	{
		GRBLinExpr obj = 0;


		for (auto i : colidx) {
			y[i.second].set(GRB_DoubleAttr_Start, x[i.first] + b[i.first]);
			m[i.second].set(GRB_DoubleAttr_Start, fabs(x[i.first]));

			obj += m[i.second] * w[i.first];


			//model.addGenConstrAbs(m[i.second], y[i.second] - b[i.first]);
			model.addConstr(m[i.second] - y[i.second] + b[i.first] >= 0);
			model.addConstr(m[i.second] + y[i.second] - b[i.first] >= 0);
		}


		for (auto index : viscons) {
			int i, j, k;
			k = index % n;
			j = index / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = colidx[ij(i, j, n)];
			idx[1] = colidx[ij(i, k, n)];
			idx[2] = colidx[ij(j, k, n)];
			model.addConstr(y[idx[0]] - y[idx[1]] - y[idx[2]] <= 0);
			model.addConstr(y[idx[1]] - y[idx[0]] - y[idx[2]] <= 0);
			model.addConstr(y[idx[2]] - y[idx[1]] - y[idx[0]] <= 0);
		}
		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		for (auto i : colidx) {
			x[i.first] = y[i.second].get(GRB_DoubleAttr_X) - b[i.first];
		}
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}

	delete[] y;
	delete[] m;
}

void gurobi::n2(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx)
{
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar* y = model.addVars(viscons.size(), GRB_CONTINUOUS);
	try
	{
		GRBQuadExpr obj = 0;


		for (auto i : colidx) {
			y[i.second].set(GRB_DoubleAttr_Start, x[i.first] + b[i.first]);

			obj += w[i.first]*(y[i.second]-b[i.first])*(y[i.second] - b[i.first]);


			//model.addGenConstrAbs(m[i.second], y[i.second] - b[i.first]);
		}


		for (auto index : viscons) {
			int i, j, k;
			k = index % n;
			j = index / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = colidx[ij(i, j, n)];
			idx[1] = colidx[ij(i, k, n)];
			idx[2] = colidx[ij(j, k, n)];
			model.addConstr(y[idx[0]] - y[idx[1]] - y[idx[2]] <= 0);
			model.addConstr(y[idx[1]] - y[idx[0]] - y[idx[2]] <= 0);
			model.addConstr(y[idx[2]] - y[idx[1]] - y[idx[0]] <= 0);
		}
		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		for (auto i : colidx) {
			x[i.first] = y[i.second].get(GRB_DoubleAttr_X) - b[i.first];
		}
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}

	delete[] y;

}
