#include"gurobi.h"
#include<chrono>
#include"mytime.h"
#include"tri.h"
#include<thread>
#include<mutex>
#include"cons_generate.h"
#include"norm_class.h"
void cons_generate(int start, int end, int n,double&kkt1,double&kkt2, Eigen::VectorXd& x, Eigen::VectorXd& b, std::mutex* mu, std::unordered_map<int, int>& colidx, std::set<int64_t>& viscons)
{
	double sk1 = 0, sk2 = 0;
	for (int i = start; i < end; i++) {
		int64_t index;
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				index = static_cast<int64_t>(i) * n * n + j * n + k;
				double tol = x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]];
				if (tol > 0)sk1 += tol * tol;
				else sk2 += tol * tol;
				tol = x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]];
				if (tol > 0)sk1 += tol * tol;
				else sk2 += tol * tol;
				tol = x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]];
				if (tol > 0)sk1 += tol * tol;
				else sk2 += tol * tol;

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
	kkt1 = sk1; kkt2 = sk2;


}
gurobi::gurobi(string graphname,int type) :metric(graphname,type)
{
	kkt_value[0] = kkt_value[1] = 0;
}

gurobi::~gurobi()
{
}

double gurobi::run(int type,int maxstep)
{
	cons_generator generator(n,32,cgflg);
	cout << generator.cgflag << endl;
	string pwd("skip");
	//std::cout << "plz input pwd:";
	//cin >> pwd;
	auto i_start = std::chrono::high_resolution_clock::now();

	double rlim = 1e-4;

	int mthread = 24;
	//std::vector<int64_t>cons;
	//std::set<int64_t>viscons;
	double c_t = 0;
	double bnorm = calc_bnorm(mthread, b, n);
	int lastcons = 0;
	timer time1, time2, time3, time4;
	double answer=0;
	bool initflag = true;
	for (int iter = 0; iter < maxstep;iter++) {
#ifdef FULLCONS
		if (!initflag) {
			
			return answer;
		}
#endif
		auto c_start = std::chrono::high_resolution_clock::now();
		double kt1 = 0, total_ax=0,kt2=0;
		generator.generate(x,b,1e-4);
		auto& viscons(generator.viscons);
		auto& colidx(generator.colidx);
		int colidx_num = generator.colidx_num;
		//if (false) {

		//	for (int i = 0; i < n; i++) {
		//		int64_t index;
		//		for (int j = i + 1; j < n; j++) {
		//			for (int k = j + 1; k < n; k++) {
		//				int idx[3];
		//				idx[0] = ij(i, j, n);
		//				idx[1] = ij(i, k, n);
		//				idx[2] = ij(j, k, n);
		//				index = static_cast<int64_t>(i) * n * n + j * n + k;

		//				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 1e-3 ||
		//					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 1e-3 ||
		//					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 1e-3) {

		//					if (viscons.find(index) == viscons.end()) {
		//						for (int tp = 0; tp < 3; tp++) {
		//							if (colidx.find(idx[tp]) == colidx.end()) {
		//								int tmmp = colidx.size();
		//								colidx[idx[tp]] = tmmp;
		//							}
		//						}

		//						viscons.insert(index);

		//					}
		//				}
		//			}
		//		}
		//	}

		//}
		//else {
		//	std::vector<thread>vt;
		//	std::vector<double>vkkt1(mthread, 0), vkkt2(mthread, 0);
		//	mutex mu[2];
		//	for (int i = 0; i < mthread; i++) {
		//		int ta, tb;
		//		ta = n - (int)(cbrt((i + 1) * 1.0 / mthread) * n);
		//		tb = n - (int)(cbrt(i * 1.0 / mthread) * 1.0 * n);
		//		vt.push_back(std::thread(cons_generate, ta, tb, n, std::ref(vkkt1[i]),std::ref(vkkt2[i]),std::ref(x), std::ref(b), mu, std::ref(colidx), std::ref(viscons)));
		//	}
		//	for (int i = 0; i < mthread; i++) {
		//		vt[i].join();
		//		kt1 += vkkt1[i];
		//		kt2 += vkkt2[i];
		//	}
		//}
		if (generator.consvio) {
			return answer;

		}
		if (!initflag) {
#ifdef NEWFM
			double kkt = sqrt(kt1 + kkt_value[0]) / (1 +  sqrt(kkt_value[1]) +sqrt(total_ax) +sqrt(bnorm));
#else
			cout << kkt_value[0] << "\t" << kkt_value[1] << "\t" << kt1 << "\t" << kt2 << endl;
			double kkt = sqrt(kt1 + kkt_value[0]) / (1 + sqrt(kt2+kkt_value[1]));
			//double kkt = sqrt(kt1 + kkt_value[0]);
#endif
			cout << "totalkkt:" << kkt << endl;
			if (kkt< 1e-3&&maxstep==9999) {
				return answer;
			}
		}
		else {
			initflag = false;
		}
		//if (old_size == viscons.size())
		//	break;
		//old_size = viscons.size();
		//std::cout << "size:" << cons.size() << "vis:" << viscons.size() << endl;

		//if (viscons.size() < 0 * var_num / 3) {
		//	for (int i = 0; i < n && viscons.size() < 10 * var_num / 3; i++) {
		//		int64_t index;
		//		for (int j = i + 1; j < n && viscons.size() < 10 * var_num / 3; j++) {
		//			for (int k = j + 1; k < n && viscons.size() < 10 * var_num / 3; k++) {

		//				index = static_cast<int64_t>(i) * n * n + j * n + k;



		//				if (viscons.find(index) == viscons.end()) {
		//					int idx[3];
		//					idx[0] = ij(i, j, n);
		//					idx[1] = ij(i, k, n);
		//					idx[2] = ij(j, k, n);

		//					for (int tp = 0; tp < 3; tp++) {
		//						if (colidx.find(idx[tp]) == colidx.end()) {
		//							int tmmp = colidx.size();
		//							colidx[idx[tp]] = tmmp;
		//						}
		//					}

		//					viscons.insert(index);
		//				}


		//			}
		//		}
		//	}
		//}


		lastcons = viscons.size();
		std::cout << "consnum:" << viscons.size() << endl;

		std::cout << "total num of variables:" << var_num << " num of current variables:" << colidx_num << endl;

		auto c_end = std::chrono::high_resolution_clock::now();
		auto c_eld = std::chrono::duration_cast<std::chrono::nanoseconds>(c_end - c_start);
		c_t += double(c_eld.count()) * 1.0 * 1e-9;
		if (type == 1) {
			answer= n1(viscons, colidx,colidx_num);
		}
		else if (type == 2) {
			answer= n2(viscons, colidx,colidx_num);
		}
		else if (type == 0) {
			answer = n_inf(viscons, colidx, colidx_num);
		}
		
		//savex("grbpower1");
		//return -1;
	}
return answer;
}

double gurobi::n1(std::set<int64_t>& viscons, std::vector<int>& colidx, int colidx_num)
{

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar* y = model.addVars(colidx_num, GRB_CONTINUOUS);
	GRBVar* m = model.addVars(colidx_num, GRB_CONTINUOUS);
	try
	{
		GRBLinExpr obj = 0;

		int64_t conscount = 0;
		for (int i = 0; i < var_num;i++) {
			if (colidx[i] != -1) {


				y[colidx[i]].set(GRB_DoubleAttr_Start, x[i] + b[i]);
				m[colidx[i]].set(GRB_DoubleAttr_Start, fabs(x[i]));

				obj += m[colidx[i]] * w[i];


				//model.addGenConstrAbs(m[i.second], y[i.second] - b[i.first]);
				model.addConstr(-m[colidx[i]] + y[colidx[i]] - b[i] <= 0,string("p")+to_string(i));
				model.addConstr(-m[colidx[i]] - y[colidx[i]] + b[i] <= 0, string("n") + to_string(i));
				conscount += 2;
			}
		}






		add_constr_A(model, viscons, colidx, y);
		//model.update();




		model.set(GRB_IntParam_Method,2);
		model.set(GRB_DoubleParam_FeasibilityTol, 0.01);
		model.set(GRB_DoubleParam_BarConvTol, 0.0001);
		model.set(GRB_IntParam_Crossover, 0);
		model.set(GRB_IntParam_Threads,32);
		//model.set(GRB_IntParam_Presolve, 0);
		model.set(GRB_DoubleParam_OptimalityTol, 0.0001);

		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		double kkt1=0,kkt2=0;
		auto square = [](double x) {return x * x;  };
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {
				double y_value = y[colidx[i]].get(GRB_DoubleAttr_X);
				double m_value = m[colidx[i]].get(GRB_DoubleAttr_X);

				double up_value = model.getConstrByName(string("p") + to_string(i)).get(GRB_DoubleAttr_Pi);
				double un_value = model.getConstrByName(string("n") + to_string(i)).get(GRB_DoubleAttr_Pi);
				if (up_value -m_value + y_value - b[i] > 0) {
					kkt1 += square(m_value - y_value + b[i]);
				}
				else {
					kkt1 += square(up_value);
				}
				if (un_value - m_value - y_value + b[i] > 0) {
					kkt1 += square(m_value + y_value - b[i]);
				}
				else {
					kkt1 += square(un_value);
				}
				x[i] = y_value - b[i];
				double kt1 = m_value - std::abs(y_value - b[i]);
				if (kt1 < 0) {

				}
				else {
					kkt2 += kt1 * kt1;
				}
			}
		}
		double kkt1a = 0;
		for (auto index : viscons) {
			int tst = index % 4;
			index /= 4;
			int i, j, k;
			k = index % n;
			j = index / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = colidx[ij(i, j, n)];
			idx[1] = colidx[ij(i, k, n)];
			idx[2] = colidx[ij(j, k, n)];

			std::swap(idx[tst], idx[0]);
			double u_value = model.getConstrByName(to_string(index)).get(GRB_DoubleAttr_Pi);
			double y_value = y[idx[0]].get(GRB_DoubleAttr_X) - y[idx[1]].get(GRB_DoubleAttr_X) - y[idx[2]].get(GRB_DoubleAttr_X);
			if (y_value + u_value > 0) {
				kkt1a += y_value * y_value;

			}
			else {
				kkt1a += u_value * u_value;
			}
			//model.addConstr(y[idx[1]] - y[idx[0]] - y[idx[2]] <= 0);
			//model.addConstr(y[idx[2]] - y[idx[1]] - y[idx[0]] <= 0);
		}





		kkt_value[0] = kkt1 + kkt1a;

		cout << "u:" << get_unorm(model) << endl;
		//kkt_value[1]=get_unorm(model);
		kkt_value[1] = kkt2;
		
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}

	delete[] y;
	delete[] m;
	return model.get(GRB_DoubleAttr_ObjVal);
}

double gurobi::n2(std::set<int64_t>& viscons, std::vector<int>& colidx, int colidx_num)
{
	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar* y = model.addVars(colidx_num, GRB_CONTINUOUS);
	try
	{
		GRBQuadExpr obj = 0;


		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {

				y[colidx[i]].set(GRB_DoubleAttr_Start, x[i] + b[i]);

				obj += w[i] * (y[colidx[i]] - b[i]) * (y[colidx[i]] - b[i]);
			}

			//model.addGenConstrAbs(m[i.second], y[i.second] - b[i.first]);
		}


		add_constr_A(model, viscons, colidx, y);
		model.set(GRB_DoubleParam_BarQCPConvTol,1e-4);
		model.set(GRB_DoubleParam_OptimalityTol, 1e-4);
		model.set(GRB_DoubleParam_FeasibilityTol, 0.01);
		model.set(GRB_IntParam_Threads, 32);
		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;



		auto square = [](double x) {return x * x; };
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {
				double y_value = y[colidx[i]].get(GRB_DoubleAttr_X);

				x[i] = y_value - b[i];

			}
		}
		double kkt1a = 0;
		for (auto index : viscons) {
			int tst = index % 4;
			index /= 4;
			int i, j, k;
			k = index % n;
			j = index / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = colidx[ij(i, j, n)];
			idx[1] = colidx[ij(i, k, n)];
			idx[2] = colidx[ij(j, k, n)];

			std::swap(idx[tst], idx[0]);
			double u_value = model.getConstrByName(to_string(index)).get(GRB_DoubleAttr_Pi);
			double y_value = y[idx[0]].get(GRB_DoubleAttr_X) - y[idx[1]].get(GRB_DoubleAttr_X) - y[idx[2]].get(GRB_DoubleAttr_X);
			if (y_value + u_value > 0) {
				kkt1a += y_value * y_value;

			}
			else {
				kkt1a += u_value * u_value;
			}
			//model.addConstr(y[idx[1]] - y[idx[0]] - y[idx[2]] <= 0);
			//model.addConstr(y[idx[2]] - y[idx[1]] - y[idx[0]] <= 0);
		}

		kkt_value[0] = kkt1a;
		kkt_value[1] = 0;



	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}

	delete[] y;
	return model.get(GRB_DoubleAttr_ObjVal);
}

double gurobi::n2s(std::set<int64_t>& viscons, std::vector<int>& colidx)
{

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar* y = model.addVars(colidx.size(), GRB_CONTINUOUS);
	GRBVar* m = model.addVars(colidx.size(), GRB_CONTINUOUS);
	try
	{
		GRBQuadExpr obj = 0;


		for (int i = 0; i < var_num;i++) {
			y[colidx[i]].set(GRB_DoubleAttr_Start, x[i] + b[i]);
			m[colidx[i]].set(GRB_DoubleAttr_Start, fabs(x[i]));

			obj += m[colidx[i]] * w[i];
			obj += m[colidx[i]] * m[colidx[i]] * w[i]/20;
			obj+= (y[colidx[i]]-b[i] )* (y[colidx[i]] - b[i]) * w[i]/20;
			//model.addGenConstrAbs(m[i.second], y[i.second] - b[i.first]);
			model.addConstr(m[colidx[i]] - y[colidx[i]] + b[i] >= 0);
			model.addConstr(m[colidx[i]] + y[colidx[i]] - b[i] >= 0);
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
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {
				x[i] = y[colidx[i]].get(GRB_DoubleAttr_X) - b[i];
			}
		}
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}

	delete[] y;
	delete[] m;
	return model.get(GRB_DoubleAttr_ObjVal);
}
double gurobi::n_inf(std::set<int64_t>& viscons, std::vector<int>& colidx, int colidx_num)
{

	GRBEnv env = GRBEnv();
	GRBModel model = GRBModel(env);
	GRBVar* y = model.addVars(colidx_num, GRB_CONTINUOUS);
	GRBVar* m = model.addVars(1, GRB_CONTINUOUS);
	try
	{
		GRBLinExpr obj = 0;
		obj = m[0];
		m[0].set(GRB_DoubleAttr_Start, 0);
		for (int i = 0; i < var_num;i++) {
			if (colidx[i] != -1) {


				y[colidx[i]].set(GRB_DoubleAttr_Start, x[i] + b[i]);

				//model.addGenConstrAbs(m[i.second], y[i.second] - b[i.first]);
				model.addConstr( w[i] * (y[colidx[i]] - b[i])-m[0] <= 0,string("p")+to_string(i));
				model.addConstr(- w[i] * (y[colidx[i]] - b[i])-m[0] <= 0, string("n") + to_string(i));
			}
		}

		add_constr_A(model, viscons, colidx, y);


		model.set(GRB_IntParam_Method, 2);
		model.set(GRB_DoubleParam_FeasibilityTol, 0.01);
		model.set(GRB_DoubleParam_BarConvTol, 0.0001);
		model.set(GRB_IntParam_Crossover, 0);
		model.set(GRB_IntParam_Threads, 32);
		model.set(GRB_DoubleParam_OptimalityTol, 0.0001);

		model.setObjective(obj, GRB_MINIMIZE);
		model.optimize();

		cout << "Obj: " << model.get(GRB_DoubleAttr_ObjVal) << endl;
		cout << "m0:" << m[0].get(GRB_DoubleAttr_X) << endl;
		double maxres = 0;
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {
				double y_value = y[colidx[i]].get(GRB_DoubleAttr_X);

				maxres = std::max(maxres, std::abs(y_value - b[i]) * w[i]);

			}
		}
		cout << "trueobj:" << maxres << endl;
		double kkt1 = 0,kkt2=0;
		double m_value = m[0].get(GRB_DoubleAttr_X);
		auto square = [](double x) {return x * x; };
		for (int i = 0; i < var_num; i++) {
			if (colidx[i] != -1) {
				double y_value = y[colidx[i]].get(GRB_DoubleAttr_X);
				
				double up_value = model.getConstrByName(string("p") + to_string(i)).get(GRB_DoubleAttr_Pi);
				double un_value = model.getConstrByName(string("n") + to_string(i)).get(GRB_DoubleAttr_Pi);
				if (up_value-m_value+ (y_value - b[i]) * w[i] >0) {
					kkt1 += square(m_value - (y_value - b[i]) * w[i]);
				}
				else {
					kkt1 += square(up_value);
				}
				if (un_value-m_value - (y_value - b[i]) * w[i] > 0) {
					kkt1 += square(m_value + (y_value - b[i]) * w[i]);
				}
				else {
					kkt1 += square(un_value);
				}
				x[i] = y_value - b[i];
				double kt1 = m_value - std::abs(y_value - b[i])*w[i];
				if (kt1 < 0) {
					
				}
				else {
					kkt2 += kt1 * kt1;
				}
			}
		}
		double kkt1a = 0;
		for (auto index : viscons) {
			int tst = index % 4;
			index /= 4;
			int i, j, k;
			k = index % n;
			j = index / n;
			i = j / n;
			j = j % n;
			int idx[3];
			idx[0] = colidx[ij(i, j, n)];
			idx[1] = colidx[ij(i, k, n)];
			idx[2] = colidx[ij(j, k, n)];

			std::swap(idx[tst], idx[0]);
			double u_value = model.getConstrByName(to_string(index)).get(GRB_DoubleAttr_Pi);
			double y_value = y[idx[0]].get(GRB_DoubleAttr_X) - y[idx[1]].get(GRB_DoubleAttr_X) - y[idx[2]].get(GRB_DoubleAttr_X);
			if (y_value +u_value> 0) {
				kkt1a += y_value * y_value;
				
			}
			else {
				kkt1a += u_value * u_value;
			}
			//model.addConstr(y[idx[1]] - y[idx[0]] - y[idx[2]] <= 0);
			//model.addConstr(y[idx[2]] - y[idx[1]] - y[idx[0]] <= 0);
		}





		kkt_value[0] = kkt1+kkt1a;





		cout << "u:" << get_unorm(model) << endl;
		kkt_value[1] = kkt2;
	}
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;

	}

	delete[] y;
	delete[] m;
	return model.get(GRB_DoubleAttr_ObjVal);
}

void gurobi::add_constr_A(GRBModel& model, std::set<int64_t>& viscons, std::vector<int>& colidx, GRBVar*y)
{
	for (auto index : viscons) {
		int tst = index % 4;
		index /= 4;
		int i, j, k;
		k = index % n;
		j = index / n;
		i = j / n;
		j = j % n;
		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];

		std::swap(idx[tst], idx[0]);
		model.addConstr(y[idx[0]] - y[idx[1]] - y[idx[2]] <= 0,to_string(index));
		//model.addConstr(y[idx[1]] - y[idx[0]] - y[idx[2]] <= 0);
		//model.addConstr(y[idx[2]] - y[idx[1]] - y[idx[0]] <= 0);
	}
}

double gurobi::get_unorm(GRBModel& model)
{
	double kkt = 0;
	auto conspointer = model.getConstrs();
	int consnum = model.get(GRB_IntAttr_NumConstrs);
#pragma omp parallel for reduction(+:kkt)
	for (int i = 0; i < consnum; i++) {
		double u_value = conspointer[i].get(GRB_DoubleAttr_Pi);
		kkt += u_value * u_value;
	}
	delete[] conspointer;
	return kkt;
}
