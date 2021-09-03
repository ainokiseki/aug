#include"norm_class.h"
#include"dk.h"
#include<chrono>
#include<string>
#include"gurobi.h"
#include"admm_class.h"
#include"mytime.h"
//#include"admm_class.h"
typedef std::pair<double, double> pdd;
typedef std::pair<pdd, pdd>ppdd;
std::pair<double,double> test(Eigen::VectorXd mat,int n) {

	double dist = 0;
	double maxdist = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {

				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);

				double res2;
				res2 = mat[idx[0]] - mat[idx[1]] - mat[idx[2]];
				dist += fmax(0, res2);
				maxdist = fmax(maxdist, res2);
				res2 = mat[idx[1]] - mat[idx[0]] - mat[idx[2]];
				dist += fmax(0, res2);
				maxdist = fmax(maxdist, res2);
				res2 = mat[idx[2]] - mat[idx[1]] - mat[idx[0]];
				dist += fmax(0, res2);
				maxdist = fmax(maxdist, res2);
			}
		}
	}
	cout << "sum:" << dist << endl;
	cout << "max:" << maxdist << endl;
	return std::pair<double,double>(dist, maxdist);
}
int main() {

	//vector<double>vx = { 1,5,20};
	//vector<double>vt = { 1e-5};
	//vector<string>vf = { "100","150","200","300","500" };
	//for (auto i : vx) {

	//	for (auto j : vt) {
	//		for (auto k : vf) {

	//			dk nm("graph"+k+".txt", i);
	//			int k = nm.run_test(j);
	//			cout << "gam:" << i << "\tthreshold:" << j << "\ttimes:" << k << endl;
	//		}

	//	}
	//}

	//
	//auto i_start = std::chrono::high_resolution_clock::now();
	//dk nm("graph1.txt");
	//nm.run_test(1e-5);
	//auto i_end = std::chrono::high_resolution_clock::now();
	//auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	//double tm = double(eld.count()) * 1.0 * 1e-9;
	//cout << "time is:" << tm << endl;
	//cin >> tm;
	//dk s("graph500.txt",5);
	//s.run_test(1e-2);
	// 
	// 
	vector<string>target = {"graph1k.txt"};
	vector<ppdd>res;
	timer tm;
	for (auto i : target) {
		dk s1(i,1);
		//admm_basic s2(i);

		double a1, a2, a3;
		for (int j = 0; j < s1.n * (s1.n - 1) / 2; j++) {
			s1.w[j]  = 1;
		}

		cout << s1.w.size() << endl;
		tm.tic();
		a1=s1.run_test(1e-3);
		tm.toc();
		auto c1 = pdd(a1, tm.blink);
		auto b1=ppdd(c1,test(s1.x+s1.b,s1.n));



		res.push_back(b1);		//res.push_back(b2);		

	}
	//vector<double>tau = { 1.618,1.9 };
	//vector<double>sigma = { 1,0.1,0.001 };
	//vector<double>scale_c = { 1,9 };
	//for (auto i : tau) {
	//	for (auto j : sigma) {
	//		for (auto k : scale_c) {

	//			admm_basic s2("graph500.txt");
	//			tm.tic();
	//			double rs=s2.run(i, j, 1, k);
	//			tm.toc();
	//			res.push_back(ppdd(pdd(rs, tm.blink), test(s2.x + s2.b, s2.n)));
	//		}
	//	}
	//}
	//admm_basic s2("graph500.txt");

	//	for (int j = 0; j < s2.n * (s2.n - 1) / 2; j++) {
	//	s2.w[j] = 1;
	//}

	//tm.tic();
	//double tp=s2.run();
	//tm.toc();
	//res.push_back(ppdd(pdd(tp, tm.blink), test(s2.x + s2.b, s2.n)));
	for (auto i : res) {
		cout << "res:" << i.first.first<<"\t"<<i.first.second<< "\t" << "sum:" << i.second.first << " " << "max:" << i.second.second << endl;
	}
	//norm_ori s2(s);
	//s2.run_dense();
	//s2.run(0.001, 0, 20);

	//if (s.substr(s.length() - 3, s.length()) == "dat") {
	//	string ssubs = s.substr(s.find('/') + 1);
	//	string subs=ssubs.substr(ssubs.find('_') + 1);
	//	subs = subs.substr(0,subs.find('_'));

	//	cout << ssubs << "\t" << subs << "\t" << endl;
	//	cout << "cst=delta/sigma,plz input delta(do not input 1e-3, use 0.001 instead):";
	//	double cst; int cons;
	//	cin >> cst;
	//	cout << "number of constraints is no less than cons*n*(n-1)/2,plz input cons:";
	//	cin >> cons;

	//	norm_ori s2(subs, s);
	//	s2.run(cst,cons,20);
	//}
	//else {
	// 
		//cout << "cst=delta/sigma,plz input delta(do not input 1e-3, use 0.001 instead):";
		//double cst; int cons; int sub_iter;
		//cin >> cst;
		//cout << "number of constraints is no less than cons*n*(n-1)/2,plz input cons:";
		//cin >> cons;
		//cout << "plz input maxium steps of subiteration:";
		//cin >> sub_iter;
		//norm_ori s2(s);
		//s2.run(cst, cons,sub_iter);


	//}
	//cin >> s;

	//cin >> s;
	// 
	// 
	// 
	//dk s2(s);
	//s2.run_test(1e-3);
	// 
	// 
	// 
	// 
	// 
	// 
	//string dataname;
	//cin >> dataname;
	//norm_2 s(dataname+string(".txt"));
	//s.run();
	//string tem;
	//cin >> tem;
    return 0;
}
