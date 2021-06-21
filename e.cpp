#include"norm_class.h"
#include"dk.h"
#include"admm_class.h"
#include<chrono>
#include<string>
//#include"admm_class.h"


int main() {

	//vector<double>vx = { 1,5,20};
	//vector<double>vt = { 1e-5};
	//vector<string>vf = { "100","150","200","300","500" };
	//for (auto i : vx) {

	//	for (auto j : vt) {
	//		for (auto k : vf) {

	//			dk nm("dfsample"+k+".txt", i);
	//			auto k = nm.run_test(j);
	//			cout << "gam:" << i << "\tthreshold:" << j << "\ttimes:" << k.second << endl;
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


	//string s;
	//cin >> s;
	//if (s.substr(s.length() - 3, s.length()) == "dat") {
	//	string ssubs = s.substr(s.find('/') + 1);
	//	string subs=ssubs.substr(ssubs.find('_') + 1);
	//	subs = subs.substr(0,subs.find('_'));

	//	cout << ssubs << "\t" << subs << "\t" << endl;
	//	norm_ori s2(subs, s);
	//	s2.run();
	//}
	//else {
	//	norm_ori s2(s);
	//	s2.run();
	//}
	norm_ori a("graph500.txt");
	int n = a.n;
	int var_num = a.var_num;
	Eigen::MatrixXd w(n, n),b(n,n);

	int ct = 0;
	for (int i = 0; i < n; i++) {
		for (int j = i+1; j < n; j++,ct++) {
			w(i, j) = w(j, i) = a.w[ct];
			b(i, j) = b(j, i) = a.b[ct];
		}
	}
	std::vector<int>t1, t2;
	for (int i = 0; i < 250; i++)
		t1.push_back(i);
	for (int i = 250; i < n; i++)
		t2.push_back(i);
	std::sort(t1.begin(), t1.end());
	std::sort(t2.begin(), t2.end());
	//{
	//	Eigen::VectorXd w1(t1.size() * (t1.size() - 1) / 2), b1(t1.size() * (t1.size() - 1) / 2);
	//	ct = 0;
	//	for (int i = 0; i < t1.size(); i++) {
	//		for (int j = i + 1; j < t1.size(); j++, ct++) {
	//			w1[ct] = w(t1[i], t1[j]);
	//			b1[ct] = b(t1[i], t1[j]);
	//		}
	//	}
	//	norm_ori t(t1.size(), b1, w1);
	//	t.run();
	//	ct = 0;
	//	for (int i = 0; i < t1.size(); i++) {
	//		for (int j = i + 1; j < t1.size(); j++, ct++) {
	//			a.x[ij(t1[i], t1[j], n)] = t.x[ct];
	//			a.v[ij(t1[i], t1[j], n)] = t.v[ct];
	//		}
	//	}
	//	int n1 = t1.size();
	//	for (auto i : t.u) {
	//		int ui, uj, uk;
	//		uk = i.first % n1;
	//		uj = i.first / n1;
	//		ui = uj / n1;
	//		uj = uj % n1;
	//		int64_t res = static_cast<int64_t>(ui) * n * n + static_cast<int64_t>(uj) * n + uk;
	//		a.u.insert(std::pair<int64_t, double>(res, i.second));
	//	}
	//}
	//{
	//	Eigen::VectorXd w1(t2.size() * (t2.size() - 1) / 2), b1(t2.size() * (t2.size() - 1) / 2);
	//	ct = 0;
	//	for (int i = 0; i < t2.size(); i++) {
	//		for (int j = i + 1; j < t2.size(); j++, ct++) {
	//			w1[ct] = w(t2[i], t2[j]);
	//			b1[ct] = b(t2[i], t2[j]);
	//		}
	//	}
	//	norm_ori t(t2.size(), b1, w1);
	//	t.run();
	//	ct = 0;
	//	for (int i = 0; i < t2.size(); i++) {
	//		for (int j = i + 1; j < t2.size(); j++, ct++) {
	//			a.x[ij(t2[i], t2[j], n)] = t.x[ct];
	//			a.v[ij(t2[i], t2[j], n)] = t.v[ct];
	//		}
	//	}
	//	int n1 = t2.size();
	//	for (auto i : t.u) {
	//		int ui, uj, uk;
	//		uk = i.first % n1;
	//		uj = i.first / n1;
	//		ui = uj / n1;
	//		uj = uj % n1;
	//		int64_t res = static_cast<int64_t>(ui) * n * n + static_cast<int64_t>(uj) * n + uk;
	//		a.u.insert(std::pair<int64_t, double>(res, i.second));
	//	}
	//}
	a.run();
	//string s;
	//ofstream series("series.txt", std::ios::out);
	//for (int i = 100; i <= 1000; i += 50) {
	//	s = string("dfsample") + to_string(i) + string(".txt");
	//	{
	//		norm_ori s2(s);
	//		auto tm = s2.run();
	//		series << s << ":" << tm.first <<"\t"<<tm.second.first<<"\t"<<tm.second.second<< endl;
	//	}
	//	{
	//		dk s3(s);
	//		auto k=s3.run_test(1e-5);
	//		series << s << ":" << k.first << "\t" << k.second << endl;
	//	}
	//}


	//string s;
	//ofstream series("series.txt", std::ios::app);
	//for (int i = 100; i <= 1000; i += 50) {
	//	s = string("dfsample") + to_string(i) + string(".txt");
	//	{
	//		double gama[3] = { 1,5,10 };
	//		for (int j = 0; j < 3; j++) {


	//			dk s3(s,gama[j]);

	//			auto k = s3.run_test(1e-6);
	//			series << s << ":" <<gama[j]<<"\t"<< k.first << "\t" << k.second << endl;
	//		}
	//	}
	//}



	//dk s2("dfsample100.txt");
	//s2.run_test(1e-5);
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
