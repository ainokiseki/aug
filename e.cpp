#include"norm_class.h"
#include"dk.h"
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
	// 
	// 
	string s;
	cin >> s;
	if (s.substr(s.length() - 3, s.length()) == "dat") {
		string ssubs = s.substr(s.find('/') + 1);
		string subs=ssubs.substr(ssubs.find('_') + 1);
		subs = subs.substr(0,subs.find('_'));

		cout << ssubs << "\t" << subs << "\t" << endl;
		norm_ori s2(subs, s);
		s2.run();
	}
	else {
		norm_ori s2(s);
		s2.run();
	}
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
