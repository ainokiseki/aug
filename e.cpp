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
	norm_ori s2("graph8k.txt");
	s2.run();
	
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
