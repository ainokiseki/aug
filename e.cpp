#include"norm_class.h"
#include"dk.h"
#include<chrono>
//#include"admm_class.h"


int main() {

	//vector<double>vx = { 1,2,5,10,20,40,60,80,100,150,200 };
	//vector<double>vt = { 1e-2};
	//for (auto i : vx) {

	//	for (auto j : vt) {
	//		dk nm("graph300.txt", i);
	//		int k=nm.run_test(j);
	//		cout << "gam:" << i << "\tthreshold:" << j <<"\ttimes:"<<k<< endl;
	//	}
	//}

	/*
	auto i_start = std::chrono::high_resolution_clock::now();
	dk nm("graph200.txt");
	nm.run();
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	cout << "time is:" << tm << endl;
	*/
	//dk s("graph500.txt",5);
	//s.run_test(1e-2);
	// 
	//dk s2("graph500.txt", 5);
	//s2.run_test(1e-3);
	// 
	norm_ori s("graph1.txt");
	s.run();
    return 0;
}
