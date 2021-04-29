#include"norm_class.h"
#include"dk.h"
#include<chrono>
//#include"admm_class.h"


int main() {
	auto i_start = std::chrono::high_resolution_clock::now();
	dk nm("graph150.txt");
	nm.run();
	auto i_end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	double tm = double(eld.count()) * 1.0 * 1e-9;
	cout << "time is:" << tm << endl;
	//norm_ori s("graph150.txt");
	//s.run();
    return 0;
}
