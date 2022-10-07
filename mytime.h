#pragma once
#include <chrono>

class timer
{
public:
	timer()=default;
	~timer()=default;
	void tic() {
		start = std::chrono::high_resolution_clock::now();
	}
	void toc() {
		end = std::chrono::high_resolution_clock::now();
		auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
		tm += double(eld.count()) * 1.0 * 1e-9;
		blink = double(eld.count()) * 1.0 * 1e-9;
	
	}
	double time() {
		return tm;
	}
	double blink = 0;
private:
	std::chrono::high_resolution_clock::time_point start,end;
	double tm=0;

};



