#pragma once
#include <chrono>

class timer
{
public:
	timer();
	~timer();
	void tic();
	void toc();
	double time();
private:
	std::chrono::steady_clock::time_point start,end;
	double tm=0;
};

