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
	double blink = 0;
private:
	std::chrono::steady_clock::time_point start,end;
	double tm=0;

};

