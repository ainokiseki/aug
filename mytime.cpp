#include "mytime.h"

timer::timer()
{
}

timer::~timer()
{
}

void timer::tic()
{
	start = std::chrono::high_resolution_clock::now();
}

void timer::toc()
{
	end = std::chrono::high_resolution_clock::now();
	auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	tm += double(eld.count()) * 1.0 * 1e-9;

}

double timer::time()
{
	return tm;
}

