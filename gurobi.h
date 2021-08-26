#pragma once
#include"metric.h"
#include"gurobi_c++.h"
class gurobi:public metric
{
public:
	gurobi(string graphname);
	~gurobi();
	void run();
	void n1(std::set<int64_t>&viscons , std::unordered_map<int, int>&colidx);
	void n2(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx);
private:

};

