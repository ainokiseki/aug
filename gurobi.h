#pragma once
#include"metric.h"
#include"gurobi_c++.h"
class gurobi:public metric
{
public:
	gurobi(string graphname);
	~gurobi();
	double run(int type);
	double n1(std::set<int64_t>&viscons , std::unordered_map<int, int>&colidx);
	double n2(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx);
private:

};

