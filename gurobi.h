#pragma once
#include"metric.h"
#include"gurobi_c++.h"
class gurobi:public metric
{
public:
	gurobi(string graphname,int type=0);
	~gurobi();
	double run(int type);
	double n1(std::set<int64_t>&viscons , std::unordered_map<int, int>&colidx);
	double n2(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx);
	double n2s(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx);
	double n_inf(std::set<int64_t>& viscons, std::unordered_map<int, int>& colidx);
private:

};

