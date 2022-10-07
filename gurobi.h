#pragma once
#include"metric.h"
#include"cons_generate.h"
#include"gurobi_c++.h"
class gurobi:public metric
{
public:
	gurobi(string graphname,int type=0);
	~gurobi();
	double run(int type,int maxstep=9999);
	double n1(std::set<int64_t>&viscons , std::vector<int>& colidx,int colidx_num);
	double n2(std::set<int64_t>& viscons, std::vector<int>& colidx, int colidx_num);
	double n2s(std::set<int64_t>& viscons, std::vector<int>& colidx);
	double n_inf(std::set<int64_t>& viscons, std::vector<int>& colidx, int colidx_num);
	double kkt_value[2];
	void add_constr_A(GRBModel&model, std::set<int64_t>& viscons, std::vector<int>& colidx, GRBVar* y);
	double get_unorm(GRBModel& model);
private:

};

class cons_generator_pro_gurobi
{
public:
	timer ava;
	std::vector<int>colidx;
	gurobi* gurobi_obj;
	Eigen::VectorXd diag;
	int n;
	int var_num;
	int colidx_num;
	int initconsnum = 0;
	int mthread;
	bool initflag = true;

	float limit = 1e-2;
	std::set<int64_t>viscons;
	cons_generator_pro_gurobi(gurobi* gurobi_obj);
	~cons_generator_pro_gurobi() = default;


	inline std::array<int, 3> unzip(int64_t i_, int& i, int& j, int& k, int& tst) {
		tst = i_ % 4;
		i_ /= 4;
		k = i_ % n;
		j = i_ / n;
		i = j / n;
		j = j % n;
		std::array<int, 3> res{ ij(i, j, n),ij(i,k,n),ij(j,k,n) };
		std::swap(res[0], res[tst]);
		return res;
	}
	inline float product_with_a(const Eigen::VectorXd& x, std::array<int, 3>& idx) {
		return x[idx[0]] - x[idx[1]] - x[idx[2]];
	}
	void generate_cheka();
	void generate_init() {

		std::fill(colidx.begin(), colidx.end(), -1);
		colidx_num = 0;
		//viscons.clear();
		//data.clear();
	}
	void generate_loop(double& ttk, double& ttk2, double& total_ax);
	void generate(double& ttk, double& ttk2, double& total_ax) {
		ava.tic();
#ifndef NOEXIT
		generate_init();

		generate_cheka();
#endif
		if (initflag) {
			gurobi_generate_set_para();
#ifdef NOEXIT
			initconsnum = viscons.size();
#else
			initconsnum = var_num / 2;
#endif
			initflag = false;
		}
		else {

			generate_loop(ttk, ttk2, total_ax);
#ifdef NOEXIT
			initconsnum *= 2;
#endif
		}
		ava.toc();
		std::cout << "cg time:" << ava.time() << std::endl;
	}
	void gurobi_generate_set_para();

};