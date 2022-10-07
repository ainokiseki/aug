#pragma once
#include "metric.h"
#include<thread>
#include<mutex>
#include"tri.h"
#include<concurrent_vector.h>
#include<queue>
typedef std::pair<int64_t, double> pid;
class looplist
{
public:
	int istart;
	int iend;
	int kstart;
	int kend;
	int cnum;

	looplist(int istart, int iend, int kstart, int kend, int cnum) :
		istart(istart), iend(iend), kstart(kstart), kend(kend), cnum(cnum)
	{


	}
	~looplist() {}
};

class project_method :public metric {
public:
	project_method(string graph_name, int type = 0, double g = 1);

	~project_method();
	void multi_triloop2(int threadnum, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction);
	void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, mutex& mut, vector<queue<pid>>& vc, vector<queue<pid>>& vnc,int* tnum);
	double gam;
	double kkt_2;
	double kkt_1;
	Eigen::VectorXd y, m, inv_w, p;
	bool mconstrain_check(double tol = 1e-2);
	void constrain_check(double tol, atomic_int& ant, int istart, int iend);
};

class jk :public project_method {

public:

	jk(string graph_name, int type = 0, double g = 1);
	~jk();
	//void multi_triloop2(int threadnum, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction);
	void multi_triloop(int threadnum, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction);

	

	//void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, atomic_int& home, mutex& mut, vector<queue<pid>>& vc, vector<queue<pid>>& vnc);

	//void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, atomic_int& home, mutex& mut);
	//void generator(int threadnum, queue<looplist> q, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction, std::mutex& mut, std::condition_variable& cond);
	void triloop(int threadnum, std::queue<pid>& correction, std::queue<pid>& next_correction, mutex* mut);

	//void triloop2(int istart, int iend, int kstart, int blocksize, std::queue<pid>& correction, std::queue<pid>& next_correction);


	void biloop();
	void multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction);

	void multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction, int te);




	//void constrain_check();

	double run(double gaptol,double contol,int statusfreq,int maxits);
	void init();



};

class jk_square:public project_method
{
public:

	jk_square(string graph_name, int type=0);
	~jk_square();

	double run(double gaptol, double contol, int statusfreq, int maxits);
	const double gam = 0.5;
private:

};

class jk_inf :public project_method {

public:

	jk_inf(string graph_name, int type = 0, double g =500);
	~jk_inf();
	//void multi_triloop2(int threadnum, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction);
	//void multi_triloop(int threadnum, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction);



	//void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, atomic_int& home, mutex& mut, vector<queue<pid>>& vc, vector<queue<pid>>& vnc);

	//void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, atomic_int& home, mutex& mut);
	//void generator(int threadnum, queue<looplist> q, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction, std::mutex& mut, std::condition_variable& cond);
	//void triloop(int threadnum, std::queue<pid>& correction, std::queue<pid>& next_correction, mutex* mut);

	//void triloop2(int istart, int iend, int kstart, int blocksize, std::queue<pid>& correction, std::queue<pid>& next_correction);


	void biloop();
	//void multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction);

	//void multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction, int te);




	//void constrain_check();

	double run(double gaptol, double contol, int statusfreq, int maxits);
	//void init();



};