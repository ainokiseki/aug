#pragma once
#include "metric.h"
#include<thread>
#include<mutex>
#include"tri.h"
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
class jk :public metric {

public:

	jk(string graph_name, int type = 0, double g = 10);
	~jk();
	void multi_triloop(int threadnum, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction);

	void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, atomic_int& home, mutex& mut, vector<queue<pid>>& vc, vector<queue<pid>>& vnc);

	//void worker(std::queue<looplist>& q, condition_variable& cond, atomic_int& ant, atomic_int& home, mutex& mut);
	//void generator(int threadnum, queue<looplist> q, vector<queue<pid>>& correction, vector<queue<pid>>& next_correction, std::mutex& mut, std::condition_variable& cond);
	void triloop(int threadnum, std::queue<pid>& correction, std::queue<pid>& next_correction, mutex* mut);

	void triloop2(int istart, int iend, int kstart, int blocksize, std::queue<pid>& correction, std::queue<pid>& next_correction);


	void biloop();
	void multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction);

	void multi_triloop(std::vector<std::queue<pid> >& correction, std::vector<std::queue<pid> >& next_correction, int te);






	double run(double tol);
	void init();
	double gam;
	double kkt_2;
	double kkt_1;
	Eigen::VectorXd y, m, inv_w, p;
	int mthread;

};





class dk : public metric
{
public:

	dk(string graph_name, int type = 0, double g = 10);
	~dk();
	void triloop(std::unordered_map<int64_t, double>* correction, std::unordered_map<int64_t, double>* next_correction);
	void biloop();
	void multi_triloop(std::queue<pid>* correction, std::queue<pid>* next_correction, int te, int wtf);
	void multi_triloop(std::unordered_map<int64_t, double>* correction, std::unordered_map<int64_t, double>* next_correction, int te);
	void multi_triloop(std::unordered_map<int64_t, double>* correction, std::unordered_map<int64_t, double>* next_correction, std::vector<int64_t>* cons);
	void multi_triloop(std::unordered_map<int64_t, double>* correction, std::unordered_map<int64_t, double>* next_correction);






	int run_cg(double threshod);

	double run_test(double threshod);





	void run();
	void init();
	double gam;
	double kkt_2;
	double kkt_1;
	Eigen::MatrixXd y, m, inv_w, mat_b, p;

};

class dk_triloop
{
public:

	dk_triloop(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::unordered_map<int64_t, double>* c1, std::unordered_map<int64_t, double>* c2);
	~dk_triloop();
	void operator()(int istart, int iend, int kstart, int blocksize);


	int n;
	Eigen::MatrixXd* y;
	Eigen::MatrixXd* m;
	Eigen::MatrixXd* inv_w;
	Eigen::MatrixXd* mat_b;
	std::unordered_map<int64_t, double>* correction;
	std::unordered_map<int64_t, double>* next_correction;
	std::mutex mut;
};

class dk_mutex :public dk_triloop
{
public:
	dk_mutex(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::unordered_map<int64_t, double>* c1, std::unordered_map<int64_t, double>* c2);
	~dk_mutex();
	virtual void operator()(int istart, int iend);
	std::mutex* mute;


};
class dk_mutex_test
{
public:

	dk_mutex_test(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::queue<pid>* c1, std::queue<pid>* c2);

	~dk_mutex_test();
	virtual void operator()(int istart, int iend);
	std::mutex* mute;
	int n;
	Eigen::MatrixXd* y;
	Eigen::MatrixXd* m;
	Eigen::MatrixXd* inv_w;
	Eigen::MatrixXd* mat_b;
	std::queue<pid>* correction;
	std::queue<pid>* next_correction;
	std::mutex mut;

};
class dk_mutex_cg :public dk_mutex
{
public:
	dk_mutex_cg(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::unordered_map<int64_t, double>* c1, std::unordered_map<int64_t, double>* c2, std::vector<int64_t>* cons);
	~dk_mutex_cg();
	void operator()(int istart, int iend);
	std::vector<int64_t>* cons;


};

