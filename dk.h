#pragma once
#include "metric.h"
#include<thread>
#include<mutex>
#include"tri.h"
typedef std::pair<int64_t, double> pid;

class jk :public metric {

public:

	jk(string graph_name, int type = 0, double g = 10);
	~jk();
	void triloop(std::unordered_map<int64_t, double>* correction, std::unordered_map<int64_t, double>* next_correction);
	void biloop();
	void multi_triloop(std::vector<pid>* correction, std::vector<pid>* next_correction);

	int run_cg(double threshod);

	double run_test(double threshod);

	void run();
	void init();
	double gam;
	double kkt_2;
	double kkt_1;
	Eigen::MatrixXd y, m, inv_w, mat_b, p;




};





class dk : public metric
{
public:

	dk(string graph_name, int type=0,double g=10);
	~dk();
	void triloop(std::unordered_map<int64_t, double>* correction, std::unordered_map<int64_t, double>* next_correction);
	void biloop();
	void multi_triloop(std::vector<pid>* correction, std::vector<pid>* next_correction, int te, int wtf);
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
	Eigen::MatrixXd y, m,inv_w,mat_b,p;

};

class dk_triloop
{
public:

	dk_triloop(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::unordered_map<int64_t, double>* c1, std::unordered_map<int64_t, double>* c2);
	~dk_triloop();
	void operator()(int istart,int iend, int kstart, int blocksize);
	

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

	dk_mutex_test(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::vector<pid>* c1, std::vector<pid>* c2);

	~dk_mutex_test();
	virtual void operator()(int istart, int iend);
	std::mutex* mute;
	int n;
	Eigen::MatrixXd* y;
	Eigen::MatrixXd* m;
	Eigen::MatrixXd* inv_w;
	Eigen::MatrixXd* mat_b;
	std::vector<pid>* correction;
	std::vector<pid>* next_correction;
	std::mutex mut;

};
class dk_mutex_cg :public dk_mutex
{
public:
	dk_mutex_cg(int n, Eigen::MatrixXd* y, Eigen::MatrixXd* m, Eigen::MatrixXd* inv_w, Eigen::MatrixXd* mat_b, std::unordered_map<int64_t, double>* c1, std::unordered_map<int64_t, double>* c2,std::vector<int64_t>*cons);
	~dk_mutex_cg();
	void operator()(int istart, int iend);
	std::vector<int64_t>*cons;


};

