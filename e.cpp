#include"norm_class.h"

#include<chrono>
#include<string>
#include"dk.h"

#include"gurobi.h"
#include"clicpp.h"
#include"tri.h"
#include"mytime.h"
//#define SINGLETEST
//#define ADMMTEST
//#define ADMMTEST
#define MULTITEST
//#define REL
//#include"admm_class.h"
typedef std::pair<double, double> pdd;
typedef std::pair<pdd, pdd>ppdd;
template<typename T1,typename T2>
std::pair<double,double> test( T1 &x, T2 &b,int n) {
	//return std::pair<double, double>(0, 0);
	double dist = 0;
	double maxdist = 0;
#ifndef REL
	auto check = [&x,&b,n](double&distt,double&maxdistt,int istart,int iend) {
		double d1=0;
		double m1=0;
		for (int i = istart; i < iend; i++) {
			for (int j = i + 1; j < n; j++) {
				for (int k = j + 1; k < n; k++) {

					int idx[3];
					idx[0] = ij(i, j, n);
					idx[1] = ij(i, k, n);
					idx[2] = ij(j, k, n);

					double res2;
					res2 = x[idx[0]] - x[idx[1]] - x[idx[2]]+  b[idx[0]] - b[idx[1]] - b[idx[2]];
					d1 += std::max(0.0, res2);
					m1 = max(m1, res2);
					res2 = x[idx[1]] - x[idx[0]] - x[idx[2]] + b[idx[1]] - b[idx[0]] - b[idx[2]];
					d1 += std::max(0.0, res2);
					m1 = max(m1, res2);
					res2 = x[idx[2]] - x[idx[1]] - x[idx[0]] + b[idx[2]] - b[idx[1]] - b[idx[0]];
					d1 += std::max(0.0, res2);
					m1 = max(m1, res2);
				}
			}
		}
		distt = d1;
		maxdistt = m1;
		//cout << istart << " " << iend << endl;
	};
	vector<thread>vt;
	vector<double>vdist(32, 0), vmax(32, 0);
	for (int i = 0; i < 32; i++) {
		//cout << n - (int)(cbrt((i + 1) * 1.0 / 32) * n) << " " << n - (int)(cbrt(i * 1.0 / 32) * 1.0 * n) << endl;
		vt.push_back(thread(check,std::ref(vdist[i]),std::ref(vmax[i]), n - (int)(cbrt((i + 1) * 1.0 / 32) * n), n - (int)(cbrt(i * 1.0 / 32) * 1.0 * n)));
	}
	for (int i = 0; i < 32; i++) {
		vt[i].join();
		dist+= vdist[i];
		maxdist = fmax(maxdist, vmax[i]);
	}

	cout << "sum:" << dist << endl;
	cout << "max:" << maxdist << endl;
	return std::pair<double,double>(dist, maxdist);
#endif
	return std::pair<double, double>(0, 0);
}


Eigen::VectorXd admm_inf_pi_c(Eigen::VectorXd x, double sigma) {
	Eigen::VectorXd tmp = x.cwiseAbs();
	double gap = tmp.sum() - 1 / sigma;

	if (gap <= 0)
		return x;

	sort(tmp.data(), tmp.data() + tmp.size());
	double act;
	int len = tmp.size();
	for (int i = 0; i < len; i++) {

		if ((len - i) * tmp[i] >= gap) {
			act = gap / (len - i);
			break;
		}
		else {
			gap -= tmp[i];
		}
	}
	return x.cwiseSign().cwiseProduct(x - act * x.cwiseSign()).cwiseMax(0.0).cwiseProduct(x.cwiseSign());


}

int main(int argc,char**argv) {

	//vector<double>vx = { 1,5,20};
	//vector<double>vt = { 1e-5};
	//vector<string>vf = { "100","150","200","300","500" };
	//for (auto i : vx) {

	//	for (auto j : vt) {
	//		for (auto k : vf) {

	//			dk nm("graph"+k+".txt", i);
	//			int k = nm.run_test(j);
	//			cout << "gam:" << i << "\tthreshold:" << j << "\ttimes:" << k << endl;
	//		}

	//	}
	//}

	//
	//auto i_start = std::chrono::high_resolution_clock::now();
	//dk nm("graph1.txt");
	//nm.run_test(1e-5);
	//auto i_end = std::chrono::high_resolution_clock::now();
	//auto eld = std::chrono::duration_cast<std::chrono::nanoseconds>(i_end - i_start);
	//double tm = double(eld.count()) * 1.0 * 1e-9;
	//cout << "time is:" << tm << endl;
	//cin >> tm;
	//dk s("graph500.txt",5);
	//s.run_test(1e-2);
	// 
	// 
	//vector<string>target = { "celegansneural.txt","celegans_metabolic.txt","email.txt","Erdos991.txt","Harvard500.txt","jazz.txt","netscience.txt","Roget.txt","SmaGri.txt","SmallW.txt","USAir97.txt"};
	//vector<string>target = { "power.txt" };
	//vector<string>target = { "ca-GrQc.txt","ca-HepTh.txt","polblogs.txt","power.txt" };

	vector<int>mxstep;
	//int gurobi_type = stoi(argv[1]);


	 
	//vector<string>target = { "graph500.txt" };
	vector<ppdd>res;
	timer tm;

	//norm_ori s2(string("./flit/celegansneural.txt"));
	//s2.run_small(1e-3,0,20);
	// 	   =====================================almtest=============================
	//for (auto i : target) {
	//	norm_ori s1(string("./flit/")+i+string(".txt"));   
	//	norm_ori s2(string("./flit/")+i+string(".txt"));

	//	double a1, a2, a3;
	//	double tst = 1;

	//	//for (int j = 0; j < s1.n * (s1.n - 1) / 2; j++) {
	//	//	s1.w[j]  = std::fmax(s1.w[j],1);
	//	//}


	//	cout << "tst:" << tst << endl;
	//	cout << s1.w.size() << endl;
	//	tm.tic();
	//	a1=s1.run(0.001,0,30);
	//	tm.toc();
	//	auto c1 = pdd(a1, tm.blink);
	//	auto b1=ppdd(c1,test(s1.x+s1.b,s1.n));
	//	res.push_back(b1);

	//	cout << "tst:" << tst << endl;
	//	cout << s2.w.size() << endl;
	//	tm.tic();
	//	a2 = s2.run_small(0.001, 0, 20);
	//	tm.toc();
	//	auto c2 = pdd(a2, tm.blink);
	//	auto b2 = ppdd(c2, test(s2.x + s2.b, s2.n));

	//	res.push_back(b2);		//res.push_back(b2);		

	//}
	//  =====================================almtest=============================
	////metric m(string("graph14.txt"));
	////for (int i = 0; i < 7 * 13; i++) {
	////	cout << m.w[i] << " " << m.b[i] << endl;
	////}
	////vector<double>tau = { 1.618,1.9 };
	////vector<double>sigma = { 1,0.1,0.001 };
	////vector<double>scale_c = { 1,9 };
	////for (auto i : tau) {
	////	for (auto j : sigma) {
	////		for (auto k : scale_c) {

	////			admm_basic s2("graph500.txt");
	////			tm.tic();
	////			double rs=s2.run(i, j, 1, k);
	////			tm.toc();
	////			res.push_back(ppdd(pdd(rs, tm.blink), test(s2.x + s2.b, s2.n)));
	////		}
	////	}
	////}
#ifdef SINGLETEST
	//admm_basic s2("./graph200.txt");
	vector<double>timerec;
	//vector<double>sb{ 1,5,10,15 };
	//vector<double>sa{ 3,5,7,9 };
	//vector<double>st{ 1.618,1.9 };
#ifdef ADMMTEST
	vector<double>sb{ 1 };
	vector<double>sa{ 3,5,10,15 };
	vector<double>st{ 1.618};
	for (auto i:st) {

		for(auto j:sa){
			for (auto k : sb) {
				MATFile* pmatFile = matOpen(string("./random/random500.mat").c_str(), "r");

				admm_basic s2(pmatFile);
				s2.mthread = 32;
				tm.tic();
				double rs = s2.run_worker_simple(i,j,k);
				tm.toc();
				cout << tm.blink << endl;
				timerec.push_back(tm.blink);

				res.push_back(ppdd(pdd(rs, tm.blink), test(s2.x + s2.b, s2.n)));
				cout << i << " " << j << " " << k << endl;
			}
		}
	}
#else
	MATFile* pmatFile = matOpen(string("./random/random500.mat").c_str(), "r");
	metric m(pmatFile);
	norm_ori s2(m);
	s2.run_ori(0.001,30);

#endif
	//for (int i = 0; i < s2.w.size(); i++) {
	//	cout << s2.w[i] << endl;
	//}
	//s2.run();
#endif 

#ifdef MULTITEST
		

	//for (auto i : target) {
		//cout << i << endl;
	double cst_const;
	double min_const_num;
	int max_subiter_num;
#ifdef REL
	fstream tpara("./para.txt");
	try {
		string t1;
		tpara >> t1;
		if (t1.substr(0, 9) == string("cst_const")) {
			cst_const = stod(t1.substr(10));

		}
		else {
			throw t1;
		}
		tpara >> t1;
		if (t1.substr(0, 12) == string("min_cons_num")) {
			min_const_num = stod(t1.substr(13));
		}
		else {
			throw t1;
		}
		tpara >> t1;
		if (t1.substr(0, 15) == string("max_subiter_num")) {
			max_subiter_num = stoi(t1.substr(16));
		}
		else {
			throw t1;
		}

	}
	catch (string msg) {
	
		cout << "parameter error:"<<msg;
		return 0;
	}

#endif
	//44454


	vector<double>vkkt;
	vector<string> target;
	vector<int>grb_iters;
	string input_dir(".");
	int norm,maxits=-1,threadnum;
	bool cgflag = false;
	threadnum = std::thread::hardware_concurrency();
	enum class algorithm { palm, dykstra,grb};
	algorithm mode = algorithm::palm;
	int  statusfreq = 10;
	double contol = 1e-2, gaptol = 1e-2, gam = 1;


	auto palm_mode = (clipp::command("palm").set(mode, algorithm::palm),
		clipp::integer("norm", norm).doc("0 for l_\inf, 1 for l_1, 2 for l_2"),
		clipp::values("filename", target),
		clipp::option("--maxits") & clipp::integer("max_iter", maxits).doc("maximum value for sub iterations, default:30 for l_1 and l_2, 20 for l_inf"),
		clipp::option("--nocg").set(cgflag).doc("do not use constrain generate"),
		clipp::option("--tnum") & clipp::integer("thread_num", threadnum).doc("kkt tolerance of inner iterations, default is the return value of std::thread::hardware_concurrency()")
		);
	auto dykstra_mode = (clipp::command("dykstra").set(mode, algorithm::dykstra),
		clipp::integer("norm", norm).doc("0 for l_\inf, 1 for l_1, 2 for l_2"),
		clipp::values("filename", target),
		clipp::option("--maxits") & clipp::integer("max_iter", maxits).doc("maximum value of iterations, default:10000"),
		clipp::option("-ct", "--contol") & clipp::value("contol", contol).doc("tolerance of violation for constraints, default:0.01"),
		clipp::option("-gt", "--gaptol") & clipp::value("gaptol", gaptol).doc("tolerance of relative gap, default:0.01"),
		clipp::option("-sf", "--statusfreq") & clipp::integer("statusfreq", statusfreq).doc("frequence of checking status and exit, default:10"),
		clipp::option("--gam") & clipp::value("gamma", gam).doc("parameter gamma for l_1 and l_\inf, note that for l_2 gamma is fixed to 0.5, for l_1 and l_\inf default value=1")
		);
	auto gurobi_mode = (clipp::command("grb").set(mode, algorithm::grb),
		clipp::integer("norm", norm).doc("0 for l_\inf, 1 for l_1, 2 for l_2"),
		clipp::values("filename", target),
		clipp::option("--nocg").set(cgflag,true).doc("do not use constrain generate"),
		clipp::required("--maxits")&clipp::values("max_iter", grb_iters)
		);


	auto cli = (
		(palm_mode | dykstra_mode| gurobi_mode),
		clipp::option("-i") & clipp::value("input directory", input_dir)
		);
	pdd testres;
	double rs;
	vector<map<string, double> >v_runhist;

	if (clipp::parse(argc, argv, cli)) {
		switch (mode) {
		case algorithm::grb:

			for (int i = 0; i < target.size();i++) {
				gurobi s2(input_dir + string("/") + target[i] + string(".txt"));
				cout << "cgflag:" << cgflag << endl;
				s2.cgflg = cgflag;
				tm.tic();
				try {
					rs = s2.run(norm,grb_iters[i]);
					testres = test(s2.x, s2.b, s2.n);
					tm.toc();
				}
				catch (...) {
					cout << "memory not enough" << endl;
				}
				res.push_back(ppdd(pdd(rs, tm.blink), testres));
				cout << std::fixed << "total time:" << tm.blink << std::scientific << "\tpobj:" << rs << endl;

			}

			
			
			
			break;

		case algorithm::dykstra:
			if (maxits == -1)maxits = 1000;
			for (auto i : target) {
				cout << "dataset:" << i << endl;
				if (norm == 0) {
					jk_inf s2(input_dir + string("/") + i + string(".txt"));
					s2.gam = gam;
					tm.tic();
					try {
						rs = s2.run(gaptol, contol, statusfreq, maxits);
						testres = test(s2.x, s2.b, s2.n);
						tm.toc();
					}
					catch (...) {
						cout << "memory not enough" << endl;
					}

				}
				else if (norm == 1) {
					jk s2(input_dir + string("/") + i + string(".txt"));
					s2.gam = gam;
					tm.tic();
					try {
						rs = s2.run(gaptol, contol, statusfreq, maxits);
						testres = test(s2.x, s2.b, s2.n);
						tm.toc();
					}
					catch (...) {
						cout << "memory not enough" << endl;
					}

				}
				else if (norm == 2) {
					jk_square s2(input_dir + string("/") + i + string(".txt"));
					tm.tic();
					try {
						rs = s2.run(gaptol, contol, statusfreq, maxits);
						testres = test(s2.x, s2.b, s2.n);
						tm.toc();
					}
					catch (...) {
						cout << "memory not enough" << endl;
					}

				}
				res.push_back(ppdd(pdd(rs, tm.blink), testres));
				cout << std::fixed << "total time:" << tm.blink << std::scientific << "\tpobj:" << rs << endl;
			}
			break;

		case algorithm::palm:
			for (auto i : target) {
				if (norm == 2) {


					norm_2 s2(input_dir + string("/") + i + string(".txt"));
					s2.mthread = threadnum;
					s2.cgflg = cgflag;
					tm.tic();
					if (maxits == -1)maxits = 30;
					try {

						double rs = s2.run_axbu(1, 0, maxits);
						//double rs = s2.run_run_run();
						//double rs = s2.run(gurobi_type,mxstep[i]);
						//double rs = s2.run(1);
						//double rs = s2.ru_run_runn(0);
						tm.toc();
						res.push_back(ppdd(pdd(rs, tm.blink), test(s2.x, s2.b, s2.n)));
						vkkt.push_back(s2.kkt);
						cout << "time:" << tm.blink << endl;
						v_runhist.push_back(s2.runhist);
					}
					catch (...) {
						cout << "memory not enough" << endl;
					}
				}
				else if (norm == 1) {


					norm_ori s2(input_dir + string("/") + i + string(".txt"));
					s2.mthread = threadnum;
					s2.cgflg = cgflag;
					tm.tic();
					if (maxits == -1)maxits = 30;
					try {

						double rs = s2.run_axbu(1, 0, maxits);
						//double rs = s2.run_run_run();
						//double rs = s2.run(gurobi_type,mxstep[i]);
						//double rs = s2.run(1);
						//double rs = s2.ru_run_runn(0);
						tm.toc();
						res.push_back(ppdd(pdd(rs, tm.blink), test(s2.x, s2.b, s2.n)));
						vkkt.push_back(s2.kkt);
						cout << "time:" << tm.blink << endl;
						v_runhist.push_back(s2.runhist);
					}
					catch (...) {
						cout << "memory not enough" << endl;
					}
				}
				else if (norm == 0) {

					if (maxits == -1)maxits = 20;
					norm_inf s2(0,input_dir + string("/") + i + string(".txt"));
					s2.mthread = threadnum;
					s2.cgflg = cgflag;
					tm.tic();

					try {

						double rs = s2.run_axbu_sparse(1, 0, maxits);
						//double rs = s2.run_run_run();
						//double rs = s2.run(gurobi_type,mxstep[i]);
						//double rs = s2.run(1);
						//double rs = s2.ru_run_runn(0);
						tm.toc();
						res.push_back(ppdd(pdd(rs, tm.blink), test(s2.sv_x, s2.sv_b, s2.n)));
						//res.push_back(ppdd(pdd(rs, tm.blink), pdd(-1.0,-1.0)));
						vkkt.push_back(s2.kkt);
						cout << "time:" << tm.blink << endl;
						for (auto i : s2.itervec) {
							cout << i << " ";
						}
						cout << endl;
						v_runhist.push_back(s2.runhist);
					}
					catch (...) {
						cout << "memory not enough" << endl;
					}
				}





			}
		}
	}
	else {
		cout << clipp::make_man_page(cli, "") << endl;
	}

	//}




#endif


	//tm.tic();
	//double tp=s2.run();
	//tm.toc();
	//res.push_back(ppdd(pdd(tp, tm.blink), test(s2.x + s2.b, s2.n)));
	int tnum = 0;
	for (int i = 0; i < res.size();i++) {
		//cout << "res:" << i.first.first<<"\t"<<i.first.second<< "\t" << "sum:" << i.second.first << " " << "max:" << i.second.second << endl;
		if (mode == algorithm::palm) {
			cout << "dataset:" << target[i] << endl;
			cout << "primal:" << res[i].first.first << endl;
			cout << "total_time:" << res[i].first.second << endl;
			for (auto& j : v_runhist[i]) {
				cout << j.first << ":" << j.second << endl;
			}
			cout << "==========\n";
		}
		else {
			cout << "dataset:" << target[i] << "\tprimal:" << res[i].first.first << "\t" << "time:" << res[i].first.second << " kkt" << endl;

		}
		
	}
	cout << "over,input 'exit' to exit";
	string over;
	//while (cin >> over) {
	//	if (over == string("exit"))return 0;
	//}
	//norm_ori s2(s);
	//s2.run_dense();
	//s2.run(0.001, 0, 20);

	//if (s.substr(s.length() - 3, s.length()) == "dat") {
	//	string ssubs = s.substr(s.find('/') + 1);
	//	string subs=ssubs.substr(ssubs.find('_') + 1);
	//	subs = subs.substr(0,subs.find('_'));

	//	cout << ssubs << "\t" << subs << "\t" << endl;
	//	cout << "cst=delta/sigma,plz input delta(do not input 1e-3, use 0.001 instead):";
	//	double cst; int cons;
	//	cin >> cst;
	//	cout << "number of constraints is no less than cons*n*(n-1)/2,plz input cons:";
	//	cin >> cons;

	//	norm_ori s2(subs, s);
	//	s2.run(cst,cons,20);
	//}
	//else {
	// 
		//cout << "cst=delta/sigma,plz input delta(do not input 1e-3, use 0.001 instead):";
		//double cst; int cons; int sub_iter;
		//cin >> cst;
		//cout << "number of constraints is no less than cons*n*(n-1)/2,plz input cons:";
		//cin >> cons;
		//cout << "plz input maxium steps of subiteration:";
		//cin >> sub_iter;
		//norm_ori s2(s);
		//s2.run(cst, cons,sub_iter);


	//}
	//cin >> s;

	//cin >> s;
	// 
	// 
	// 
	//dk s2(s);
	//s2.run_test(1e-3);
	// 
	// 
	// 
	// 
	// 
	// 
	//string dataname;
	//cin >> dataname;
	//norm_2 s(dataname+string(".txt"));
	//s.run();
	//string tem;
	//cin >> tem;
    return 0;
}
