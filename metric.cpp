#include"metric.h"



metric::metric(std::string graph_name) :graph_name(graph_name)
{
	ifstream graph(graph_name, ios::in);
	map<int, set<int> > g;
	int tmp_a, tmp_b;
	int input_ct = 0;
	while (graph >> tmp_a >> tmp_b) {
		g[tmp_a].insert(tmp_b);
		g[tmp_b].insert(tmp_a);
		input_ct++;
	}

	cout << "line:" << input_ct << endl;
	//point list
	vector<int> target;
	//int n= 50;
	for (auto it = g.begin(); it != g.end(); it++) {
		target.push_back(it->first);
	}
	n = target.size();
	cout << "n:" << n << endl;
	var_num = n * (n - 1) / 2;
	cons_num = n * (n - 1) * (n - 2) / 6;
	w = Eigen::VectorXd::Zero(var_num);
	b = Eigen::VectorXd::Zero(var_num);

	int tmp_it = 0;
	const double delta = 0.05;
	const double epsi = 0.01;
	//calculate w
	for (int i = 0; i < n; i++)
	{
		int t_i = target[i];
		for (int j = i + 1; j < n; j++, tmp_it++) {
			//calculate w
			set<int>tmp1;
			set<int>tmp2;
			int t_j = target[j];
			set_intersection(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp1, tmp1.begin()));
			set_union(g[t_i].begin(), g[t_i].end(), g[t_j].begin(), g[t_j].end(), inserter(tmp2, tmp2.begin()));
			double tmp = 1.0 * tmp1.size() / tmp2.size();
			tmp = log((1 + tmp - delta) / (1 - tmp + delta));
			if (tmp > 0) {
				tmp += epsi;
			}
			else if (tmp < 0) {
				tmp -= epsi;
			}
			else {
				if (g[t_i].find(t_j) != g[t_i].end())
					tmp = epsi;
				else
					tmp = -epsi;
			}
			b[tmp_it] = (tmp < 0);
			w[tmp_it] = abs(tmp);
			//calculate d
		}
	}

	cout << "w:" << w.sum() << endl;
	cout << "b:" << b.sum() << endl;
	v = Eigen::VectorXd::Zero(var_num);
	x = Eigen::VectorXd::Zero(var_num);
	cout << "x,v" << x[0] << "\t" << v[0] << endl;
}

metric::metric(int n, Eigen::VectorXd b, Eigen::VectorXd w):n(n),b(b),w(w)
{
	var_num = n * (n - 1) / 2;
	x = Eigen::VectorXd::Zero(var_num);
	v = Eigen::VectorXd::Zero(var_num);
}

metric::~metric()
{
}

void metric::init()
{


	return;
}
