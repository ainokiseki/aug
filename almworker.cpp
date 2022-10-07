#include"almworker.h"



alm::alm(std::vector<int> &colidx) :colidx(colidx)
{

}

alm::~alm()
{
}


void alm::cons_generate_admm(Eigen::VectorXd& x, Eigen::VectorXd& b, int& colidx_num, std::vector<Eigen::Triplet<double> >&acons,std::set<int64_t>&vviscons,Eigen::VectorXd&diag, vector<alm_total*>&worker, ReaderWriterLatch&rw,double limit) {
	for (int i = istart; i < iend; i++) {
		int64_t index;
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				index = static_cast<int64_t>(i) * n * n + j * n + k;

				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 3e-2) {
					rw.RLock();
					bool flg = (vviscons.find(index) == vviscons.end());
					rw.RUnlock();
					if (flg) {
						mu->lock();
						for (int tp = 0; tp < 3; tp++) {
							if (colidx[idx[tp]] == -1) {

								colidx[idx[tp]] = colidx_num;
								colidx_num++;

							}
						}


						

						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[0]], 1));
						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[1]], -1));
						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[2]], -1));

						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[0]], -1));
						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[1]], 1));
						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[2]], -1));

						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[0]], -1));
						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[1]], -1));
						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[2]], 1));

						acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[2]], -1.0));
						acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[1]], -1.0));


						acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[2]], -1.0));
						acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[0]], -1.0));


						acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[1]], -1.0));
						acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[0]], -1.0));

						diag[colidx[idx[0]]] += 3; diag[colidx[idx[1]]] += 3; diag[colidx[idx[2]]] += 3;

						mu->unlock();

						rw.WLock();
						worker[vviscons.size() % worker.size()]->cons.push_back(index);
						vviscons.insert(index);
						rw.WUnlock();
					}
				}
			}
		}
	}
}
void alm::check_violence(Eigen::VectorXd & sub_x,Eigen::VectorXd&sub_b,atomic<bool> &flag,std::pair<double,int64_t>&vionum,double limit) {

	int uidx = 0;
	for (int istart = 0; flag.load()&&istart < cons.size(); istart++) {

		int64_t lindex = cons[istart];
		int i, j, k;
		k = lindex % n;
		j = lindex / n;
		i = j / n;
		j = j % n;

		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];

		if (sub_x[idx[0]] + sub_b[idx[0]] - sub_x[idx[1]] - sub_b[idx[1]] - sub_x[idx[2]] - sub_b[idx[2]] > limit ||
			sub_x[idx[1]] + sub_b[idx[1]] - sub_x[idx[0]] - sub_b[idx[0]] - sub_x[idx[2]] - sub_b[idx[2]] > limit ||
			sub_x[idx[2]] + sub_b[idx[2]] - sub_x[idx[1]] - sub_b[idx[1]] - sub_x[idx[0]] - sub_b[idx[0]] > limit) {
			//mu->lock();
			//vionum.second = lindex;
			//vionum.first = std::fmax(sub_x[idx[0]] + sub_b[idx[0]] - sub_x[idx[1]] - sub_b[idx[1]] - sub_x[idx[2]] - sub_b[idx[2]], sub_x[idx[1]] + sub_b[idx[1]] - sub_x[idx[0]] - sub_b[idx[0]] - sub_x[idx[2]] - sub_b[idx[2]]);
			//vionum.first = std::fmax(vionum.first, sub_x[idx[2]] + sub_b[idx[2]] - sub_x[idx[1]] - sub_b[idx[1]] - sub_x[idx[0]] - sub_b[idx[0]]);
			//cout << vionum.first << "  " << vionum.second << endl;
			//mu->unlock();
			
			flag.store(false);
			return;

		}



	}


}
void alm::cons_generate(Eigen::VectorXd& x, Eigen::VectorXd& b, int& colidx_num) {
	for (int i = istart; i < iend; i++) {
		int64_t index;
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				index = static_cast<int64_t>(i) * n * n + j * n + k;

				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 1e-4) {

					if (viscons.find(index) == viscons.end()) {

						for (int tp = 0; tp < 3; tp++) {
							mu->lock();
							if (colidx[idx[tp]] == -1) {

								colidx[idx[tp]] = colidx_num;
								colidx_num++;

							}
							mu->unlock();
						}


						cons.push_back(index);
						viscons.insert(index);

					}
				}
			}
		}
	}
}




void alm::axbu(func_eigen& f) {

	int unum = 0;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int istart = 0; istart < cons.size(); istart++) {
		int64_t longindex = cons[istart];
		int i, j, k;
		k = longindex % n;
		j = longindex / n;
		i = j / n;
		j = j % n;
		//std::cout << i << "\t" << j << "\t" << k <<std::endl;


		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];


		double axbu = 0;
		for (int time = 0; time < 3; time++) {
			axbu -= (*x)[idx[time]] + (*b)[idx[time]];
		}
		for (int time = 0; time < 3; time++) {
			int64_t index;
			switch (time)
			{
			case 0:
				index = longindex;
				break;
			case 1:
				index = static_cast<int64_t>(i) * n * n + k * n + j;
				break;
			case 2:
				index = static_cast<int64_t>(j) * n * n + k * n + i;
			}

			double axbu2 = 0;

			if (unum < u.size() && u[unum].first == index) {
				axbu2 += u[unum].second / sigma;
				unum++;
			}
			for (int tt = 0; tt < 3; tt++) {
				if (tt == time)
					axbu2 += ((*x)[idx[tt]] + (*b)[idx[tt]]);
				else
					axbu2 -= ((*x)[idx[tt]] + (*b)[idx[tt]]);
			}
			//axbu2 += 2 * ((*x)[idx[time]] + (*b)[idx[time]]);
			f(axbu2, idx, time, index);
		}



	}
	f.integrate();
}



void alm::atu(Eigen::VectorXd* au)
{

	int unum = 0;
	std::vector< std::pair<int64_t, double> >& tu(u);
	int var_num = n * (n - 1) / 2;
	//int64_t index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int istart = 0; istart < cons.size(); istart++) {
		int64_t longindex = cons[istart];
		int i, j, k;
		k = longindex % n;
		j = longindex / n;
		i = j / n;
		j = j % n;


		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];

		for (int time = 0; time < 3; time++) {
			int64_t index;
			switch (time)
			{
			case 0:
				index = longindex;
				break;
			case 1:
				index = static_cast<int64_t>(i) * n * n + k * n + j;
				break;
			case 2:
				index = static_cast<int64_t>(j) * n * n + k * n + i;
			}
			if (unum >= tu.size()) {
				break;
			}
			if (tu[unum].first != index)
				continue;
			double tem = tu[unum].second;
			(*au)[idx[0]] -= tem;
			(*au)[idx[1]] -= tem;
			(*au)[idx[2]] -= tem;
			(*au)[idx[time]] += 2 * tem;
			unum++;
		}


	}

}
void alm::all_in_all(Eigen::VectorXd& x, Eigen::VectorXd& b, Eigen::VectorXd& atu, Eigen::VectorXd& old_x, double& axby, double& y_norm, double sigma,  double tau, double scale_a) {
	std::vector< std::pair<int64_t, double> > unew;
	auto square = [](double x) {return x * x; };
	std::vector<std::pair<int64_t, double> >ins;
	std::vector<std::map<int64_t, double>::iterator>del;

	int uidx = 0;
	for (int istart = 0; istart < cons.size(); istart++) {

		int64_t lindex = cons[istart];
		int i, j, k;
		k = lindex % n;
		j = lindex / n;
		i = j / n;
		j = j % n;

		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];


		int64_t longindex[3];
		longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
		longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
		longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;
		vector<double> restar(3);

		{
			double u_value = 0;

			if (uidx < u.size() && longindex[0] == u[uidx].first) {
				u_value = u[uidx].second;
				uidx++;
			}

			double axb, axbu;

			axbu = (old_x[idx[0]] + b[idx[0]] - old_x[idx[1]] - b[idx[1]] - old_x[idx[2]] - b[idx[2]]) / scale_a;
			//calc y,ax
			axb = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;


			axbu += u_value / sigma;
			if (axbu > 0)
				axbu = 0;
			y_norm += square(axbu);
			axb -= axbu;
			axby += square(axb);
			double u_new = axb * tau * sigma + u_value;

			if (u_new != 0) {

				unew.push_back(std::pair<int64_t, double>(longindex[0], u_new));
			}
			restar[0] = u_new / scale_a;


		}

		{
			double u_value = 0;
			if (uidx < u.size() && longindex[1] == u[uidx].first) {
				u_value = u[uidx].second;
				uidx++;
			}

			double axb, axbu;


			//calc y,ax

			axbu = (old_x[idx[1]] + b[idx[1]] - old_x[idx[0]] - b[idx[0]] - old_x[idx[2]] - b[idx[2]]) / scale_a;
			//calc y,ax
			axb = (x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]]) / scale_a;


			axbu += u_value / sigma;
			if (axbu > 0)
				axbu = 0;
			y_norm += square(axbu);
			axb -= axbu;
			axby += square(axb);
			double u_new = axb * tau * sigma + u_value;


			if (u_new != 0) {

				unew.push_back(std::pair<int64_t, double>(longindex[1], u_new));
			}
			restar[1] = u_new / scale_a;


		}

		{
			double u_value = 0;
			if (uidx < u.size() && longindex[2] == u[uidx].first) {
				u_value = u[uidx].second;
				uidx++;
			}

			double axb, axbu;


			//calc y,ax
			axbu = (old_x[idx[2]] + b[idx[2]] - old_x[idx[1]] - b[idx[1]] - old_x[idx[0]] - b[idx[0]]) / scale_a;
			//calc y,ax
			axb = (x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]]) / scale_a;


			axbu += u_value / sigma;
			if (axbu > 0)
				axbu = 0;
			y_norm += square(axbu);
			axb -= axbu;
			axby += square(axb);
			double u_new = axb * tau * sigma + u_value;

			if (u_new != 0) {

				unew.push_back(std::pair<int64_t, double>(longindex[2], u_new));
			}
			restar[2] = u_new / scale_a;


		}

		mu[idx[0] % n].lock();
		atu[idx[0]] += restar[0] - restar[1] - restar[2];
		mu[idx[0] % n].unlock();
		mu[idx[1] % n].lock();
		atu[idx[1]] += restar[1] - restar[0] - restar[2];
		mu[idx[1] % n].unlock();
		mu[idx[2] % n].lock();
		atu[idx[2]] += restar[2] - restar[1] - restar[0];
		mu[idx[2] % n].unlock();



	}
	u = std::move(unew);
	unew.~vector();



}
void alm::atyu( Eigen::VectorXd& x, Eigen::VectorXd& b,  Eigen::VectorXd& res, double sigma,  double scale_a ) {

	int uidx = 0;
	for (int istart = 0; istart < cons.size(); istart++) {
		vector<double> restar(3);
		int64_t lindex = cons[istart];
		int i, j, k;
		k = lindex % n;
		j = lindex / n;
		i = j / n;
		j = j % n;

		int idx[3];
		idx[0] = colidx[ij(i, j, n)];
		idx[1] = colidx[ij(i, k, n)];
		idx[2] = colidx[ij(j, k, n)];

		double axbu;

		double flag;
		double tem;
		int64_t longindex[3];
		longindex[0] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(j) * n + k;
		longindex[1] = static_cast<int64_t>(i) * n * n + static_cast<int64_t>(k) * n + j;
		longindex[2] = static_cast<int64_t>(j) * n * n + static_cast<int64_t>(k) * n + i;

		axbu = (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]]) / scale_a;
		if (uidx < u.size() && u[uidx].first == longindex[0]) {
			flag = u[uidx].second;
			uidx++;
		}
		else {
			flag = 0;
		}




		axbu += flag / sigma;



		if (axbu > 0) {
			axbu = 0;
		}
		tem = axbu - (b[idx[0]] - b[idx[1]] - b[idx[2]]) / scale_a;

		tem -= flag / sigma;

		restar[0] = tem / scale_a;
		//res[idx[0]] += tem / scale_a;
		//res[idx[1]] -= tem / scale_a;
		//res[idx[2]] -= tem / scale_a;

		//second scene
		axbu = (x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]]) / scale_a;
		if (uidx < u.size() && u[uidx].first == longindex[1]) {
			flag = u[uidx].second;
			uidx++;
		}
		else {
			flag = 0;
		}
		axbu += flag / sigma;

		if (axbu > 0) {
			axbu = 0;
		}
		tem = axbu - (b[idx[1]] - b[idx[0]] - b[idx[2]]) / scale_a;

		tem -= flag / sigma;

		restar[1] = tem / scale_a;

		//third scene
		axbu = (x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]]) / scale_a;
		if (uidx < u.size() && u[uidx].first == longindex[2]) {
			flag = u[uidx].second;
			uidx++;
		}
		else {
			flag = 0;
		}
		axbu += flag / sigma;

		if (axbu > 0) {
			axbu = 0;
		}
		tem = axbu - (b[idx[2]] - b[idx[1]] - b[idx[0]]) / scale_a;

		tem -= flag / sigma;

		restar[2] = tem / scale_a;
		mu[idx[0] % n].lock();
		res[idx[0]] += restar[0] - restar[1] - restar[2];
		mu[idx[0] % n].unlock();
		mu[idx[1] % n].lock();
		res[idx[1]] += restar[1] - restar[0] - restar[2];
		mu[idx[1] % n].unlock();
		mu[idx[2] % n].lock();
		res[idx[2]] += restar[2] - restar[1] - restar[0];
		mu[idx[2] % n].unlock();

	}
}

alm_small::alm_small() 
{

}

alm_small::~alm_small()
{
}




void alm_small::cons_generate(Eigen::VectorXd& x, Eigen::VectorXd& b,int& colidx_num) {
	for (int i = istart; i < iend; i++) {
		int64_t index;
		for (int j = i + 1; j < n; j++) {
			for (int k = j + 1; k < n; k++) {
				int idx[3];
				idx[0] = ij(i, j, n);
				idx[1] = ij(i, k, n);
				idx[2] = ij(j, k, n);
				index = static_cast<int64_t>(i) * n * n + j * n + k;

				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 1e-4 ||
					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 1e-4) {

					if (viscons.find(index) == viscons.end()) {
						cons.push_back(index);
						viscons.insert(index);

					}
				}
			}
		}
	}



}











void alm_small::axbu(func_eigen& f) {

	int unum = 0;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int istart = 0; istart < cons.size(); istart++) {
		int64_t longindex = cons[istart];
		int i, j, k;
		k = longindex % n;
		j = longindex / n;
		i = j / n;
		j = j % n;
		//std::cout << i << "\t" << j << "\t" << k <<std::endl;


		int idx[3];
		idx[0] = ij(i, j, n);
		idx[1] = ij(i, k, n);
		idx[2] = ij(j, k, n);


		double axbu = 0;
		for (int time = 0; time < 3; time++) {
			axbu -= (*x)[idx[time]] + (*b)[idx[time]];
		}
		for (int time = 0; time < 3; time++) {
			int64_t index;
			switch (time)
			{
			case 0:
				index = longindex;
				break;
			case 1:
				index = static_cast<int64_t>(i) * n * n + k * n + j;
				break;
			case 2:
				index = static_cast<int64_t>(j) * n * n + k * n + i;
			}

			double axbu2 = 0;

			if (unum < u.size() && u[unum].first == index) {
				axbu2 += u[unum].second / sigma;
				unum++;
			}
			for (int tt = 0; tt < 3; tt++) {
				if (tt == time)
					axbu2 += ((*x)[idx[tt]] + (*b)[idx[tt]]);
				else
					axbu2 -= ((*x)[idx[tt]] + (*b)[idx[tt]]);
			}
			//axbu2 += 2 * ((*x)[idx[time]] + (*b)[idx[time]]);
			f(axbu2, idx, time, index);
		}



	}
	f.integrate();
}



void alm_small::atu(Eigen::VectorXd* au)
{

	int unum = 0;
	std::vector< std::pair<int64_t, double> >& tu(u);
	int var_num = n * (n - 1) / 2;
	//int64_t index = n * (n - 1) * (n - 2) / 2 - (n - init) * (n - init - 1) * (n - init - 2) / 2;
	//std::cout << init << "\t" << end << "\t" << index << std::endl;
	for (int istart = 0; istart < cons.size(); istart++) {
		int64_t longindex = cons[istart];
		int i, j, k;
		k = longindex % n;
		j = longindex / n;
		i = j / n;
		j = j % n;


		int idx[3];
		idx[0] = ij(i, j, n);
		idx[1] = ij(i, k, n);
		idx[2] = ij(j, k, n);

		for (int time = 0; time < 3; time++) {
			int64_t index;
			switch (time)
			{
			case 0:
				index = longindex;
				break;
			case 1:
				index = static_cast<int64_t>(i) * n * n + k * n + j;
				break;
			case 2:
				index = static_cast<int64_t>(j) * n * n + k * n + i;
			}
			if (unum >= tu.size()) {
				break;
			}
			if (tu[unum].first != index)
				continue;
			double tem = tu[unum].second;
			(*au)[idx[0]] -= tem;
			(*au)[idx[1]] -= tem;
			(*au)[idx[2]] -= tem;
			(*au)[idx[time]] += 2 * tem;
			unum++;
		}


	}

}
//
//
//void alm_small::cons_generate_admm(Eigen::VectorXd& x, Eigen::VectorXd& b, std::vector<Eigen::Triplet<double> >& acons, std::set<int64_t>& vviscons, Eigen::VectorXd& diag, vector<alm>& worker, ReaderWriterLatch& rw, double limit) {
//	for (int i = istart; i < iend; i++) {
//		int64_t index;
//		for (int j = i + 1; j < n; j++) {
//			for (int k = j + 1; k < n; k++) {
//				int idx[3];
//				idx[0] = ij(i, j, n);
//				idx[1] = ij(i, k, n);
//				idx[2] = ij(j, k, n);
//				index = static_cast<int64_t>(i) * n * n + j * n + k;
//
//				if (x[idx[0]] + b[idx[0]] - x[idx[1]] - b[idx[1]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
//					x[idx[1]] + b[idx[1]] - x[idx[0]] - b[idx[0]] - x[idx[2]] - b[idx[2]] > 3e-2 ||
//					x[idx[2]] + b[idx[2]] - x[idx[1]] - b[idx[1]] - x[idx[0]] - b[idx[0]] > 3e-2) {
//					rw.RLock();
//					bool flg = (vviscons.find(index) == vviscons.end());
//					rw.RUnlock();
//					if (flg) {
//
//						mu->lock();
//						worker[viscons.size() % worker.size()].cons.push_back(index);
//
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[0]], 1));
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[1]], -1));
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3, colidx[idx[2]], -1));
//
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[0]], -1));
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[1]], 1));
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 1, colidx[idx[2]], -1));
//
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[0]], -1));
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[1]], -1));
//						//acons.push_back(Eigen::Triplet<double>(viscons.size() * 3 + 2, colidx[idx[2]], 1));
//
//						acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[2]], -1.0));
//						acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[1]], -1.0));
//
//
//						acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[2]], -1.0));
//						acons.push_back(Eigen::Triplet<double>(colidx[idx[2]], colidx[idx[0]], -1.0));
//
//
//						acons.push_back(Eigen::Triplet<double>(colidx[idx[0]], colidx[idx[1]], -1.0));
//						acons.push_back(Eigen::Triplet<double>(colidx[idx[1]], colidx[idx[0]], -1.0));
//
//						diag[colidx[idx[0]]] += 3; diag[colidx[idx[1]]] += 3; diag[colidx[idx[2]]] += 3;
//
//						mu->unlock();
//
//						rw.WLock();
//
//						vviscons.insert(index);
//						rw.WUnlock();
//					}
//				}
//			}
//		}
//	}
//}