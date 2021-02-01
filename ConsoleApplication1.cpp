// ConsoleApplication1.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include <iostream>
#include <fstream>
#include<vector>
#include<set>
#include<map>
#include<algorithm>
#include <iterator> 
#include<math.h>
#include<numeric>
#include <omp.h>
using namespace std;   
typedef map<int, set<int> > mis;

map<int, set<int> > g;





double gmma =2;
double epsi = 0.01;

void prt(double* x, int var_num) {
    for (int i = 0; i < var_num; i++)
        cout << x[i] << "\t" << x[i + var_num] << endl;
}

inline long long ij(long long x, long long y, long long n) {
    return n * (n - 1) / 2 - n - (n - x - 1) * (n - 2 - x) / 2 + y;
}
inline long long ij2(long long x, long long y, long long n) {
    return n * (n - 1) / 2 - n - (n - x) * (n - 1 - x) / 2 + y - 1;
}
int main()
{

    ifstream graph("graph150.txt", ios::in);
    ofstream prt_res("prt_res.txt", ios::out);
    int tmp_a, tmp_b;
    //build graph
    int input_ct = 0;
    while (graph >> tmp_a >> tmp_b) {
        g[tmp_a].insert(tmp_b);
        g[tmp_b].insert(tmp_a);
        input_ct++;
        //if (input_ct > 300)
        //    break;
    }
    cout << "line:" << input_ct << endl;
    //point list
    vector<int> target;
    //int n= 50;
    for (mis::iterator it = g.begin(); it != g.end(); it++) {
        target.push_back(it->first);
        //if (target.size() >= n)
        //   break;
    }
    int n = target.size();
    cout << "n:" << n << endl;
    int var_num = n * (n - 1) / 2;
    long long cons_num = n * (n - 1) * (n - 2) / 6;
    double* w = new double[var_num];
    std::fill(w, w + var_num, 0);

    bool* b = new bool[var_num];
    double* z = new double[var_num];

    int tmp_it = 0;
    double delta = 0.05;
    //calculate w
    cout << gmma << endl;
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
            w[tmp_it] = 1/abs(tmp);
            //calculate d


        }
    }

    cout << "w:" << accumulate(w, w + var_num, 0.0) << endl;
    cout << "b:" << accumulate(b, b + var_num, 0.0) << endl;
    //map for y
    map<long long, double>y;
    //default x
    double* x = new double[2*var_num];
    fill(x, x + var_num, 0.0);
    fill( x + var_num,x+2*var_num, -1.0*gmma);


    for (int main_i = 0; main_i < 6000; main_i++) {
        long long ct = 0;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                for (int k = j + 1; k < n; k++) {

                    long long idx[3];
                    idx[0] = ij(i, j, n);
                    idx[1] = ij(i, k, n);
                    idx[2] = ij(j, k, n);
                    for (int tmp = 0; tmp < 3; tmp++, ct++) {
                        if (y.find(ct) != y.end())
                        {
                            //update x
                            x[idx[0]] += y[ct] * w[idx[0]] * (tmp == 0 ? 1 : -1);
                            x[idx[1]] += y[ct] * w[idx[1]] * (tmp == 1 ? 1 : -1);
                            x[idx[2]] += y[ct] * w[idx[2]] * (tmp == 2 ? 1 : -1);
                        }
                        //calculate theta

                        int b_ijk = -2 * b[idx[tmp]] + (b[idx[0]] + b[idx[1]] + b[idx[2]]);
                        double ax;
                        ax = 2 * x[idx[tmp]] - (x[idx[0]] + x[idx[1]] + x[idx[2]]);
                        if (ax - b_ijk > 0) {

                            //cout << "wtf\n";
                            double theta = (ax - b_ijk) / ( (w[idx[0]] + w[idx[1]] + w[idx[2]]));
                            //cout << "theta=" << theta << endl;
                            y[ct] = theta;
                            x[idx[0]] -= y[ct] *  w[idx[0]] * (tmp == 0 ? 1 : -1);
                            x[idx[1]] -= y[ct] *  w[idx[1]] * (tmp == 1 ? 1 : -1);
                            x[idx[2]] -= y[ct] *  w[idx[2]] * (tmp == 2 ? 1 : -1);
                        }
                        else {
                            if (y.find(ct) != y.end())
                                y.erase(ct);
                        }
                    }
                }
            
            }
        }

        for (int i = 0; i < var_num; i++, ct++) {

            if (y.find(ct) != y.end()) {
                //update x
                x[i] += y[ct]  * w[i];
                x[i + var_num] -= y[ct] *  w[i];
            }
            if (x[i] - x[i + var_num] >= 0) {

                //update y
                y[ct] = (x[i] - x[i + var_num]) / ( (2 * w[i]));
                //correct x
                x[i] -= y[ct] * w[i];
                x[i + var_num] += y[ct]  * w[i];
            }
            else {
                if (y.find(ct) != y.end())
                    y.erase(ct);
            }
        }

        for (int i = 0; i < var_num; i++, ct++) {

            if (y.find(ct) != y.end()) {
                //update x
                x[i] -= y[ct]  * w[i];
                x[i + var_num] -= y[ct]  * w[i];
            }
            if (-x[i] - x[i + var_num] >= 0) {
                //update y


                y[ct] = (-x[i] - x[i + var_num]) / ( (2 * w[i]));
                //correct x
                x[i] += y[ct] * w[i];
                x[i + var_num] += y[ct]  * w[i];
            }
            else {
                if (y.find(ct) != y.end())
                    y.erase(ct);
            }
        }
        
        cout << "time:" << main_i << endl;

        if (main_i % 100 == 0) {
            bool flag = true;
            for (int i = 0; i < n && flag; i++) {
                for (int j = i + 1; j < n && flag; j++) {
                    for (int k = j + 1; k < n && flag; k++) {

                        long long idx[3];
                        idx[0] = ij(i, j, n);
                        idx[1] = ij(i, k, n);
                        idx[2] = ij(j, k, n);
                        for (int tmp = 0; tmp < 3; tmp++, ct++) {

                            //calculate theta

                            int b_ijk = -2 * b[idx[tmp]] + (b[idx[0]] + b[idx[1]] + b[idx[2]]);
                            double ax;
                            ax = 2 * x[idx[tmp]] - (x[idx[0]] + x[idx[1]] + x[idx[2]]);
                            if (ax - b_ijk > 1e-6) {
                                flag = false;
                                goto X;
                            }
                        }
                    }
                }
            }
            X:
            if (flag)
                break;
        }
    }

    double res = 0;
    for (int i = 0; i < var_num; i++) {
        res += abs(x[i + var_num])/w[i];

    }
    cout << "res:"<<res << endl;

    delete []b;
    delete []w;
    cout << "succeed" << endl;

    delete[]x;
    return 0;
}

// 运行程序: Ctrl + F5 或调试 >“开始执行(不调试)”菜单
// 调试程序: F5 或调试 >“开始调试”菜单

// 入门使用技巧: 
//   1. 使用解决方案资源管理器窗口添加/管理文件
//   2. 使用团队资源管理器窗口连接到源代码管理
//   3. 使用输出窗口查看生成输出和其他消息
//   4. 使用错误列表窗口查看错误
//   5. 转到“项目”>“添加新项”以创建新的代码文件，或转到“项目”>“添加现有项”以将现有代码文件添加到项目
//   6. 将来，若要再次打开此项目，请转到“文件”>“打开”>“项目”并选择 .sln 文件
