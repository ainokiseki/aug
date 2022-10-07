#include<fstream>
#include<iostream>
#include<concurrent_vector.h>
#include"metric.h"
#include"norm_class.h"



static void save4(void* v, ofstream& fout) {
	fout.write((char*)v, sizeof(int));
}
static void save8(void* v, ofstream& fout) {
	fout.write((char*)v, sizeof(double));
}

template<typename T>
static void saveVector(concurrency::concurrent_vector<T>& v, ofstream& fout) {
	int sz = v.size();
	fout.write((char*)(&sz), sizeof(int));
	fout.write((char*)v.data(), sizeof(T) * sz);
}
template<typename T>
static void saveVector(vector<T>& v, ofstream& fout) {
	int sz = v.size();
	fout.write((char*)(&sz), sizeof(int));
	fout.write((char*)v.data(), sizeof(T) * sz);
}

template<typename T>
static void saveSet(set<T>& s,ofstream&fout) {
	int sz = s.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i), sizeof(T));
	}
}
template<typename T>
static void saveSet(unordered_set<T>& s, ofstream& fout) {
	int sz = s.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i), sizeof(T));
	}
}

template<typename T1,typename T2>
static void saveMap(map<T1,T2>& m, ofstream& fout) {
	int sz = m.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i.first), sizeof(T1));
		fout.write((char*)(&i.second), sizeof(T2));
	}
}
template<typename T1, typename T2>
static void saveMap(unordered_map<T1, T2>& m, ofstream& fout) {
	int sz = m.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i.first), sizeof(T1));
		fout.write((char*)(&i.second), sizeof(T2));
	}
}


void palmSparseSave(metric&m, cons_generator_sparse&c,string s) {
	ofstream fout(s, ios::binary);
	save4(&m.sv_b.default_value,fout);
	saveSet(m.sparse_b, fout);

	save8(&m.sv_x.default_value, fout);
	saveMap(m.sparse_x, fout);
	save8(&m.sv_v.default_value, fout);
	saveMap(m.sparse_v, fout);
	save8(&m.sv_w.default_value, fout);
	saveMap(m.sparse_w, fout);
}

