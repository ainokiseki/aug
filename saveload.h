#pragma once
#include<fstream>
#include<iostream>
#include<concurrent_vector.h>
#include"metric.h"
#include"norm_class.h"



void save4(void* v, ofstream& fout) {
	fout.write((char*)v, sizeof(int));
}
void save8(void* v, ofstream& fout) {
	fout.write((char*)v, sizeof(double));
}
void read4(void* v, ifstream& fin) {
	fin.read((char*)v, sizeof(int));
}
void read8(void* v, ifstream& fin) {
	fin.read((char*)v, sizeof(double));
}

template<typename T>
void saveVector(concurrency::concurrent_vector<T>& v, ofstream& fout) {
	int sz = v.size();
	fout.write((char*)(&sz), sizeof(int));
	fout.write((char*)v.data(), sizeof(T) * sz);
}
template<typename T>
void readVector(concurrency::concurrent_vector<T>& v, ifstream& fin) {
	int sz;
	read4(&sz,fin);
	v.resize(sz);
	fin.read((char*)v.data(),sizeof(T)*sz);
}
template<typename T>
void saveVector(vector<T>& v, ofstream& fout) {
	int sz = v.size();
	fout.write((char*)(&sz), sizeof(int));
	fout.write((char*)v.data(), sizeof(T) * sz);
}
template<typename T>
void readVector(vector<T>& v, ifstream& fin) {
	int sz;
	read4(&sz, fin);
	v.resize(sz);
	fin.read((char*)v.data(), sizeof(T) * sz);
}


template<typename T>
void saveSet(set<T>& s, ofstream& fout) {
	int sz = s.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i), sizeof(T));
	}
}
template<typename T>
void readSet(set<T>& s, ifstream& fin) {
	int sz;
	read4(&sz, fin);
	T tmp;
	for (int i = 0; i < sz; i++) {
		fin.read((char*)&tmp, sizeof(T));
		s.insert(tmp);
	}
}
template<typename T>
void saveSet(unordered_set<T>& s, ofstream& fout) {
	int sz = s.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i), sizeof(T));
	}
}
template<typename T>
void readSet(unordered_set<T>& s, ifstream& fin) {
	int sz;
	read4(&sz, fin);
	T tmp;
	for (int i = 0; i < sz; i++) {
		fin.read((char*)&tmp, sizeof(T));
		s.insert(tmp);
	}
}

template<typename T1, typename T2>
void saveMap(map<T1, T2>& m, ofstream& fout) {
	int sz = m.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i.first), sizeof(T1));
		fout.write((char*)(&i.second), sizeof(T2));
	}
}
template<typename T1, typename T2>
void readMap(map<T1, T2>& m, ifstream& fin) {
	int sz;
	read4(&sz, fin);
	T1 tmp1; T2 tmp2;
	for (int i = 0; i < sz; i++) {
		fin.read((char*)tmp1, sizeof(T1));
		fin.read((char*)tmp2, sizeof(T2));
		m.insert(pair<T1, T2>(tmp1, tmp2));
	}
}
template<typename T1, typename T2>
void saveMap(unordered_map<T1, T2>& m, ofstream& fout) {
	int sz = m.size();
	fout.write((char*)(&sz), sizeof(int));
	for (auto i : s) {
		fout.write((char*)(&i.first), sizeof(T1));
		fout.write((char*)(&i.second), sizeof(T2));
	}
}
template<typename T1, typename T2>
void readMap(unordered_map<T1, T2>& m, ifstream& fin) {
	int sz;
	read4(&sz, fin);
	T1 tmp1; T2 tmp2;
	for (int i = 0; i < sz; i++) {
		fin.read((char*)tmp1, sizeof(T1));
		fin.read((char*)tmp2, sizeof(T2));
		m.insert(pair<T1, T2>(tmp1, tmp2));
	}
}
