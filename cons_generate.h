#pragma once
#include<vector>
#include<concurrent_unordered_set.h>
#include<concurrent_vector.h>
#include<concurrent_priority_queue.h>
#include<unordered_set>
#include<atomic>
#include<set>
#include"mytime.h"
#include"tri.h"
#define NOEXIT
//#define FULLCONS
//#define MDEBUG

class cons_full
{
public:
	cons_full(int a1, int a2, int a3,int64_t index=0,float u_value=0):index(index) {
		idx[0] = a1; idx[1] = a2; idx[2] = a3;
		u = u_value;
		y =0;
	}
	~cons_full() = default;
	int idx[3];
	float u;
	double y;
	int64_t index;
private:

};



