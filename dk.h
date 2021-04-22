#pragma once
#include "metric.h"

class dk : public metric
{
public:
	dk(string graph_name);
	~dk();

	void run();


};

dk::dk(string graph_name) :metric(graph_name)
{
}

dk::~dk()
{
}

void dk::run()
{

}
