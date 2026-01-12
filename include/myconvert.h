#pragma once
#include "myio.h"
#include <Eigen/Dense>

class myconvert : public myio
{
public:
	myconvert(int);
	~myconvert(void);
	void bal2pba();
	void outlierdet();
	int colmap2pba(string p1r, string p2r, string p3r, string p1w, string p2w, string p3w, string p4w,\
		double* sigobs, double* sigeul, double* sigpos, double* sigxyz, double* sigcal);
	int ncams = 0, n3Dpts = 0, n2Dprojs = 0;
};

