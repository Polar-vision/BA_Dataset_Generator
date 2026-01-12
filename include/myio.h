#pragma once
#include "mycommon.h"
#include <vector>
class myio : public mycommon
{
public:
    //myio(int das);
    //~myio(void);
    void ba_readCameraPose(FILE* fp, FILE* fpw, vector<int> deleteid, double* motcam);
    vector<int> readf(const char* fname, const char* wname, double* ical, int& ncams);
    vector<int> readProjection(FILE* fp, FILE* fpw, vector<int> deleteid, double* motz, double* motcam);
};

