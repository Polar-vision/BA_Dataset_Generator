#pragma once
#include <stdio.h>
#include "myio.h"
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include <Eigen/StdVector>
#include <Eigen/Cholesky>
#define PI  3.1415926535898 
#define MAXARCHOR 0.5
using namespace std;
using namespace Eigen;

class triangulation:public myio
{
public:
	triangulation(int mag);
	~triangulation(void);
	bool ba_initialize(char* szCamera, char* szFeature, char* szCalib, char* szXYZ);
private:
	void pba_angle2xyz(double* p);
	void pba_saveInitialXYZ(const char* sz3Dpt, double* p);
	void pba_saveInitialParallax(const char* sz3Dpt, double* p);
	void pba_readProjectionAndInitilizeFeature(FILE* fp,
		double* params, double* projs, char* vmask, int ncams,
		int* archor, char* umask, int* nphoto, int* nfeature, int* archorSort);
	void pba_readAndInitialize(char* camsfname, char* ptsfname, char* calibfname, int* ncams,
		int* n3Dpts, int* n2Dprojs,
		double** motstruct, double** imgpts,
		int** archor, char** vmask,
		char** umask, int** nphoto,
		int** nfeature, int** archorSort);
	bool pba_initializeOtheArchors(
		double* imgpts,
		int* photo,
		double* camera,
		double* K,
		double* feature,
		int* archorSort,
		int nfeacout,
		int nOI,
		int FID);
	bool pba_initializeAssoArchor(
		double* imgpts,
		int* photo,
		double* camera,
		double* K,
		double* feature,
		int nMI,
		int nAI,
		int FID,
		bool bLast);
	bool pba_initializeMainArchor(
		double* imgpts,
		double* camera,
		double* K,
		double* feature,
		int nP,
		int FID,
		double* KR);
	void ba_constructP(double* P, double* K, double* p);
	void ba_readCameraPose(FILE* fp, double* params);
	void ba_updateKR(double* KR, double* KdA, double* KdB, double* KdG, double* K, double* p);
	void ba_readCameraPoseration(char* fname, double* ical);
	void ba_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams);
	void ba_saveTriangulatedxyz(const char* sz3Dpt, double* p);

	int		m_ncams, m_n3Dpts, m_n2Dprojs, m_nS;  //number of camera, 3D points, 2D projection points, non-zero element of S matrix
	int* m_archor;
	int* m_photo, * m_feature;
	double* m_motstruct, * m_imgpts;			  //6 camera pose and 3 feature parameters/PBA parameter,
	double* m_XYZ;								  //initial XYZ provided 	
	double* m_K;								  //calibration parameters
	int* m_C;//camera id

	char* m_vmask, * m_umask, * m_smask;
	int* m_imgptsSum, * m_struct, * m_pnt2main, * m_archorSort;
	double* m_KR, * m_KdA, * m_KdB, * m_KdG, * m_P;

	bool    m_bProvideXYZ=false, m_bFocal;
	char* m_szCameraInit;
	char* m_szFeatures;
	char* m_szCalibration;
	char* m_szXYZ;
	char* m_szCamePose;
	char* m_sz3Dpts;
	char* m_szReport;
	int		m_nMaxIter;
	double  m_Tau, m_e1, m_e2, m_e3, m_e4;
	bool	m_bRobustKernel;
	bool	m_bsolverLM;
	bool	m_bsolverGN;
	int		m_nRobustType;
	double  m_delt;

	int savePara = 0;
	int zu;
};

