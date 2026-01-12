#pragma once
#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <iostream>
#include <random>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <opencv2\opencv.hpp>
#include <dataPath.h>

using namespace std;
using namespace Eigen;
#define MAXSTRLEN  2048 /* 2K */
#define SKIP_LINE(f){                                                       \
	char buf[MAXSTRLEN];                                                        \
	while(!feof(f))                                                           \
	if(!fgets(buf, MAXSTRLEN-1, f) || buf[strlen(buf)-1]=='\n') break;      \
}

class mycommon
{
	//mycommon();
	//~mycommon();
public:
	void undist(std::vector<cv::Mat> K, std::vector<cv::Mat> distCoeffs, \
		std::vector<cv::Point2d> distortedPoints, \
		std::vector<cv::Point2d> & undistortedPixelPoints, \
		vector<int> camid_group);
	int readNInts(FILE* fp, int* vals, int nvals);
	int readNDoubles(FILE* fp, double* vals, int nvals);
	void QuaternionToAngleAxis(double* quaternion, double* angle_axis);
	void RotationMatrixToQuaternion(double* R, double* quaternion);
	void rotationMatrixToEulerAngles(double* R, double* eulerAngles);
	void rotationMatrixToEulerAngles_phi_omega_kappa(double *R, double *eulerAngles);
	void eulerAnglesToRotationMatrix(double* eulerAngles, double* R);
	void eulerAnglesToRotationMatrix_phi_omega_kappa(double* eulerAngles, double* R);
	void AngleAxisRotatePoint(double angle_axis[3], double pt[3], double result[3]);
	int findNcameras(FILE* fp);
	void PLY();
	int countNDoubles(FILE* fp);
	void readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp);
	double* qt2euc(double a, double b, double c, double d, double t1, double t2, double t3);
	double deg2rad(double deg);
	double rad2deg(double rad);
	Matrix3d ExpSO3(const Vector3d& phi);
	Vector3d LogSO3(const Matrix3d& R);
	bool isRotationMatrix(const Eigen::Matrix3d& R);
	typedef enum Data_Type
	{
		bal = 1,
		colmap = 2,
	}DataType;
	struct eop8 {
		int viewid;
		double eop[7];
		// ���캯��
		eop8(int d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8)
			: viewid(d1) {
			eop[0] = d2;
			eop[1] = d3;
			eop[2] = d4;
			eop[3] = d5;
			eop[4] = d6;
			eop[5] = d7;
			eop[6] = d8;
		}
	};
	struct eop8_ {
		int viewid;
		int camid;
		double eop[6];
		// ���캯��
		eop8_(int d1, double d2, double d3, double d4, double d5, double d6, double d7, int d8)
			: viewid(d1),camid(d8) {
			eop[0] = d2;
			eop[1] = d3;
			eop[2] = d4;
			eop[3] = d5;
			eop[4] = d6;
			eop[5] = d7;
		}
	};
	struct obs5 {
		int viewid;
		int uvid;
		double uv[2];
		int ptid;
		// ���캯��
		obs5(int d1, int d2, double d3, double d4, int d5)
			: viewid(d1) {
			uvid = d2;
			uv[0] = d3;
			uv[1] = d4;
			ptid = d5;
		}
	};
	struct pt4 {
		int ptid;
		double xyz[3];
		// ���캯��
		pt4(int d1, double d2, double d3, double d4)
			: ptid(d1) {
			xyz[0] = d2;
			xyz[1] = d3;
			xyz[2] = d4;
		}
	};
	struct cal3 {
		double f, cx, cy;
		cal3(double d1, double d2, double d3)
			:f(d1) {
			cx = d2;
			cy = d3;
		}
	};
	struct cal8 {
		double fx, fy, cx, cy, k1, k2, p1, p2;
		cal8(double d1, double d2, double d3, double d4, double d5, double d6, double d7, double d8)
			:fx(d1), fy(d2), cx(d3), cy(d4), k1(d5), k2(d6), p1(d7), p2(d8) {
		}
	};
	struct track5 {
		int ptid;
		size_t viewid;
		int uvid;
		double uv[2];
		// ���캯��
		track5(int d1, size_t d2, int d3, double d4, double d5)
			: ptid(d1) {
			viewid = d2;
			uvid = d3;
			uv[0] = d4;
			uv[1] = d5;
		}
	};
	std::ptrdiff_t findViewIndex(const std::vector<eop8_>& views, int targetViewId);
	void groupAndSort(std::vector<track5>& track);
	void groupAndSort(vector<pt4>& pt);
	void add_gano2xyz(vector<pt4>& points, double mean, double stddev);
	void add_gano2eul(vector<eop8_>& eop, double mean, double stddev);
	void add_gano2c(vector<eop8_>& eop, double mean, double stddev);
	void add_gano2obs(vector<obs5>& obs, double mean, double stddev);
	void add_gano2cal(vector<cal3>& cal, double mean, double stddev);
	void add_gano2cal8(vector<cal8>& cal, double mean, double stddev);
	DataType tpe = colmap;
	string model;
	int nc;
};

