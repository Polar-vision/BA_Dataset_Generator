#include "myconvert.h"
#include <time.h>
#include <string.h>
using namespace std;
myconvert::myconvert(int)
{
}
myconvert::~myconvert(void)
{
}
string findparentpath()
{
	string originalPath(dt);
	size_t pos = originalPath.find_last_of("/\\");
	string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	return parentPath;
}
void myconvert::outlierdet()
{
	FILE *fp = fopen(dt, "r");
	fscanf(fp, "%d %d %d", &ncams, &n3Dpts, &n2Dprojs);
	fclose(fp);

	string parpth = findparentpath();
	string fea = parpth + "/" + "Feature.txt";
	string cam = parpth + "/" + "Cam.txt";
	string pt = parpth + "/" + "XYZ.txt";
	string cal = parpth + "/" + "cal.txt";
	string feanew = parpth + "/" + "Feature_new.txt";
	string camnew = parpth + "/" + "Cam_new.txt";
	string ptnew = parpth + "/" + "XYZ_new.txt";
	string calnew = parpth + "/" + "cal_new.txt";

	vector<int> deleteid;
	double* motf = (double*)malloc((ncams * 1) * sizeof(double));

	deleteid = readf(cal.c_str(), calnew.c_str(), motf, ncams);



	FILE* fpc = fopen(cam.c_str(), "r");
	FILE* fpc_new = fopen(camnew.c_str(), "w");
	double* motcam = (double*)malloc((ncams * 6) * sizeof(double));

	ba_readCameraPose(fpc, fpc_new, deleteid, motcam);
	//for (int i = 0; i < ncams; i++)
	//	printf("%f %f %f %f %f %f\n", motcam[6 * i], motcam[6 * i + 1], motcam[6 * i + 2], \
	//		motcam[6 * i + 3], motcam[6 * i + 4], motcam[6 * i + 5]);


	FILE* fpz = fopen(pt.c_str(), "r");
	double* motz = (double*)malloc((n3Dpts * 3) * sizeof(double));
	for (int i = 0; i < n3Dpts; i++)
	{
		fscanf(fpz, "%lf  %lf  %lf", motz + 3 * i, motz + 3 * i + 1, motz + 3 * i + 2);
		//printf("%f %f %f\n", motz[3 * i], motz[3 * i + 1], motz[3 * i + 2]);
	}


	FILE* fpt = fopen(fea.c_str(), "r");
	FILE* fpt_new = fopen(feanew.c_str(), "w");


	vector<int> deleteptid = readProjection(fpt, fpt_new, deleteid, motz, motcam);

	FILE* fpz_new = fopen(ptnew.c_str(), "w");
	for (int i = 0; i < n3Dpts; i++)
	{
		if (!deleteptid.empty())
		{
			auto it = std::find(deleteptid.begin(), deleteptid.end(), i);
			if (it != deleteptid.end())
				continue;
		}
		//printf("%f %f %f\n", motz[3 * i], motz[3 * i + 1], motz[3 * i + 2]);
		fprintf(fpz_new, "%f %f %f\n", motz[3 * i], motz[3 * i + 1], motz[3 * i + 2]);
	}
	fclose(fpz);
	fclose(fpz_new);

	free(motz);
	free(motcam);
	free(motf);
}

int myconvert::colmap2pba(string p1r,string p2r, string p3r,\
	string p1w, string p2w, string p3w, string p4w,\
	double* sigobs, double *sigeul, double *sigpos, double *sigxyz, double *sigcal){
	FILE* fp1 = NULL, * fp2 = NULL, * fp3 = NULL;

	int line = 0;
	vector<cal3> iop;//f, cx, cy
	vector<cal8> iop8;//fx fy cx cy k1 k2 p1 p2
	vector<pt4> points;//ptid, x, y, z
	//vector<eop8> eop;
	vector<eop8_> eop;
	vector<obs5> obs;
	vector<int> camid_group, cmidG, cmieG;
	vector<track5> track;//ptid, viewid, uvid, uv
	vector<int> uvN;
	vector<int> nview;

	std::vector<cv::Mat> K;
	std::vector<cv::Mat> distCoeffs;
	//cv::Mat K;

	//cv::Mat distCoeffs;
	std::vector<cv::Point2d> distortedPoints;
	std::vector<cv::Point2d> undistortedPixelPoints;
	printf("colmap2pba: Version 1.0\n");

	fp1 = fopen(p1r.c_str(), "r");
	fp2 = fopen(p2r.c_str(), "r");
	fp3 = fopen(p3r.c_str(), "r");

	if (!fp1 || !fp2 || !fp3){
		fprintf(stderr, "colmap2pba: Missing data file! \n");
		exit(1);
	}

	line = 0;
	vector<char> buffer3(1024);
	while (!feof(fp3)) {
		line++;
		if (line < 4) {
			if (line == 3) {
				assert(buffer3.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
				if (fgets(buffer3.data(), sizeof(char) * static_cast<int>(buffer3.size()), fp3) != nullptr) {
					char* pBuffer = NULL;
					pBuffer = buffer3.data();
					char tmp1[20], tmp2[20], tmp3[20], tmp4[20];
					int n = sscanf(pBuffer, "%s %s %s %s %d", &tmp1, &tmp2, &tmp3, &tmp4, &nc);
					if (n == 5) {
						printf("%s %s %s %s %d\n", tmp1, tmp2, tmp3, tmp4, nc);
					}
					else {
						printf("%s\n", "Error happened in the read operation!");
						break;
					}
				}
				else {
					break;
				}
			}
			else {
				assert(buffer3.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
				if (fgets(buffer3.data(), sizeof(char) * static_cast<int>(buffer3.size()), fp3) != nullptr)
					continue;
				else
					break;
			}
		}
		else {
			char* pBuffer = NULL;
			assert(buffer3.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
			if (fgets(buffer3.data(), sizeof(char) * static_cast<int>(buffer3.size()), fp3) != nullptr)
				pBuffer = buffer3.data();
			else
				break;

			int cameraid, width, height;
			double f, fx, fy, cx, cy, k1, k2, p1, p2;
			int pos;
			char modelBuf[100]; 
			int n = sscanf(pBuffer, "%d %s %n", &cameraid, &modelBuf, &pos);
			model = modelBuf;
			cmidG.push_back(cameraid);
			//printf("%d %s\n", cameraid, model.c_str());
			pBuffer += pos;
			if (model == "SIMPLE_PINHOLE") {
				n = sscanf(pBuffer, "%d %d %lf %lf %lf %n", \
					& width, &height, &f, &cx, &cy, &pos);
				iop.push_back(cal3(f, cx, cy));
				//printf("%d %d %lf %d %d\n", width, height, f, cx, cy);
			}if (int t = model.compare("OPENCV") == 0) {
				n = sscanf(pBuffer, "%d %d %lf %lf %lf %lf %lf %lf %lf %lf %n", \
					& width, &height, &fx, &fy, &cx, &cy, &k1, &k2, &p1, &p2, &pos);
				iop8.push_back(cal8(fx, fy, cx, cy, k1, k2, p1, p2));
				cv::Mat tmp_K = (cv::Mat_<double>(3, 3) <<
					fx, 0, cx,
					0, fy, cy,
					0, 0, 1);
				K.push_back(tmp_K);
				cv::Mat tmp_distCoeffs = (cv::Mat_<double>(1, 4) <<
					k1, k2, p1, p2);
				distCoeffs.push_back(tmp_distCoeffs);
				//printf("%d %d %lf %lf %d %d %lf %lf %lf %lf\n", width, height, fx, fy, cx, cy, k1, k2, p1, p2);
			}
		}
	}
	fclose(fp3);

	for (int i = 0; i < cmidG.size(); i++) {
		cmieG.push_back(cmidG[i] - i);
	}

	double t1 = clock();
	vector<char> buffer(1024*5000);
	line = 0;
	int mgh = 0;
	while (!feof(fp2)){
		line++;
		if (line < 5) {
			assert(buffer.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
			if (fgets(buffer.data(), sizeof(char) * static_cast<int>(buffer.size()), fp2) != nullptr)
				continue;
			else
				break;
		}else{
			char* pBuffer = NULL;
			assert(buffer.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
			if (fgets(buffer.data(), sizeof(char) * static_cast<int>(buffer.size()), fp2) != nullptr)
				pBuffer = buffer.data();
			else
				break;

			int viewid, cameraid;
			double w, x, y, z, t1, t2, t3;
			int n = sscanf(pBuffer, "%d %lf %lf %lf %lf %lf %lf %lf %d\n", \
				& viewid, &w, &x, &y, &z, &t1, &t2, &t3, &cameraid);

			for (int i = 0; i < cmidG.size(); i++) {
				if (cameraid == cmidG[i]) {
					cameraid -= (cmieG[i] - 1);
					break;
				}
			}

			double* ptr = qt2euc(w, x, y, z, t1, t2, t3);
			//printf("%f %f %f %f %f %f\n", ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5]);
			eop.push_back(eop8_(viewid, ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5],cameraid));
			assert(buffer.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
			if (fgets(buffer.data(), sizeof(char) * static_cast<int>(buffer.size()), fp2) != nullptr)
				pBuffer = buffer.data();
			else
				break;

			double u, v;
			int ptid;
			int uvid = 0;
			int pos;
			if (model == "SIMPLE_PINHOLE") {
				while (sscanf(pBuffer, "%lf %lf %d %n", &u, &v, &ptid, &pos) >= 3) {
					obs.push_back(obs5(viewid, uvid, u, v, ptid));
					uvid++;
					pBuffer += pos;
				}
			}if (model == "OPENCV") {
				while (sscanf(pBuffer, "%lf %lf %d %n", &u, &v, &ptid, &pos) >= 3) {
					distortedPoints.emplace_back(u, v);
					camid_group.push_back(cameraid);
					obs.push_back(obs5(viewid, uvid, u, v, ptid));
					uvid++;
					pBuffer += pos;
				}
			}
			mgh += uvid;
			uvN.push_back(uvid);
		}
	}
	if (model == "OPENCV") {
		undist(K, distCoeffs, distortedPoints, undistortedPixelPoints, camid_group);

		//for (size_t i = 0; i < distortedPoints.size(); ++i) {
		//	std::cout << "Original (distorted) point: " << distortedPoints[i]
		//		<< " -> Undistorted point: " << undistortedPixelPoints[i] << std::endl;
		//} 
		int i = 0;
		for (auto& val : obs) {
			//printf("%f %f  ", val.uv[0], val.uv[1]);
			val.uv[0] = undistortedPixelPoints[i].x;
			val.uv[1] = undistortedPixelPoints[i].y;
			//printf("%f %f\n", val.uv[0], val.uv[1]);
			i++;
		}
	}

	add_gano2obs(obs, sigobs[0], sigobs[1]);

	fclose(fp2);
	double t2 = clock();
	double dt1 = (t2 - t1) * 0.001;
	printf("time1:%f\n", dt1);

	double t3 = clock();
	vector<char> buffer1(1024*20);
	line = 0;
	
	while (!feof(fp1)){
		line++;
		if (line < 4) {
			assert(buffer1.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
			if (fgets(buffer1.data(), sizeof(char) * static_cast<int>(buffer1.size()), fp1) != nullptr)
				continue;
			else
				break;
		}else{
			char* pBuffer1 = NULL;
			assert(buffer1.size() <= static_cast<size_t>(std::numeric_limits<int>::max()));
			if (fgets(buffer1.data(), sizeof(char) * static_cast<int>(buffer1.size()), fp1) != nullptr)
				pBuffer1 = buffer1.data();
			else
				break;
			
			int ptid;
			double x, y, z, r, g, b, error;
			int pos;
			if (sscanf(pBuffer1, "%d %lf %lf %lf %lf %lf %lf %lf %n", \
				& ptid, &x, &y, &z, &r, &g, &b, &error, &pos) >= 8)
				pBuffer1 += pos;
			points.push_back(pt4(ptid, x, y, z));

			//if (ptid == 63) {
			//	printf("%s\n", buffer1.data());
			//}

			int viewid, uvid;
			int nv = 0;
			int ds = 3;
			while (sscanf(pBuffer1, "%d %d %n", &viewid, &uvid, &pos) >= 2) {
				//ds = sscanf(pBuffer1, "%d %d %n", &viewid, &uvid, &pos);
				//if (ptid == 63) {
				//	printf("%d %d %d %d\n", viewid, uvid, pos, ds);
				//	if (ds < 2) {
				//		printf("%s\n", pBuffer1);
				//		printf("%d %d %d %d\n", viewid, uvid, pos, ds);
				//	}
				//}
				nv++;
				pBuffer1 += pos;

				for (int i = 0, j = 0, sum; i < obs.size(); sum = uvN[j],i += sum, j++) {
					if (obs[i].viewid == viewid){
						obs5 ob = obs[i + uvid];
						if (ob.uvid == uvid && ob.ptid == ptid) {
		
							size_t viewidNew = findViewIndex(eop, viewid);
							if (eop[viewidNew].viewid == viewid) {
								track.push_back(track5(ptid, viewidNew, uvid, obs[i + uvid].uv[0], obs[i + uvid].uv[1]));
								break;
							}
							else {
								printf("error!");
							}
						}

					}
				}
			}
			nview.push_back(nv);
		}
	}
	
	fclose(fp1);
	double t4 = clock();
	double dt2 = (t4 - t3) * 0.001;
	printf("time2:%f\n", dt2);

	FILE* fp1w = NULL, * fp2w = NULL, * fp3w = NULL, * fp4w = NULL;
	//printf("s\n", "sssssssssssssssssssssssssss");
	groupAndSort(points);
	add_gano2xyz(points, sigxyz[0], sigxyz[1]);
	double t5 = clock();
	fp1w = fopen(p1w.c_str(), "w");
	for (const auto& it : points) {

		fprintf(fp1w, "%f %f %f\n", it.xyz[0], it.xyz[1], it.xyz[2]);
	}
	fclose(fp1w);

	add_gano2eul(eop, deg2rad(sigeul[0]), deg2rad(sigeul[1]));
	add_gano2c(eop, sigpos[0], sigpos[1]);
	fp2w = fopen(p2w.c_str(), "w");
	for (const auto& it : eop) {
		fprintf(fp2w, "%f %f %f %f %f %f %d\n", \
			it.eop[0], it.eop[1], it.eop[2], \
			it.eop[3], it.eop[4], it.eop[5],\
		    it.camid);
	}  
	fclose(fp2w);

	fp3w = fopen(p3w.c_str(), "w");
	if (model == "SIMPLE_PINHOLE") {
		add_gano2cal(iop, sigcal[0], sigcal[1]);
		for (int i = 0; i < nc; i++) {
			fprintf(fp3w, "%f %d %f\n", iop[i].f, 0, iop[i].cx);
			fprintf(fp3w, "%d %f %f\n", 0, iop[i].f, iop[i].cy);
			fprintf(fp3w, "%d %d %d\n", 0, 0, 1);
		}
	}
	if (model == "OPENCV") {
		add_gano2cal8(iop8, sigcal[0], sigcal[1]);
		for (int i = 0; i < nc; i++) {
			fprintf(fp3w, "%f %d %f\n", iop8[i].fx, 0, iop8[i].cx);
			fprintf(fp3w, "%d %f %f\n", 0, iop8[i].fy, iop8[i].cy);
			fprintf(fp3w, "%d %d %d\n", 0, 0, 1);
		}
	}
	fclose(fp3w);

	groupAndSort(track);
	nview.clear();
	int last_ptid = -1;
	int nv = 0;
	for (int i = 0; i < track.size(); i++) {
		if (track[i].ptid != last_ptid) {
			if (i != 0) {
				nview.push_back(nv);
				nv = 0;
			}
			last_ptid = track[i].ptid;
		}
		nv++;
		if (i == track.size() - 1) {
			nview.push_back(nv);
		}
	}
	fp4w = fopen(p4w.c_str(), "w");
	last_ptid = -1;
	for (int i = 0, j = 0, sum = 0; i < track.size(); i++) {
		if (track[i].ptid != last_ptid) {
			if (i == 0) 
				fprintf(fp4w, "%d", nview[j]);
			else
				fprintf(fp4w, "\n%d", nview[j]);
			j++;
			last_ptid = track[i].ptid;
		}
		fprintf(fp4w, "  %zu %f %f", track[i].viewid, track[i].uv[0], track[i].uv[1]);
	}
	fclose(fp4w);
	double t6 = clock();
	double dt3 = (t6 - t5) * 0.001;
	printf("time3:%f\n", dt3);
	return 0;
}
void myconvert::bal2pba()
{
	FILE* fp1 = NULL, * fp2 = NULL, * fp3 = NULL, * fp4 = NULL;
	string parpth = findparentpath();
	string newPath1 = parpth + "/" + "Feature.txt";
	string newPath2 = parpth + "/" + "Cam.txt";
	string newPath3 = parpth + "/" + "XYZ.txt"; 
	string newPath4 = parpth + "/" + "cal.txt";
	printf("bal2pba: Version 1.0\n");
	FILE* fp;
	int line = 0, n = 0, nframes = 0, ptno = 0, last_ptno = 0, sum_nframes = 0, sum_nframeP = 0, sum_n3DptsP = 0, k = 0;
	double* motstructs = NULL, * imgpts = NULL;
	int* frameno = NULL;
	double* ptr1 = NULL, * ptr2 = NULL, * ptr5 = NULL, * ptr6 = NULL, * ptr7 = NULL;
	int* ptr3 = NULL, * ptr4 = NULL;
	fp = fopen(dt, "r");
	double nx = 0, ny = 0, nz = 0, t1 = 0, t2 = 0, t3 = 0;
	Eigen::Matrix3d MatInv;
	double R[9];
	if (fp == NULL)
	{
		fprintf(stderr, "DataPreprocess: Missing data file! \n");
		exit(1);
	}
	else
	{
		while (!feof(fp))
		{
			line++;
			if (line == 1)
			{
				int num = fscanf(fp, "%d %d %d", &ncams, &n3Dpts, &n2Dprojs);
				motstructs = (double*)malloc((ncams * 9 + n3Dpts * 3) * sizeof(double));//���ڴ洢R t f k1 k2 X
				if (motstructs == NULL)
				{
					fprintf(stderr, "bal2pba error: Memory allocation for 'motstruct' failed \n");
					exit(1);
				}

				imgpts = (double*)malloc(n2Dprojs * 2 * sizeof(double));//���ڴ洢����ͶӰ���2D����
				frameno = (int*)malloc(n2Dprojs * sizeof(int));
				if (imgpts == NULL)
				{
					fprintf(stderr, "bal2pba error: Memory allocation for 'imgpts' failed\n");
					exit(1);
				}
				ptr1 = imgpts;
				ptr3 = frameno;
				ptr5 = motstructs;
				ptr7 = motstructs;
				fp1 = fopen(newPath1.c_str(), "w");
				fp2 = fopen(newPath2.c_str(), "w");
				fp3 = fopen(newPath3.c_str(), "w");
				fp4 = fopen(newPath4.c_str(), "w");
				//fprintf(fp4, "%d %d %d\n", 1, 0, 0);
				//fprintf(fp4, "%d %d %d\n", 0, 1, 0);
				//fprintf(fp4, "%d %d %d\n", 0, 0, 1);
			}
			else
			{
				if (sum_nframes < n2Dprojs)
				{

					if ((line - 1) > n2Dprojs)
					{
						ptno += 1;
					}
					else
					{
						n = readNInts(fp, &(*ptr3++), 1);  //��ͼID
						if (n != 1)
							break;
						n = readNInts(fp, &ptno, 1); //3D��ID
						if (n != 1)
							break;
						n = readNDoubles(fp, ptr1, 2);
						if (n != 2)
						{
							fprintf(stderr, "bal2pba:reading image projections wrong!\n");
							break;
						}
						ptr1 += 2;
					}

					if (ptno > last_ptno)
					{
						ptr2 = imgpts + sum_nframes * 2;
						ptr4 = frameno + sum_nframes;
						fprintf(fp1, "%d ", nframes);
						for (int i = 0; i < nframes; i++)
						{
							fprintf(fp1, "%d %0.5f %0.5f ", *ptr4++, *ptr2, *(ptr2 + 1));
							ptr2 += 2;
						}
						fprintf(fp1, "\n");
						last_ptno = ptno;
						sum_nframes += nframes;
						nframes = 1;
					}
					else
					{
						nframes++;
					}
				}
				else if (sum_nframeP < ncams * 9)
				{
					sum_nframeP++;
					n = readNDoubles(fp, ptr5, 1);
					if (n != 1)
					{
						fprintf(stderr, "bal2pba:reading image parameters wrong!\n");
						break;
					}
					//fprintf(fp2, "%f ", *ptr5);

					if ((sum_nframeP - 1) % 9 == 0)//���1
						nx = *ptr5;
					else if ((sum_nframeP - 2) % 9 == 0)//���2
						ny = *ptr5;
					else if ((sum_nframeP - 3) % 9 == 0)//���3
					{
						nz = *ptr5;
						double theta = sqrt(nx * nx + ny * ny + nz * nz);//theta��
						double kx = nx / theta;//��λ����n
						double ky = ny / theta;
						double kz = nz / theta;

						//���תR����1
						//double K[9];
						//K[0] = 0;    K[1] = -kz;   K[2] = ky;
						//K[3] = kz;   K[4] = 0;     K[5] = -kx;
						//K[6] = -ky;  K[7] = kx;    K[8] = 0;
						//double K2[9];
						//K2[0] = K[0] * K[0] + K[1] * K[3] + K[2] * K[6];
						//K2[1] = K[0] * K[1] + K[1] * K[4] + K[2] * K[7];
						//K2[2] = K[0] * K[2] + K[1] * K[5] + K[2] * K[8];

						//K2[3] = K[3] * K[0] + K[4] * K[3] + K[5] * K[6];
						//K2[4] = K[3] * K[1] + K[4] * K[4] + K[5] * K[7];
						//K2[5] = K[3] * K[2] + K[4] * K[5] + K[5] * K[8];

						//K2[6] = K[6] * K[0] + K[7] * K[3] + K[8] * K[6];
						//K2[7] = K[6] * K[1] + K[7] * K[4] + K[8] * K[7];
						//K2[8] = K[6] * K[2] + K[7] * K[5] + K[8] * K[8];
						//double I[9];
						//I[0] = 1; I[1] = 0; I[2] = 0;
						//I[3] = 0; I[4] = 1; I[5] = 0;
						//I[6] = 0; I[7] = 0; I[8] = 1;

						//R[0] = I[0] + sin(theta) * K[0] + (1 - cos(theta)) * K2[0];
						//R[1] = I[1] + sin(theta) * K[1] + (1 - cos(theta)) * K2[1];
						//R[2] = I[2] + sin(theta) * K[2] + (1 - cos(theta)) * K2[2];
						//R[3] = I[3] + sin(theta) * K[3] + (1 - cos(theta)) * K2[3];
						//R[4] = I[4] + sin(theta) * K[4] + (1 - cos(theta)) * K2[4];
						//R[5] = I[5] + sin(theta) * K[5] + (1 - cos(theta)) * K2[5];
						//R[6] = I[6] + sin(theta) * K[6] + (1 - cos(theta)) * K2[6];
						//R[7] = I[7] + sin(theta) * K[7] + (1 - cos(theta)) * K2[7];
						//R[8] = I[8] + sin(theta) * K[8] + (1 - cos(theta)) * K2[8];

						//����
						//const Vector3d phi(nx, ny, nz);
						//printf("%f %f %f\n", nx, ny, nz);
						//printf("%f %f %f\n", phi(0), phi(1), phi(2));
						// Matrix3d R_test = ExpSO3(phi);
						//Vector3d eu = LogSO3(R_test);
						//printf("%f %f %f\n", eu(0), eu(1), eu(2));
						//printf("%f %f %f\n", R_test(0, 0), R_test(0, 1), R_test(0, 2));
						//printf("%f %f %f\n", R_test(1, 0), R_test(1, 1), R_test(1, 2));
						//printf("%f %f %f\n", R_test(2, 0), R_test(2, 1), R_test(2, 2));

						//���תR����2
						double q[4];//��Ԫ��

						q[0] = cos(theta / 2);
						q[1] = sin(theta / 2) * kx;
						q[2] = sin(theta / 2) * ky;
						q[3] = sin(theta / 2) * kz;
						//double angle_axis[3];
						//QuaternionToAngleAxis(q, angle_axis);
						double a = q[0], b = q[1], c = q[2], d = q[3];
						//����ϵ
						R[0] = 1 - 2 * c * c - 2 * d * d;
						R[1] = 2 * b * c - 2 * a * d;
						R[2] = 2 * b * d + 2 * a * c;
						R[3] = 2 * b * c + 2 * a * d;
						R[4] = 1 - 2 * b * b - 2 * d * d;
						R[5] = 2 * c * d - 2 * a * b;
						R[6] = 2 * b * d - 2 * a * c;
						R[7] = 2 * c * d + 2 * a * b;
						R[8] = 1 - 2 * b * b - 2 * c * c;
						//printf("%f %f %f\n", R[0], R[1], R[2]);
						//printf("%f %f %f\n", R[3], R[4], R[5]);
						//printf("%f %f %f\n", R[6], R[7], R[8]);
						//double quaternion[4];
						//RotationMatrixToQuaternion(R, quaternion);
						double eulerAngles[3];
						rotationMatrixToEulerAngles(R, eulerAngles);//Rתŷ���ǣ�kappa��phi��omega��
						//double newR[9];
						//eulerAnglesToRotationMatrix(eulerAngles, newR);
						//Eigen::Matrix3d matV(R);//������
						//printf("%f %f %f\n", matV(0, 0), matV(0, 1), matV(0, 2));
						//printf("%f %f %f\n", matV(1, 0), matV(1, 1), matV(1, 2));
						//printf("%f %f %f\n", matV(2, 0), matV(2, 1), matV(2, 2));
						//MatInv = matV.inverse();
						ptr7 = motstructs + 9 * k;
						*(ptr7 + 0) = eulerAngles[0];//��ŷ����д��motstructs
						*(ptr7 + 1) = eulerAngles[1];
						*(ptr7 + 2) = eulerAngles[2];
						//*(ptr7 + 0) = nx;
						//*(ptr7 + 1) = ny;
						//*(ptr7 + 2) = nz;
					}
					else if ((sum_nframeP - 4) % 9 == 0)//t1
						t1 = *ptr5;
					else if ((sum_nframeP - 5) % 9 == 0)//t2
						t2 = *ptr5;
					else if ((sum_nframeP - 6) % 9 == 0)//t3
					{
						t3 = *ptr5;
						double center[3];
						center[0] = -(R[0] * t1 + R[3] * t2 + R[6] * t3);//tתΪc
						center[1] = -(R[1] * t1 + R[4] * t2 + R[7] * t3);
						center[2] = -(R[2] * t1 + R[5] * t2 + R[8] * t3);

						*(ptr7 + 3) = center[0];//��cд��motstructs
						*(ptr7 + 4) = center[1];
						*(ptr7 + 5) = center[2];
						//*(ptr7 + 3) = t1;//��cд��motstructs
						//*(ptr7 + 4) = t2;
						//*(ptr7 + 5) = t3;
						//��Դ�ļ�����Ǻ�tת��Ϊŷ���Ǻ�c��д������ļ�
						fprintf(fp2, "%f %f %f %f %f %f ", *ptr7, *(ptr7 + 1), *(ptr7 + 2), *(ptr7 + 3), *(ptr7 + 4), *(ptr7 + 5));

						/*�����ǲ��Դ��룬���ڲ������תR��tתc����ȷ�ԣ����ò�ͬ��ʽ���и���任�õ���ͬ�������ת����ȷ*/
												//double angle_axis[3];
												//angle_axis[0] = nx;
												//angle_axis[1] = ny;
												//angle_axis[2] = nz;
												//double pt[3];//3D��
												//  
												//pt[0] = -0.612000;
												//pt[1] = 0.571759;
												//pt[2] = -1.847081;
												//double result[3];
												//AngleAxisRotatePoint(angle_axis, pt, result);//����Rodrigues' formula��pt���и���任
												//double P1[3];
												//P1[0] = result[0] + t1;
												//P1[1] = result[1] + t2;
												//P1[2] = result[2] + t3;
												//double P2[3];//����R��c��pt���и���任
												//P2[0] = R[0] * (pt[0] - center[0]) + R[1] * (pt[1] - center[1]) + R[2] * (pt[2] - center[2]);
												//P2[1] = R[3] * (pt[0] - center[0]) + R[4] * (pt[1] - center[1]) + R[5] * (pt[2] - center[2]);
												//P2[2] = R[6] * (pt[0] - center[0]) + R[7] * (pt[1] - center[1]) + R[8] * (pt[2] - center[2]);
												//double P3[3];//����R��t��pt���и���任
												//P3[0] = R[0] * pt[0] + R[1] * pt[1] + R[2] * pt[2] + t1;
												//P3[1] = R[3] * pt[0] + R[4] * pt[1] + R[5] * pt[2] + t2;
												//P3[2] = R[6] * pt[0] + R[7] * pt[1] + R[8] * pt[2] + t3;
												//int iii = 1;
					}
					else if ((sum_nframeP - 7) % 9 == 0)//f
						fprintf(fp4, "%f\n", *ptr5);
					//else
						//fprintf(fp2, "%f ", *ptr5);

					if (sum_nframeP % 9 == 0)
					{
						k++;
						fprintf(fp2, "\n");
					}

					ptr5 += 1;
				}
				else
				{
					sum_n3DptsP++;
					n = readNDoubles(fp, ptr5, 1);
					if (n != 1)
					{
						fprintf(stderr, "DataPreprocess:reading 3D points wrong!\n");
						break;
					}
					fprintf(fp3, "%f ", *ptr5);
					if (sum_n3DptsP % 3 == 0)
						fprintf(fp3, "\n");
					ptr5 += 1;
				}
			}

		}
		//printf("%d\n", sum_nframes);
		fclose(fp);
		fclose(fp1);
		fclose(fp2);
		fclose(fp3);
		fclose(fp4);
		free(imgpts);
		free(motstructs);
		free(frameno);
	}
}