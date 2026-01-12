#include "triangulation.h"
#define _CRT_SECURE_NO_WARNINGS
void triangulation::ba_saveTriangulatedxyz(const char* sz3Dpt, double* p)
{
	static int i = 0;
	double x, y, z;
	FILE* fp = NULL;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		errno_t err = fopen_s(&fp,sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			x = *(p + m_ncams * 6 + i * 3);
			y = *(p + m_ncams * 6 + i * 3 + 1);
			z = *(p + m_ncams * 6 + i * 3 + 2);
			//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1, x, y, z);
			fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", x, y, z);
			//fprintf(fp, "%0.5lf\n", parallax * 180 / PI);
		}
		fclose(fp);
	}
	free(p);
}
void triangulation::ba_readCameraPose(FILE* fp, double* params)
{
	int n, num, lineno = 0;
	double* tofilter = NULL;
	double* pPrams = params;
	int* pm_C = m_C;
	//the number of element per line is 8, it represents that focal length vary, or it is constant
	num = countNDoubles(fp);
	if (num == 8)
	{
		m_bFocal = true;
		m_K = (double*)malloc(m_ncams * 2 * sizeof(double));
		tofilter = (double*)malloc(8 * sizeof(double));
	}
	else if (num == 7) {
		tofilter = (double*)malloc(7 * sizeof(double));
	}
	else if (num == 6) {
		tofilter = (double*)malloc(6 * sizeof(double));
	}
		

	while (!feof(fp))
	{
		if (num == 6)
		{
			n = readNDoubles(fp, tofilter, 6);
			if (n == -1)
				break;
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			*pm_C = 1;
		}
		if (num == 7)
		{
			n = readNDoubles(fp, tofilter, 7);
			if (n == -1)
				break;
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];
			*pm_C = static_cast<int>(tofilter[6]);
		}
		if (num == 8)
		{
			n = readNDoubles(fp, tofilter, 8);
			if (n == -1)
				break;
			pPrams[0] = tofilter[0];	pPrams[1] = tofilter[1];	pPrams[2] = tofilter[2];
			pPrams[3] = tofilter[3];	pPrams[4] = tofilter[4];	pPrams[5] = tofilter[5];

			m_K[lineno * 2] = tofilter[6];
			m_K[lineno * 2 + 1] = tofilter[7];
		}

		pPrams += 6;
		if (num == 6 || num == 7)
			pm_C += 1;
		++lineno;
	}
	if (tofilter != NULL) {
		free(tofilter);
		tofilter = NULL;
	}
}
triangulation::triangulation(int mag)
{
}
triangulation::~triangulation(void)
{
}
bool triangulation::ba_initialize(char* szCamera, char* szFeature, char* szCalib, char* szXYZ)
{
	printf("triangulation: Version 1.0\n");
	FILE* fp;
	
	//must input initial initial camera pose file and projection image points file
	errno_t err = fopen_s(&fp,szCamera, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "triangulation: Missing initial camera poses file! \n");
		exit(1);
	}
	else
		fclose(fp);

	err = fopen_s(&fp,szFeature, "r");
	if (fp == NULL)
	{
		fprintf(stderr, "triangulation: Missing feature projection points file! \n");
		exit(1);
	}
	else
		fclose(fp);

	if (tpe == colmap) {
		if (szCalib != NULL)
		{
			m_bFocal = false;
			m_K = (double*)malloc(9 * sizeof(double) * nc);
			ba_readCameraPoseration(szCalib, m_K);
		}
	}
	else if (tpe == bal) {

	}
	
	if (szXYZ != NULL)
		m_bProvideXYZ = true;
	zu = 2;
	savePara = 0;
	//read camera pose & features images projs, and initialize features points( three kinds of angle )
	pba_readAndInitialize(szCamera, szFeature, szCalib, &m_ncams, &m_n3Dpts, &m_n2Dprojs, &m_motstruct,//number of camera, 3D points, 2D projection points,6 camera pose and 3 feature parameters
		&m_imgpts, &m_archor, &m_vmask, &m_umask, &m_photo, &m_feature, &m_archorSort);

	//for (int i = m_n3Dpts-1; i < m_n3Dpts; i++) {
	//	double* ptr = m_motstruct + 6 * m_ncams + 3 * i;
	//	printf("%f %f %f\n", *ptr, *(ptr + 1), *(ptr + 2));
	//}
	//double* ptr = m_motstruct + 6 * m_ncams + 3 * 67425;
	//printf("%f %f %f\n", *ptr, *(ptr + 1), *(ptr + 2));
	//double* ptr = m_motstruct;
	//printf("%f %f %f %f %f %f\n", *ptr, *(ptr + 1), *(ptr + 2), *(ptr + 3), *(ptr + 4), *(ptr + 5));

	string originalPath(szCamera);
	size_t pos = originalPath.find_last_of("/\\");
	string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	string newPath1 = parentPath + "/" + "InitialParallax.txt";
	string newPath2 = parentPath + "/" + "XYZ1.txt";
	string newPath3 = parentPath + "/" + "XYZ1.txt";

	if (savePara == 1)
		pba_saveInitialParallax(newPath1.c_str(), m_motstruct);
	if (zu == 1)
		ba_saveTriangulatedxyz(newPath2.c_str(), m_motstruct);
	else if (zu == 2) {
		// printf("entered\n");
		pba_angle2xyz(m_motstruct);
		pba_saveInitialXYZ(newPath3.c_str(), m_motstruct);//ת��Ϊxyz
	}
	else
		;//Other codes

	printf("Number of cameras: %d\n", m_ncams);
	printf("Number of points: %d\n", m_n3Dpts);
	printf("Number of projections: %d\n", m_n2Dprojs);
	return true;
}
void triangulation::pba_angle2xyz(double* p)
{
	//double* n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));
	static int i, j;
	// double* pAngle;
	double xj[3], xk[3];
	double Tik[3];
	int nM, nA;
	double Dik;
	double w, w2;
	double az, el;
	double* ptr = NULL, *ptr1 = NULL, *ptr4 = NULL;
	double* ptr2 = NULL, * ptr3 = NULL;
	int cnp = 6, pnp = 3;
	for (i = 0; i < m_n3Dpts; i++)
	{
		ptr4 = p + m_ncams * 6 + i * 3;
		//pAngle = ptr;//azimuth angle
		az = *ptr4;
		el = *(ptr4 + 1);
		w = *(ptr4 + 2);//parallax angle

		xj[0] = sin(az) * cos(el);//x
		xj[1] = sin(el);//y
		xj[2] = cos(az) * cos(el);//z

		nM = *(m_archor + i * 3 + 1);
		nA = *(m_archor + i * 3 + 2);
		
		ptr2 = p + nM * 6 + 3;
		ptr3 = p + nA * 6 + 3;


		Tik[0] = *ptr3 - *ptr2;
		Tik[1] = *(ptr3 + 1) - *(ptr2 + 1);
		Tik[2] = *(ptr3 + 2) - *(ptr2 + 2);
		

		Dik = sqrt(Tik[0] * Tik[0] + Tik[1] * Tik[1] + Tik[2] * Tik[2]);
		

		w2 = acos((xj[0] * Tik[0] + xj[1] * Tik[1] + xj[2] * Tik[2]) / Dik);


		xk[0] = (Dik * sin(w2 + w) * xj[0]) / sin(w);
		xk[1] = (Dik * sin(w2 + w) * xj[1]) / sin(w);
		xk[2] = (Dik * sin(w2 + w) * xj[2]) / sin(w);


		ptr = p + 6 * m_ncams + i * 3;
		ptr1 = p + nM * 6 + 3;

		if (ptr != nullptr) {
			*ptr = *ptr1 + xk[0];
			*(ptr + 1) = *(ptr1 + 1) + xk[1];
			*(ptr + 2) = *(ptr1 + 2) + xk[2];

			//if (i == 384) {
			//	printf("%f %f %f\n",\
			//		*ptr, \
			//		*(ptr + 1),\
			//		*(ptr + 2));
			//	printf("%d %d\n", nM, nA);
			//	printf("%f %f %f\n", *ptr2, *(ptr2 + 1), *(ptr2 + 2));
			//	printf("%f %f %f\n", *ptr3, *(ptr3 + 1), *(ptr3 + 2));
			//	printf("%f %f %f\n", Tik[0], Tik[1], Tik[2]);
			//	printf("%f\n", Dik);
			//}
		}
		else {
			throw std::runtime_error("��ָ�����?: ptr ����Ϊ NULL");
		}
	}

	//return n3DptsT;
}
void triangulation::pba_saveInitialXYZ(const char* sz3Dpt, double* p)
{
	static int i = 0;
	double dx, dy, dz;
	FILE* fp = NULL;
	// int nM, nA;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		errno_t err = fopen_s(&fp,sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			dx = *(p + m_ncams * 6 + i * 3);
			dy = *(p + m_ncams * 6 + i * 3 + 1);
			dz = *(p + m_ncams * 6 + i * 3 + 2);
			//nM = *(m_archor + i * 3 + 1);
			//nA = *(m_archor + i * 3 + 2);
			//fprintf(fp, "%0.5lf %0.5lf %0.5lf %d %d\n", dx, dy, dz,nM, nA);
			//fprintf(fp, "%d %0.5lf %0.5lf %0.5lf\n", i+1,dx, dy, dz);
			fprintf(fp, "%0.5lf %0.5lf %0.5lf\n", dx, dy, dz);
		}
		fclose(fp);
	}
	free(p);
}
void triangulation::pba_saveInitialParallax(const char* sz3Dpt, double* p)
{
	//double* n3DptsT = (double*)malloc((m_n3Dpts * 3) * sizeof(double));

	static int i = 0;
	double azimuth, elevation, parallax;
	FILE* fp = NULL;
	double* pAngle;
	//save features xyz
	if (sz3Dpt != NULL)
	{
		errno_t err = fopen_s(&fp,sz3Dpt, "w");
		for (i = 0; i < m_n3Dpts; i++)
		{
			pAngle = p + m_ncams * 6 + i * 3;

			azimuth = *pAngle;
			elevation = *(pAngle + 1);
			parallax = *(pAngle + 2);
			fprintf(fp, "%0.5lf     %0.5lf     %0.5lf\n", azimuth, elevation, parallax);
			//fprintf(fp, "%0.5lf\n", parallax*180/PI);
		}
		fclose(fp);
	}
	free(p);
}
void triangulation::pba_readProjectionAndInitilizeFeature(FILE* fp,
	double* params, double* projs, char* vmask, int ncams,
	int* archor, char* umask, int* nphoto, int* nfeature, int* archorSort)
{
	int n;
	int nframes;
	int ptno = 0, cur;

	int nproj2D = 0;

	int count = 0;

	int frameno;
	int feastart = 0;

	int nP, nP2;

	double* ptr1 = projs;

	int i, j;
	int  sum, cnp = 6;

	int nFlag;

	int* ptr2;
	bool bAdjust;

	m_smask = (char*)malloc(m_ncams * m_ncams * sizeof(char));
	memset(m_smask, 0, m_ncams * m_ncams * sizeof(char));

	int* archorEx = new int[m_n3Dpts * 2];

	bool bM, bN;
	int max_nframes = -1;

// printf("%s\n","log1");
	//read all projection point, initialize three feature angle at the same time
	while (!feof(fp))
	{
		//if (ptno == 67425) {
		//	printf("%d\n", ptno);
		//}
		nFlag = 0;
		n = readNInts(fp, &nframes, 1);
		if (n != 1)
			break;
// printf("%d\n",nframes);
		archor[ptno * 3] = nframes;
		cur = 0;
		bM = bN = false;
		for (i = 0, sum = 0; i < nframes; ++i)
		{
			n = readNInts(fp, &frameno, 1);
			nphoto[nproj2D] = frameno;
			nfeature[nproj2D] = ptno;
			nproj2D++;

			if (frameno >= ncams)
			{
				fprintf(stderr, "ParallaxBA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles(fp, ptr1, 2);
			ptr1 += 2;
			if (n != 3)
			{
				fprintf(stderr, "ParallaxBA:reading image projections wrong!\n");
				return;
			}

			if (bM && bN)
			{
				ptr2 = archorSort + ptno * 2;
				//bAdjust = pba_initializeOtheArchors_Mindw( //_Mindw
				bAdjust = pba_initializeOtheArchors( //Maxdw
					projs + feastart * 2,
					nphoto + feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams * cnp + ptno * 3,
					ptr2,
					sum,
					i,
					ptno);
				if (bAdjust)
				{
					archor[ptno * 3 + 1] = *(nphoto + feastart + ptr2[0]);
					archor[ptno * 3 + 2] = *(nphoto + feastart + ptr2[1]);

					archorEx[ptno * 2] = ptr2[0];
					archorEx[ptno * 2 + 1] = ptr2[1];
				}
				sum++;
			}

			if (bM && !bN)
			{
				bool bLast = (i == nframes - 1);
				bool bT = pba_initializeAssoArchor(
					projs + feastart * 2,
					nphoto + feastart,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams * cnp + ptno * 3,
					0,
					1,
					ptno,
					bLast);

				if (bT)
				{
					archorSort[ptno * 2 + 1] = i;
					archor[ptno * 3 + 2] = nphoto[count];
					sum++;

					archorEx[ptno * 2 + 1] = i;

					bN = true;
				}
			}

			if (!bM)
			{
				// printf("%d\n",nframes);
				bool bLast = (i == nframes - 2);
				bool bT = pba_initializeMainArchor(
					projs + feastart * 2,
					m_motstruct,
					m_K,
					m_motstruct + m_ncams * cnp + ptno * 3,
					nphoto[count],
					ptno,
					m_KR);

				archorSort[ptno * 2] = i;
				archor[ptno * 3 + 1] = nphoto[count];
				sum++;

				archorEx[ptno * 2] = i;
				bM = true;
			}
			count++;
		}
// printf("%s\n","log2");
		//set masks for U and S matrix             
		for (i = 0; i < nframes; i++)
		{
			nP = nphoto[feastart + i];
			int nM_ = archor[ptno * 3 + 1];
			int nA_ = archor[ptno * 3 + 2];
			int tmp3 = archor[ptno * 3 + 3];

			if (nM_ < nP)
				umask[nM_ * (ncams)+nP] = 1;
			else
				umask[nP * (ncams)+nM_] = 1;


			if (nA_ < nP)
				umask[nA_ * (ncams)+nP] = 1;
			else
				umask[nP * (ncams)+nA_] = 1;

			umask[nP * ncams + nP] = 1;

			for (j = i; j < nframes; j++)
			{
				nP2 = nphoto[feastart + j];

				if (nP == nP2)
					m_smask[nP * m_ncams + nP2] = 1;
				else if (nP < nP2)
					m_smask[nP * m_ncams + nP2] = 1;
				else
				{
					m_smask[nP2 * m_ncams + nP] = 1;
				}

			}
		}
		feastart += nframes;
		ptno++;
	}
	//count number of non-zero element in S matrix
	m_nS = 0;
	for (i = 0; i < m_ncams; i++)
	{
		for (j = 0; j < m_ncams; j++)
		{
			if (m_smask[i * m_ncams + j] == 1)
			{
				m_nS++;
			}
		}
	}
}
void triangulation::pba_readAndInitialize(char* camsfname, char* ptsfname, char* calibfname, int* ncams,
	int* n3Dpts, int* n2Dprojs,
	double** motstruct, double** imgpts,
	int** archor, char** vmask,
	char** umask, int** nphoto,
	int** nfeature, int** archorSort)
{
	FILE* fpc, * fpp;
	// int i, tmp1, tmp2;
	// double ptMain[3], ptA[3];
	// double dW1, dW2;

	//calculate number of cameras, 3D points and projection points
	errno_t err = fopen_s(&fpc,camsfname, "r");
	*ncams = findNcameras(fpc);
	m_ncams = *ncams;
	m_C = (int*)malloc(sizeof(int) * m_ncams);
	if (tpe == bal) {
		if (calibfname != NULL)
		{
			m_bFocal = false;
			m_K = (double*)malloc(*ncams * 9 * sizeof(double));
			ba_readCameraPoseration(calibfname, m_K);
		}
	}
	else if (tpe == colmap) {

	}

	// printf("%s\n","log1");

	err = fopen_s(&fpp,ptsfname, "r");//��ƥ����?�����ļ�
	readNpointsAndNprojections(fpp, n3Dpts, 3, n2Dprojs, 2);//��ȡ3D������ĸ�����ͶӰ�����

	*motstruct = (double*)malloc((*ncams * 6 + *n3Dpts * 3) * sizeof(double));//���ڴ洢���е��ⷽλԪ�� �� �����������?3D����
	if (*motstruct == NULL)
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'motstruct' failed \n");
		exit(1);
	}

	*imgpts = (double*)malloc(*n2Dprojs * 2 * sizeof(double));//���ڴ洢����ͶӰ���?2D����
	if (*imgpts == NULL)
	{
		fprintf(stderr, "ParallaxBA error: Memory allocation for 'imgpts' failed\n");
		exit(1);
	}

	

	rewind(fpc);
	rewind(fpp);

	//allocate indicator of U
	*umask = (char*)malloc(*ncams * *ncams);
	memset(*umask, 0, *ncams * *ncams * sizeof(char));

	//allocate main and associate anchors
	*archor = (int*)malloc(*n3Dpts * 3 * sizeof(int));//
	memset(*archor, -1, *n3Dpts * 3 * sizeof(int)); //��������ê���λ�ó�ʼ���?-1

	*nphoto = (int*)malloc(*n2Dprojs * 3 * sizeof(int));//
	*nfeature = (int*)malloc(*n2Dprojs * 3 * sizeof(int));//
	*archorSort = (int*)malloc(*n3Dpts * 3 * sizeof(int));

	ba_readCameraPose(fpc, *motstruct);

	fclose(fpc);

	//Update KR
	m_KR = (double*)malloc(m_ncams * 9 * sizeof(double));//����ڲξ���? �� ��ת����ĳ˻�?
	m_KdA = (double*)malloc(m_ncams * 9 * sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	m_KdB = (double*)malloc(m_ncams * 9 * sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�
	m_KdG = (double*)malloc(m_ncams * 9 * sizeof(double));//��ת����M_x (��) M_y (��) M_z (��)��һ�׵�

	if (zu == 1)
	{
		m_P = (double*)malloc(m_ncams * 12 * sizeof(double));
		ba_constructP(m_P, m_K, *motstruct);
	}
	else
		ba_updateKR(m_KR, m_KdA, m_KdB, m_KdG, m_K, *motstruct);


	//test
	// for (int i = 0; i < m_ncams; i++)
	// {
	// 	printf("%f %f %f\n", m_KR[i * 9], m_KR[i * 9 + 1], m_KR[i * 9 + 2]);//, m_P[i * 12], m_P[i * 12 + 1], m_P[i * 12 + 2]);
	// 	printf("%f %f %f\n", m_KR[i * 9 + 3], m_KR[i * 9 + 4], m_KR[i * 9 + 5]);//, m_P[i * 12 + 4], m_P[i * 12 + 5], m_P[i * 12 + 6]);
	// 	printf("%f %f %f\n", m_KR[i * 9 + 6], m_KR[i * 9 + 7], m_KR[i * 9 + 8]);//, m_P[i * 12 + 8], m_P[i * 12 + 9], m_P[i * 12 + 10]);
	// }

	if (zu == 1)
		ba_readProjectionAndTriangulateFeature(fpp, *imgpts, *ncams);
	else
	{
		pba_readProjectionAndInitilizeFeature(fpp,
			*motstruct + *ncams * 6,
			*imgpts,
			*vmask,
			*ncams,
			*archor,
			*umask,
			*nphoto,
			*nfeature,
			*archorSort);
	}

	fclose(fpp);
}


bool triangulation::pba_initializeMainArchor(
	double* imgpts,
	double* camera,
	double* K,
	double* feature,
	int nP,
	int FID,
	double* KR)
{
	//solve  KRX = x
	Vector3d x;
	if (m_bProvideXYZ)
	{
		x(0) = m_XYZ[FID * 3] - *(camera + nP * 6 + 3);
		x(1) = m_XYZ[FID * 3 + 1] - *(camera + nP * 6 + 4);
		x(2) = m_XYZ[FID * 3 + 2] - *(camera + nP * 6 + 5);
	}
	else
	{
		double* ptr = m_KR + nP * 9;
		Matrix3d  A;
		A << ptr[0], ptr[1], ptr[2], ptr[3], ptr[4], ptr[5], ptr[6], ptr[7], ptr[8];
		
		double matx[3];
		matx[0] = imgpts[0];
		matx[1] = imgpts[1];
		matx[2] = 1;
		
		Vector3d  b(matx);
		
		x = A.colPivHouseholderQr().solve(b);
		// printf("%f %f %f\n",x(0),x(1),x(2));
	}

	double* pKR = KR + nP * 9;
	double t = pKR[6] * x(0) + pKR[7] * x(1) + pKR[8] * x(2);

	//compute azimuth and elevation angle
	double dDAngle = atan2(x(0), x(2));			
	double dHAngle = atan2(x(1), sqrt(x(0) * x(0) + x(2) * x(2)));

	feature[0] = dDAngle;
	feature[1] = dHAngle;
	feature[2] = 0;

	// printf("%f %f %f\n",feature[0],feature[1],feature[2]);

	// printf("%f\n",t);
	if (t < 0)
		return true;
	else
		return false;
}

bool triangulation::pba_initializeAssoArchor(
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int nMI,
	int nAI,
	int FID,
	bool bLast)
{
	int nM = photo[nMI];                              
	int nA = photo[nAI];                              

	Vector3d  xM, xA;

	if (m_bProvideXYZ)
	{
		xM[0] = m_XYZ[FID * 3] - *(camera + nM * 6 + 3);
		xM[1] = m_XYZ[FID * 3 + 1] - *(camera + nM * 6 + 4);
		xM[2] = m_XYZ[FID * 3 + 2] - *(camera + nM * 6 + 5);

		xA[0] = m_XYZ[FID * 3] - *(camera + nA * 6 + 3);
		xA[1] = m_XYZ[FID * 3 + 1] - *(camera + nA * 6 + 4);
		xA[2] = m_XYZ[FID * 3 + 2] - *(camera + nA * 6 + 5);
	}
	else
	{
		//Main anchor ray
		double* ptr1 = m_KR + nM * 9;
		Matrix3d  AM;
		AM << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

		double matxM[3];
		matxM[0] = *(imgpts + 2 * nMI);
		matxM[1] = *(imgpts + 2 * nMI + 1);
		matxM[2] = 1;

		Vector3d  bM(matxM);
		xM = AM.colPivHouseholderQr().solve(bM);			
		//Associate archor ray
		double* ptr2 = m_KR + nA * 9;
		Matrix3d  AA;
		AA << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

		double matxA[3];
		matxA[0] = *(imgpts + 2 * nAI);
		matxA[1] = *(imgpts + 2 * nAI + 1);
		matxA[2] = 1;

		Vector3d  bA(matxA);
		xA = AA.colPivHouseholderQr().solve(bA);			
	}

	//Parallax Angle
	double dDot = xM(0) * xA(0) + xM(1) * xA(1) + xM(2) * xA(2);			

	double dDisM = sqrt(xM(0) * xM(0) + xM(1) * xM(1) + xM(2) * xM(2));		
	double dDisA = sqrt(xA(0) * xA(0) + xA(1) * xA(1) + xA(2) * xA(2));		

	if (dDot / (dDisM * dDisA) > 1)
		feature[2] = 0;
	else if (dDot / (dDisM * dDisA) < -1)
		feature[2] = PI;
	else
	{
		double dw = acos(dDot / (dDisM * dDisA));
		feature[2] = dw;
	}


	double pti2k[3];
	pti2k[0] = *(camera + nA * 6 + 3) - *(camera + nM * 6 + 3);
	pti2k[1] = *(camera + nA * 6 + 4) - *(camera + nM * 6 + 4);
	pti2k[2] = *(camera + nA * 6 + 5) - *(camera + nM * 6 + 5);

	double dDot1 = xM[0] * pti2k[0] + xM[1] * pti2k[1] + xM[2] * pti2k[2];
	double dDisi2k = sqrt(pti2k[0] * pti2k[0] + pti2k[1] * pti2k[1] + pti2k[2] * pti2k[2]);
	double tmp = dDot1 / (dDisM * dDisi2k);
	double dW2;
	if (tmp > 1)
		dW2 = 0;
	if (tmp < -1)
		dW2 = PI;
	else
		dW2 = acos(tmp);
	return true;
}

bool triangulation::pba_initializeOtheArchors(
	double* imgpts,
	int* photo,
	double* camera,
	double* K,
	double* feature,
	int* archorSort,
	int nfeacout,
	int nOI,
	int FID)
{
	static int i = 0;
	double dw = feature[2];                   
	double dwNew;
	double dmaxw = dw;
	int   nNewI = 0;
	bool bAdjust = false;
	double dDot, dDisM, dDisA;

	if (dw < MAXARCHOR)
	{
		//current archor vector 
		int nO = photo[nOI];

		Vector3d  xO;
		if (m_bProvideXYZ)
		{
			xO(0) = m_XYZ[FID * 3] - *(camera + nO * 6 + 3);
			xO(1) = m_XYZ[FID * 3 + 1] - *(camera + nO * 6 + 4);
			xO(2) = m_XYZ[FID * 3 + 2] - *(camera + nO * 6 + 5);
		}
		else
		{
			double* ptr1 = m_KR + nO * 9;
			Matrix3d  AO;
			AO << ptr1[0], ptr1[1], ptr1[2], ptr1[3], ptr1[4], ptr1[5], ptr1[6], ptr1[7], ptr1[8];

			double matxO[3];
			matxO[0] = *(imgpts + nOI * 2);
			matxO[1] = *(imgpts + nOI * 2 + 1);
			matxO[2] = 1;

			Vector3d  bO(matxO);
			xO = AO.colPivHouseholderQr().solve(bO);
		}

		double dDAngle = atan2(xO(0), xO(2));
		double dHAngle = atan2(xO(1), sqrt(xO(0) * xO(0) + xO(2) * xO(2)));

		for (i = 0; i < nfeacout; i++)
		{
			//Main Archor Vector
			int nM = photo[i];
			Vector3d  xM;

			if (m_bProvideXYZ)
			{
				xM(0) = m_XYZ[FID * 3] - *(camera + nM * 6 + 3);
				xM(1) = m_XYZ[FID * 3 + 1] - *(camera + nM * 6 + 4);
				xM(2) = m_XYZ[FID * 3 + 2] - *(camera + nM * 6 + 5);
			}
			else
			{
				double* ptr2 = m_KR + nM * 9;
				Matrix3d  AM;
				AM << ptr2[0], ptr2[1], ptr2[2], ptr2[3], ptr2[4], ptr2[5], ptr2[6], ptr2[7], ptr2[8];

				double matxM[3];
				matxM[0] = *(imgpts + i * 2);
				matxM[1] = *(imgpts + i * 2 + 1);
				matxM[2] = 1;

				Vector3d  bM(matxM);
				xM = AM.colPivHouseholderQr().solve(bM);
			}

			//Parallax angle between current archor and main archor
			dDot = xM(0) * xO(0) + xM(1) * xO(1) + xM(2) * xO(2);
			dDisM = sqrt(xM(0) * xM(0) + xM(1) * xM(1) + xM(2) * xM(2));
			dDisA = sqrt(xO(0) * xO(0) + xO(1) * xO(1) + xO(2) * xO(2));

			if (dDot / (dDisM * dDisA) > 1)
				dwNew = 0;
			else if (dDot / (dDisM * dDisA) < -1)
				dwNew = PI;
			else
				dwNew = acos(dDot / (dDisM * dDisA));

			if (dwNew > dmaxw)
			{
				dmaxw = dwNew;
				archorSort[0] = nOI;
				archorSort[1] = i;
				feature[0] = dDAngle;
				feature[1] = dHAngle;
				feature[2] = dmaxw;
				bAdjust = true;
			}
		}
	}
	return bAdjust;
}
void triangulation::ba_constructP(double* P, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9];
		double matT[3];

		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;
			/*kappa phi omegaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];

			//test
			//double t1 = pP[0] * (2 - ptAngle[3]) + pP[1] * (2 - ptAngle[4]) + pP[2] * (2 - ptAngle[5]);
			//double t2 = pP[4] * (2 - ptAngle[3]) + pP[5] * (2 - ptAngle[4]) + pP[6] * (2 - ptAngle[5]);
			//double t3 = pP[8] * (2 - ptAngle[3]) + pP[9] * (2 - ptAngle[4]) + pP[10] * (2 - ptAngle[5]);

			//double u1 = t1 / t3;
			//double v1 = t2 / t3;
			//printf("%f %f\n", u1, v1);

			//t1 = pP[0] * 2 + pP[1] * 2 + pP[2] * 2 + pP[3];
			//t2 = pP[4] * 2 + pP[5] * 2 + pP[6] * 2 + pP[7];
			//t3 = pP[8] * 2 + pP[9] * 2 + pP[10] * 2 + pP[11];
			//double u2 = t1 / t3;
			//double v2 = t2 / t3;
			//printf("%f %f\n", u2, v2);
		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pP;
		double matR[9], matT[3];
		double K[9];
		memset(K, 0, 9 * sizeof(double));
		K[8] = 1;
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			matT[0] = -matR[0] * ptAngle[3] - matR[1] * ptAngle[4] - matR[2] * ptAngle[5];
			matT[1] = -matR[3] * ptAngle[3] - matR[4] * ptAngle[4] - matR[5] * ptAngle[5];
			matT[2] = -matR[6] * ptAngle[3] - matR[7] * ptAngle[4] - matR[8] * ptAngle[5];

			//KR

			K[0] = m_K[i * 2];
			K[4] = m_K[i * 2 + 1];

			//pP=K*[matR matT]
			pP = P + i * 12;
			pP[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pP[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pP[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pP[3] = K[0] * matT[0] + K[3] * matT[1] + K[6] * matT[2];
			pP[4] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pP[5] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pP[6] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pP[7] = K[1] * matT[0] + K[4] * matT[1] + K[7] * matT[2];
			pP[8] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pP[9] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pP[10] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];
			pP[11] = K[2] * matT[0] + K[5] * matT[1] + K[8] * matT[2];
		}
	}
}
void triangulation::ba_readCameraPoseration(char* fname, double* ical)
{
	FILE* fp;
	int  ch = EOF;

	if ((fopen_s(&fp,fname, "r")) != 0 && fp==NULL)
	{
		fprintf(stderr, "BA: Cannot open calbration file %s, exiting\n", fname);
		return;
	}
	if (tpe == bal) {
		double* ptr;
		for (int i = 0; i < m_ncams; i++)
		{
			ptr = ical + 9 * i;
			int num = fscanf_s(fp, "%lf", &ptr[0]);
			if (num != 1)
			{
				fprintf(stderr, "BA error: Format of Calibration file is wrong");
				return;
			}
			ptr[4] = ptr[0];
			ptr[1] = ptr[2] = ptr[3] = ptr[5] = ptr[6] = ptr[7] = 0;
			ptr[8] = 1;
			//printf("%f %f %f\n", ptr[0], ptr[1], ptr[2]);
			//printf("%f %f %f\n", ptr[3], ptr[4], ptr[5]);
			//printf("%f %f %f\n", ptr[6], ptr[7], ptr[8]);
		}
	}
	else if (tpe == colmap) {
		int s = 0;
		for (int i = 0; i < nc; i++) {
			//1-2-3 || 4-5-6 || 7-8-9
			int num = fscanf_s(fp, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", \
				& ical[9 * i + 0], &ical[9 * i + 3], &ical[9 * i + 6], \
				&ical[9 * i + 1], &ical[9 * i + 4], &ical[9 * i + 7],\
				&ical[9 * i + 2], &ical[9 * i + 5], &ical[9 * i + 8]);
			s += num;
			//printf("%f %f %f\n", ical[9 * i + 0], ical[9 * i + 3], ical[9 * i + 6]);
			//printf("%f %f %f\n", ical[9 * i + 1], ical[9 * i + 4], ical[9 * i + 7]);
			//printf("%f %f %f\n", ical[9 * i + 2], ical[9 * i + 5], ical[9 * i + 8]);
		}
		if (s != 9 * nc)
		{
			fprintf(stderr, "BA error: Format of Calibration file is wrong");
			return;
		}
	}

	fclose(fp);
}

void triangulation::ba_updateKR(double* KR, double* KdA, double* KdB, double* KdG, double* K, double* p)
{
	if (!m_bFocal)
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pKdA, * pKdB, * pKdG, * pK = NULL;
		double matR[9];
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];

		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;
			/*kappa phi omegaϵͳ*/
			//matR=matRG*matRB*matRA
			//ptAngle=[kappa,phi,omega]
			// double ey = ptAngle[0];
			// double ex = ptAngle[1];
			// double ez = ptAngle[2];
			// double c1 = cos(ey);   double c2 = cos(ex);   double c3 = cos(ez);
			// double s1 = sin(ey);   double s2 = sin(ex);   double s3 = sin(ez);
			// matR[0]=c1*c3-s1*s2*s3;     matR[1]=c2*s3;     matR[2]=s1*c3+c1*s2*s3;
			// matR[3]=-c1*s3-s1*s2*c3;    matR[4]=c2*c3;     matR[5]=-s1*s3+c1*s2*c3;
			// matR[6]=-s1*c2;             matR[7]=-s2;       matR[8]=c1*c2;
			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);

			//omega
			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			//phi
			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);

			//kappa
			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			//matRG
			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);

			//matRB
			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);

			//matRA
			matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//pKR=KR*matR
			pKR = KR + i * 9;
			if (tpe == bal) {
				pK = K + i * 9;
			}
			else if (tpe == colmap) {
				pK = K + (m_C[i] - 1) * 9;
			}
			//printf("%f %f %f\n", pK[0], pK[3], pK[6]);
			//printf("%f %f %f\n", pK[1], pK[4], pK[7]);
			//printf("%f %f %f\n", pK[2], pK[5], pK[8]);
			pKR[0] = pK[0] * matR[0] + pK[3] * matR[3] + pK[6] * matR[6];
			pKR[1] = pK[0] * matR[1] + pK[3] * matR[4] + pK[6] * matR[7];
			pKR[2] = pK[0] * matR[2] + pK[3] * matR[5] + pK[6] * matR[8];
			pKR[3] = pK[1] * matR[0] + pK[4] * matR[3] + pK[7] * matR[6];
			pKR[4] = pK[1] * matR[1] + pK[4] * matR[4] + pK[7] * matR[7];
			pKR[5] = pK[1] * matR[2] + pK[4] * matR[5] + pK[7] * matR[8];
			pKR[6] = pK[2] * matR[0] + pK[5] * matR[3] + pK[8] * matR[6];
			pKR[7] = pK[2] * matR[1] + pK[5] * matR[4] + pK[8] * matR[7];
			pKR[8] = pK[2] * matR[2] + pK[5] * matR[5] + pK[8] * matR[8];

			
			pKdG = KdG + i * 9;
			tmp1[0] = pK[0] * matDRG[0] + pK[3] * matDRG[3] + pK[6] * matDRG[6];
			tmp1[1] = pK[1] * matDRG[0] + pK[4] * matDRG[3] + pK[7] * matDRG[6];
			tmp1[2] = pK[2] * matDRG[0] + pK[5] * matDRG[3] + pK[8] * matDRG[6];
			tmp1[3] = pK[0] * matDRG[1] + pK[3] * matDRG[4] + pK[6] * matDRG[7];
			tmp1[4] = pK[1] * matDRG[1] + pK[4] * matDRG[4] + pK[7] * matDRG[7];
			tmp1[5] = pK[2] * matDRG[1] + pK[5] * matDRG[4] + pK[8] * matDRG[7];
			tmp1[6] = pK[0] * matDRG[2] + pK[3] * matDRG[5] + pK[6] * matDRG[8];
			tmp1[7] = pK[1] * matDRG[2] + pK[4] * matDRG[5] + pK[7] * matDRG[8];
			tmp1[8] = pK[2] * matDRG[2] + pK[5] * matDRG[5] + pK[8] * matDRG[8];

			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdG[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdG[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdG[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdG[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdG[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdG[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdG[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdG[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdG[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

	
			pKdB = KdB + i * 9;
			tmp1[0] = pK[0] * matRG[0] + pK[3] * matRG[3] + pK[6] * matRG[6];
			tmp1[1] = pK[1] * matRG[0] + pK[4] * matRG[3] + pK[7] * matRG[6];
			tmp1[2] = pK[2] * matRG[0] + pK[5] * matRG[3] + pK[8] * matRG[6];
			tmp1[3] = pK[0] * matRG[1] + pK[3] * matRG[4] + pK[6] * matRG[7];
			tmp1[4] = pK[1] * matRG[1] + pK[4] * matRG[4] + pK[7] * matRG[7];
			tmp1[5] = pK[2] * matRG[1] + pK[5] * matRG[4] + pK[8] * matRG[7];
			tmp1[6] = pK[0] * matRG[2] + pK[3] * matRG[5] + pK[6] * matRG[8];
			tmp1[7] = pK[1] * matRG[2] + pK[4] * matRG[5] + pK[7] * matRG[8];
			tmp1[8] = pK[2] * matRG[2] + pK[5] * matRG[5] + pK[8] * matRG[8];

			tmp2[0] = tmp1[0] * matDRB[0] + tmp1[3] * matDRB[3] + tmp1[6] * matDRB[6];
			tmp2[1] = tmp1[1] * matDRB[0] + tmp1[4] * matDRB[3] + tmp1[7] * matDRB[6];
			tmp2[2] = tmp1[2] * matDRB[0] + tmp1[5] * matDRB[3] + tmp1[8] * matDRB[6];
			tmp2[3] = tmp1[0] * matDRB[1] + tmp1[3] * matDRB[4] + tmp1[6] * matDRB[7];
			tmp2[4] = tmp1[1] * matDRB[1] + tmp1[4] * matDRB[4] + tmp1[7] * matDRB[7];
			tmp2[5] = tmp1[2] * matDRB[1] + tmp1[5] * matDRB[4] + tmp1[8] * matDRB[7];
			tmp2[6] = tmp1[0] * matDRB[2] + tmp1[3] * matDRB[5] + tmp1[6] * matDRB[8];
			tmp2[7] = tmp1[1] * matDRB[2] + tmp1[4] * matDRB[5] + tmp1[7] * matDRB[8];
			tmp2[8] = tmp1[2] * matDRB[2] + tmp1[5] * matDRB[5] + tmp1[8] * matDRB[8];

			pKdB[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdB[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdB[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdB[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdB[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdB[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdB[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdB[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdB[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

		
			pKdA = KdA + i * 9;
			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdA[0] = tmp2[0] * matDRA[0] + tmp2[3] * matDRA[3] + tmp2[6] * matDRA[6];
			pKdA[3] = tmp2[1] * matDRA[0] + tmp2[4] * matDRA[3] + tmp2[7] * matDRA[6];
			pKdA[6] = tmp2[2] * matDRA[0] + tmp2[5] * matDRA[3] + tmp2[8] * matDRA[6];
			pKdA[1] = tmp2[0] * matDRA[1] + tmp2[3] * matDRA[4] + tmp2[6] * matDRA[7];
			pKdA[4] = tmp2[1] * matDRA[1] + tmp2[4] * matDRA[4] + tmp2[7] * matDRA[7];
			pKdA[7] = tmp2[2] * matDRA[1] + tmp2[5] * matDRA[4] + tmp2[8] * matDRA[7];
			pKdA[2] = tmp2[0] * matDRA[2] + tmp2[3] * matDRA[5] + tmp2[6] * matDRA[8];
			pKdA[5] = tmp2[1] * matDRA[2] + tmp2[4] * matDRA[5] + tmp2[7] * matDRA[8];
			pKdA[8] = tmp2[2] * matDRA[2] + tmp2[5] * matDRA[5] + tmp2[8] * matDRA[8];

		}
	}
	else
	{
		int i = 0;
		double* ptAngle;
		double* pKR, * pKdA, * pKdB, * pKdG;
		double matR[9];
		double matRG[9], matRB[9], matRA[9];
		double matDRG[9], matDRB[9], matDRA[9];
		double tmp1[9], tmp2[9];
		double K[9];
		memset(K, 0, 9 * sizeof(double));
		K[8] = 1;
		for (i = 0; i < m_ncams; i++)
		{
			ptAngle = p + i * 6;

			matR[0] = cos(ptAngle[1]) * cos(ptAngle[0]);
			matR[1] = cos(ptAngle[1]) * sin(ptAngle[0]);
			matR[2] = -sin(ptAngle[1]);
			matR[3] = sin(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) - cos(ptAngle[2]) * sin(ptAngle[0]);
			matR[4] = sin(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) + cos(ptAngle[2]) * cos(ptAngle[0]);
			matR[5] = sin(ptAngle[2]) * cos(ptAngle[1]);
			matR[6] = cos(ptAngle[2]) * sin(ptAngle[1]) * cos(ptAngle[0]) + sin(ptAngle[2]) * sin(ptAngle[0]);
			matR[7] = cos(ptAngle[2]) * sin(ptAngle[1]) * sin(ptAngle[0]) - sin(ptAngle[2]) * cos(ptAngle[0]);
			matR[8] = cos(ptAngle[2]) * cos(ptAngle[1]);


			matRG[0] = 1;		matRG[1] = 0;				matRG[2] = 0;
			matRG[3] = 0;		matRG[4] = cos(ptAngle[2]);	matRG[5] = sin(ptAngle[2]);
			matRG[6] = 0;		matRG[7] = -sin(ptAngle[2]);	matRG[8] = cos(ptAngle[2]);

			matRB[0] = cos(ptAngle[1]);		matRB[1] = 0;		matRB[2] = -sin(ptAngle[1]);
			matRB[3] = 0;					matRB[4] = 1;		matRB[5] = 0;
			matRB[6] = sin(ptAngle[1]);		matRB[7] = 0;		matRB[8] = cos(ptAngle[1]);
	
			matRA[0] = cos(ptAngle[0]);		matRA[1] = sin(ptAngle[0]);			matRA[2] = 0;
			matRA[3] = -sin(ptAngle[0]);	matRA[4] = cos(ptAngle[0]);			matRA[5] = 0;
			matRA[6] = 0;					matRA[7] = 0;						matRA[8] = 1;

			matDRG[0] = 0;		matDRG[1] = 0;			matDRG[2] = 0;
			matDRG[3] = 0;		matDRG[4] = -sin(ptAngle[2]);	matDRG[5] = cos(ptAngle[2]);
			matDRG[6] = 0;		matDRG[7] = -cos(ptAngle[2]);	matDRG[8] = -sin(ptAngle[2]);

			matDRB[0] = -sin(ptAngle[1]);		matDRB[1] = 0;		matDRB[2] = -cos(ptAngle[1]);
			matDRB[3] = 0;						matDRB[4] = 0;		matDRB[5] = 0;
			matDRB[6] = cos(ptAngle[1]);		matDRB[7] = 0;		matDRB[8] = -sin(ptAngle[1]);

			matDRA[0] = -sin(ptAngle[0]);		matDRA[1] = cos(ptAngle[0]);		matDRA[2] = 0;
			matDRA[3] = -cos(ptAngle[0]);		matDRA[4] = -sin(ptAngle[0]);		matDRA[5] = 0;
			matDRA[6] = 0;						matDRA[7] = 0;						matDRA[8] = 0;

			//KR

			K[0] = m_K[i * 3];
			K[4] = m_K[i * 3];

			pKR = KR + i * 9;
			pKR[0] = K[0] * matR[0] + K[3] * matR[3] + K[6] * matR[6];
			pKR[1] = K[0] * matR[1] + K[3] * matR[4] + K[6] * matR[7];
			pKR[2] = K[0] * matR[2] + K[3] * matR[5] + K[6] * matR[8];
			pKR[3] = K[1] * matR[0] + K[4] * matR[3] + K[7] * matR[6];
			pKR[4] = K[1] * matR[1] + K[4] * matR[4] + K[7] * matR[7];
			pKR[5] = K[1] * matR[2] + K[4] * matR[5] + K[7] * matR[8];
			pKR[6] = K[2] * matR[0] + K[5] * matR[3] + K[8] * matR[6];
			pKR[7] = K[2] * matR[1] + K[5] * matR[4] + K[8] * matR[7];
			pKR[8] = K[2] * matR[2] + K[5] * matR[5] + K[8] * matR[8];

			//KdG
			pKdG = KdG + i * 9;
			tmp1[0] = K[0] * matDRG[0] + K[3] * matDRG[3] + K[6] * matDRG[6];
			tmp1[1] = K[1] * matDRG[0] + K[4] * matDRG[3] + K[7] * matDRG[6];
			tmp1[2] = K[2] * matDRG[0] + K[5] * matDRG[3] + K[8] * matDRG[6];
			tmp1[3] = K[0] * matDRG[1] + K[3] * matDRG[4] + K[6] * matDRG[7];
			tmp1[4] = K[1] * matDRG[1] + K[4] * matDRG[4] + K[7] * matDRG[7];
			tmp1[5] = K[2] * matDRG[1] + K[5] * matDRG[4] + K[8] * matDRG[7];
			tmp1[6] = K[0] * matDRG[2] + K[3] * matDRG[5] + K[6] * matDRG[8];
			tmp1[7] = K[1] * matDRG[2] + K[4] * matDRG[5] + K[7] * matDRG[8];
			tmp1[8] = K[2] * matDRG[2] + K[5] * matDRG[5] + K[8] * matDRG[8];

			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdG[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdG[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdG[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdG[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdG[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdG[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdG[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdG[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdG[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//KdB
			pKdB = KdB + i * 9;
			tmp1[0] = K[0] * matRG[0] + K[3] * matRG[3] + K[6] * matRG[6];
			tmp1[1] = K[1] * matRG[0] + K[4] * matRG[3] + K[7] * matRG[6];
			tmp1[2] = K[2] * matRG[0] + K[5] * matRG[3] + K[8] * matRG[6];
			tmp1[3] = K[0] * matRG[1] + K[3] * matRG[4] + K[6] * matRG[7];
			tmp1[4] = K[1] * matRG[1] + K[4] * matRG[4] + K[7] * matRG[7];
			tmp1[5] = K[2] * matRG[1] + K[5] * matRG[4] + K[8] * matRG[7];
			tmp1[6] = K[0] * matRG[2] + K[3] * matRG[5] + K[6] * matRG[8];
			tmp1[7] = K[1] * matRG[2] + K[4] * matRG[5] + K[7] * matRG[8];
			tmp1[8] = K[2] * matRG[2] + K[5] * matRG[5] + K[8] * matRG[8];

			tmp2[0] = tmp1[0] * matDRB[0] + tmp1[3] * matDRB[3] + tmp1[6] * matDRB[6];
			tmp2[1] = tmp1[1] * matDRB[0] + tmp1[4] * matDRB[3] + tmp1[7] * matDRB[6];
			tmp2[2] = tmp1[2] * matDRB[0] + tmp1[5] * matDRB[3] + tmp1[8] * matDRB[6];
			tmp2[3] = tmp1[0] * matDRB[1] + tmp1[3] * matDRB[4] + tmp1[6] * matDRB[7];
			tmp2[4] = tmp1[1] * matDRB[1] + tmp1[4] * matDRB[4] + tmp1[7] * matDRB[7];
			tmp2[5] = tmp1[2] * matDRB[1] + tmp1[5] * matDRB[4] + tmp1[8] * matDRB[7];
			tmp2[6] = tmp1[0] * matDRB[2] + tmp1[3] * matDRB[5] + tmp1[6] * matDRB[8];
			tmp2[7] = tmp1[1] * matDRB[2] + tmp1[4] * matDRB[5] + tmp1[7] * matDRB[8];
			tmp2[8] = tmp1[2] * matDRB[2] + tmp1[5] * matDRB[5] + tmp1[8] * matDRB[8];

			pKdB[0] = tmp2[0] * matRA[0] + tmp2[3] * matRA[3] + tmp2[6] * matRA[6];
			pKdB[3] = tmp2[1] * matRA[0] + tmp2[4] * matRA[3] + tmp2[7] * matRA[6];
			pKdB[6] = tmp2[2] * matRA[0] + tmp2[5] * matRA[3] + tmp2[8] * matRA[6];
			pKdB[1] = tmp2[0] * matRA[1] + tmp2[3] * matRA[4] + tmp2[6] * matRA[7];
			pKdB[4] = tmp2[1] * matRA[1] + tmp2[4] * matRA[4] + tmp2[7] * matRA[7];
			pKdB[7] = tmp2[2] * matRA[1] + tmp2[5] * matRA[4] + tmp2[8] * matRA[7];
			pKdB[2] = tmp2[0] * matRA[2] + tmp2[3] * matRA[5] + tmp2[6] * matRA[8];
			pKdB[5] = tmp2[1] * matRA[2] + tmp2[4] * matRA[5] + tmp2[7] * matRA[8];
			pKdB[8] = tmp2[2] * matRA[2] + tmp2[5] * matRA[5] + tmp2[8] * matRA[8];

			//KdA
			pKdA = KdA + i * 9;
			tmp2[0] = tmp1[0] * matRB[0] + tmp1[3] * matRB[3] + tmp1[6] * matRB[6];
			tmp2[1] = tmp1[1] * matRB[0] + tmp1[4] * matRB[3] + tmp1[7] * matRB[6];
			tmp2[2] = tmp1[2] * matRB[0] + tmp1[5] * matRB[3] + tmp1[8] * matRB[6];
			tmp2[3] = tmp1[0] * matRB[1] + tmp1[3] * matRB[4] + tmp1[6] * matRB[7];
			tmp2[4] = tmp1[1] * matRB[1] + tmp1[4] * matRB[4] + tmp1[7] * matRB[7];
			tmp2[5] = tmp1[2] * matRB[1] + tmp1[5] * matRB[4] + tmp1[8] * matRB[7];
			tmp2[6] = tmp1[0] * matRB[2] + tmp1[3] * matRB[5] + tmp1[6] * matRB[8];
			tmp2[7] = tmp1[1] * matRB[2] + tmp1[4] * matRB[5] + tmp1[7] * matRB[8];
			tmp2[8] = tmp1[2] * matRB[2] + tmp1[5] * matRB[5] + tmp1[8] * matRB[8];

			pKdA[0] = tmp2[0] * matDRA[0] + tmp2[3] * matDRA[3] + tmp2[6] * matDRA[6];
			pKdA[3] = tmp2[1] * matDRA[0] + tmp2[4] * matDRA[3] + tmp2[7] * matDRA[6];
			pKdA[6] = tmp2[2] * matDRA[0] + tmp2[5] * matDRA[3] + tmp2[8] * matDRA[6];
			pKdA[1] = tmp2[0] * matDRA[1] + tmp2[3] * matDRA[4] + tmp2[6] * matDRA[7];
			pKdA[4] = tmp2[1] * matDRA[1] + tmp2[4] * matDRA[4] + tmp2[7] * matDRA[7];
			pKdA[7] = tmp2[2] * matDRA[1] + tmp2[5] * matDRA[4] + tmp2[8] * matDRA[7];
			pKdA[2] = tmp2[0] * matDRA[2] + tmp2[3] * matDRA[5] + tmp2[6] * matDRA[8];
			pKdA[5] = tmp2[1] * matDRA[2] + tmp2[4] * matDRA[5] + tmp2[7] * matDRA[8];
			pKdA[8] = tmp2[2] * matDRA[2] + tmp2[5] * matDRA[5] + tmp2[8] * matDRA[8];

		}
	}
}
void triangulation::ba_readProjectionAndTriangulateFeature(FILE* fp, double* projs, int ncams)
{
	int n;
	int nframes;
	int ptno = 0;
	int frameno;
	double* ptr1 = projs;

	int i;
	//read all projection point, initialize three feature angle at the same time

	while (!feof(fp))
	{
		n = readNInts(fp, &nframes, 1); 
		if (n != 1)
			break;

		Eigen::MatrixXd A(2 * nframes, 4);
		//if (nframes > 3)
		//	nframes = 3;


		for (i = 0; i < nframes; ++i)
		{
			n = readNInts(fp, &frameno, 1); 

			if (frameno >= ncams){
				fprintf(stderr, "BA: the image No. of projection point is out of max image No.\n");
				return;
			}

			n += readNDoubles(fp, ptr1, 2); 

			if (n != 3){
				fprintf(stderr, "BA:reading image projections wrong!\n");
				return;
			}

			const Eigen::Vector2d pt(ptr1[0], ptr1[1]);
			double* ptr = m_P + frameno * 12;
			const Eigen::Matrix<double, 3, 4> Pmat = (Eigen::Matrix<double, 3, 4>() <<
				ptr[0], ptr[1], ptr[2], ptr[3],
				ptr[4], ptr[5], ptr[6], ptr[7],
				ptr[8], ptr[9], ptr[10], ptr[11]).finished();
			A.row(2 * i) = pt(0) * Pmat.row(2) - Pmat.row(0);
			A.row(2 * i + 1) = pt(1) * Pmat.row(2) - Pmat.row(1);
		}

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullV);
		Eigen::Vector4d X = svd.matrixV().col(3);
		X /= X(3);

		(m_motstruct + m_ncams * 6 + ptno * 3)[0] = X[0];
		(m_motstruct + m_ncams * 6 + ptno * 3)[1] = X[1];
		(m_motstruct + m_ncams * 6 + ptno * 3)[2] = X[2];
		ptno++;
	}
}