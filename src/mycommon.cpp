#include "mycommon.h"

//mycommon::mycommon()
//{
//
//}
//mycommon::~mycommon()
//{
//}
void mycommon::undist(std::vector<cv::Mat> K, std::vector<cv::Mat> distCoeffs, \
	std::vector<cv::Point2d> distortedPoints, \
	std::vector<cv::Point2d> & undistortedPixelPoints, \
	vector<int> camid_group) {

	std::vector<int> tmp, cid_gr;
	int tmp1 = 1, s = 0;
	cid_gr.push_back(tmp1);
	for (const auto& cm : camid_group) {
		if (cm != tmp1) {
			tmp1 = cm;
			tmp.push_back(s);
			cid_gr.push_back(cm);
		}
		++s;
	}
	tmp.push_back(s);

 	std::vector<cv::Point2d> undistortedNormPoints;
	std::vector<cv::Point2d> sub_distortedPoints; 
	int s1 = 0, s2 = 0;
	for (int i = 0; i < tmp.size(); i++) {
 
		int ki = cid_gr[i] - 1;

		if (i > 0) {
			s1 = tmp[i - 1];
		}
		s2 = tmp[i];

		std::copy(
			distortedPoints.begin() + s1,
			distortedPoints.begin() + s2,    
			std::back_inserter(sub_distortedPoints)
		);
		//printf("%f %f\n", sub_distortedPoints[0].x, sub_distortedPoints[0].y);
		//int hh = sub_distortedPoints.size() - 1;
		//printf("%f %f %f %f\n", sub_distortedPoints.back().x, sub_distortedPoints.back().y, sub_distortedPoints[hh].x, sub_distortedPoints[hh].y);
		cv::undistortPoints(sub_distortedPoints, undistortedNormPoints, K[ki], distCoeffs[ki]);
		sub_distortedPoints.clear();
		for (const auto& pt : undistortedNormPoints) {
			float u = static_cast<float>(K[ki].at<double>(0, 0) * pt.x + K[ki].at<double>(0, 2));
			float v = static_cast<float>(K[ki].at<double>(1, 1) * pt.y + K[ki].at<double>(1, 2));
			undistortedPixelPoints.emplace_back(u, v);
		}
		undistortedNormPoints.clear();
	}
}
int mycommon::readNInts(FILE* fp, int* vals, int nvals)
{
	int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i) {
		j = fscanf(fp, "%d", vals + i);
		if (j == EOF) return EOF;

		if (j != 1 || ferror(fp)) return EOF - 1;

		n += j;
	}
	
	return n;
}

int mycommon::readNDoubles(FILE* fp, double* vals, int nvals)
{
	int i;
	int n, j;

	for (i = n = 0; i < nvals; ++i)
	{
		j = fscanf(fp, "%lf", vals + i);
		if (j == EOF) return EOF;

		if (j != 1 || ferror(fp)) return EOF - 1;

		n += j;
	}

	return n;//????2
}

void mycommon::QuaternionToAngleAxis(double* quaternion, double* angle_axis)
{
	double& q1 = quaternion[1];
	double& q2 = quaternion[2];
	double& q3 = quaternion[3];
	double sin_squared_theta = q1 * q1 + q2 * q2 + q3 * q3;

	// For quaternions representing non-zero rotation, the conversion
	// is numerically stable.
	if (sin_squared_theta > 0.0) {
		double sin_theta = sqrt(sin_squared_theta);
		double& cos_theta = quaternion[0];

		// If cos_theta is negative, theta is greater than pi/2, which
		// means that angle for the angle_axis vector which is 2 * theta
		// would be greater than pi.
		//
		// While this will result in the correct rotation, it does not
		// result in a normalized angle-axis vector.
		//
		// In that case we observe that 2 * theta ~ 2 * theta - 2 * pi,
		// which is equivalent saying
		//
		//   theta - pi = atan(sin(theta - pi), cos(theta - pi))
		//              = atan(-sin(theta), -cos(theta))
		//
		double two_theta =
			2.0 * ((cos_theta < 0.0)
				? atan2(-sin_theta, -cos_theta)
				: atan2(sin_theta, cos_theta));
		double k = two_theta / sin_theta;
		angle_axis[0] = q1 * k;
		angle_axis[1] = q2 * k;
		angle_axis[2] = q3 * k;
	}
	else {
		// For zero rotation, sqrt() will produce NaN in the derivative since
		// the argument is zero.  By approximating with a Taylor series,
		// and truncating at one term, the value and first derivatives will be
		// computed correctly when Jets are used.
		double k(2.0);
		angle_axis[0] = q1 * k;
		angle_axis[1] = q2 * k;
		angle_axis[2] = q3 * k;
	}
}
void mycommon::RotationMatrixToQuaternion(double* R, double* quaternion)
{
	double trace = R[0] + R[4] + R[8];
	if (trace >= 0.0)
	{
		double t = sqrt(trace + 1.0);
		quaternion[0] = 0.5 * t;
		t = 0.5 / t;
		quaternion[1] = (R[7] - R[5]) * t;
		quaternion[2] = (R[2] - R[6]) * t;
		quaternion[3] = (R[3] - R[1]) * t;
	}
	else
	{
		int i = 0;
		if (R[4] > R[0]) {
			i = 1;
		}

		if (R[8] > R[4 * i]) {
			i = 2;
		}

		const int j = (i + 1) % 3;
		const int k = (j + 1) % 3;
		double t = sqrt(R[4 * i] - R[4 * j] - R[4 * k] + 1.0);
		quaternion[i + 1] = 0.5 * t;
		t = 0.5 / t;
		quaternion[0] = (R[3 * k + j] - R[3 * j + k]) * t;
		quaternion[j + 1] = (R[3 * j + i] + R[3 * i + j]) * t;
		quaternion[k + 1] = (R[3 * k + i] + R[3 * i + k]) * t;
	}
}
void mycommon::rotationMatrixToEulerAngles(double* R, double* eulerAngles)
{
	//assert(isRotationMatrix(R));
	double sy = sqrt(R[0] * R[0] + R[1] * R[1]);

	bool singular = sy < 1e-6;

	double phi, omega, kappa;
	if (!singular)
	{
		phi = atan2(-R[2], sy);
		omega = atan2(R[5], R[8]);
		kappa = atan2(R[1], R[0]);
	}
	else
	{
		phi = 0;
		omega = atan2(R[5], R[8]);
		kappa = atan2(R[1], R[0]);
	}
	eulerAngles[0] = kappa;
	eulerAngles[1] = phi;
	eulerAngles[2] = omega;
}
void mycommon::rotationMatrixToEulerAngles_phi_omega_kappa(double* R, double* eulerAngles)
{
	//assert(isRotationMatrix(R));
	double sy = sqrt(R[1] * R[1] + R[4] * R[4]);

	bool singular = sy < 1e-6;

	double phi, omega, kappa;
	if (!singular)
	{
		phi = atan2(-R[6], R[8]);
		omega = atan2(-R[7], sy);
		kappa = atan2(R[1], R[4]);
	}
	else
	{
		phi = atan2(-R[3], R[0]);
		omega = atan2(-R[7], sy);
		kappa = 0;
	}
	eulerAngles[0] = phi;
	eulerAngles[1] = omega;
	eulerAngles[2] = kappa;
}
void mycommon::eulerAnglesToRotationMatrix(double* eulerAngles, double* R)
{
	double kappa = eulerAngles[0], phi = eulerAngles[1], omega = eulerAngles[2];
	double R_x[9], R_y[9], R_z[9];
	R_x[0] = 1;   R_x[1] = 0;               R_x[2] = 0;
	R_x[3] = 0;   R_x[4] = cos(omega);      R_x[5] = sin(omega);
	R_x[6] = 0;   R_x[7] = -sin(omega);     R_x[8] = cos(omega);

	R_y[0] = cos(phi);   R_y[1] = 0;       R_y[2] = -sin(phi);
	R_y[3] = 0;          R_y[4] = 1;       R_y[5] = 0;
	R_y[6] = sin(phi);   R_y[7] = 0;       R_y[8] = cos(phi);

	R_z[0] = cos(kappa);  R_z[1] = sin(kappa); R_z[2] = 0;
	R_z[3] = -sin(kappa); R_z[4] = cos(kappa); R_z[5] = 0;
	R_z[6] = 0;           R_z[7] = 0;          R_z[8] = 1;

	double tmp[9];
	tmp[0] = R_y[0] * R_z[0] + R_y[1] * R_z[3] + R_y[2] * R_z[6];
	tmp[1] = R_y[0] * R_z[1] + R_y[1] * R_z[4] + R_y[2] * R_z[7];
	tmp[2] = R_y[0] * R_z[2] + R_y[1] * R_z[5] + R_y[2] * R_z[8];

	tmp[3] = R_y[3] * R_z[0] + R_y[4] * R_z[3] + R_y[5] * R_z[6];
	tmp[4] = R_y[3] * R_z[1] + R_y[4] * R_z[4] + R_y[5] * R_z[7];
	tmp[5] = R_y[3] * R_z[2] + R_y[4] * R_z[5] + R_y[5] * R_z[8];

	tmp[6] = R_y[6] * R_z[0] + R_y[7] * R_z[3] + R_y[8] * R_z[6];
	tmp[7] = R_y[6] * R_z[1] + R_y[7] * R_z[4] + R_y[8] * R_z[7];
	tmp[8] = R_y[6] * R_z[2] + R_y[7] * R_z[5] + R_y[8] * R_z[8];

	R[0] = R_x[0] * tmp[0] + R_x[1] * tmp[3] + R_x[2] * tmp[6];
	R[1] = R_x[0] * tmp[1] + R_x[1] * tmp[4] + R_x[2] * tmp[7];
	R[2] = R_x[0] * tmp[2] + R_x[1] * tmp[5] + R_x[2] * tmp[8];

	R[3] = R_x[3] * tmp[0] + R_x[4] * tmp[3] + R_x[5] * tmp[6];
	R[4] = R_x[3] * tmp[1] + R_x[4] * tmp[4] + R_x[5] * tmp[7];
	R[5] = R_x[3] * tmp[2] + R_x[4] * tmp[5] + R_x[5] * tmp[8];

	R[6] = R_x[6] * tmp[0] + R_x[7] * tmp[3] + R_x[8] * tmp[6];
	R[7] = R_x[6] * tmp[1] + R_x[7] * tmp[4] + R_x[8] * tmp[7];
	R[8] = R_x[6] * tmp[2] + R_x[7] * tmp[5] + R_x[8] * tmp[8];
}
void mycommon::eulerAnglesToRotationMatrix_phi_omega_kappa(double* eulerAngles, double* R){
	double phi = eulerAngles[0], omega = eulerAngles[1], kappa = eulerAngles[2];
	double ey = phi, ex = omega, ez = kappa;
	R[0]=cos(ey)*cos(ez)-sin(ey)*sin(ex)*sin(ez);
	R[1]=cos(ex)*sin(ez);
	R[2]=sin(ey)*cos(ez)+cos(ey)*sin(ex)*sin(ez);
	R[3]=-cos(ey)*sin(ez)-sin(ey)*sin(ex)*cos(ez);
	R[4]=cos(ex)*cos(ez);
	R[5]=-sin(ey)*sin(ez)+cos(ey)*sin(ex)*cos(ez);
	R[6]=-sin(ey)*cos(ex);
	R[7]=-sin(ex);
	R[8]=cos(ey)*cos(ex);
}
double DotProduct(double x[3], double y[3]) {
	return (x[0] * y[0] + x[1] * y[1] + x[2] * y[2]);
}
void mycommon::AngleAxisRotatePoint(double angle_axis[3], double pt[3], double result[3]) {
	double theta2 = DotProduct(angle_axis, angle_axis);
	if (theta2 > std::numeric_limits<double>::epsilon()) {
		// Away from zero, use the rodriguez formula
		//
		//   result = pt costheta +
		//            (w x pt) * sintheta +
		//            w (w . pt) (1 - costheta)
		//
		// We want to be careful to only evaluate the square root if the
		// norm of the angle_axis vector is greater than zero. Otherwise
		// we get a division by zero.
		//
		double theta = sqrt(theta2);
		double costheta = cos(theta);
		double sintheta = sin(theta);
		double theta_inverse = 1.0 / theta;

		double w[3] = { angle_axis[0] * theta_inverse,
						 angle_axis[1] * theta_inverse,
						 angle_axis[2] * theta_inverse };

		// Explicitly inlined evaluation of the cross product for
		// performance reasons.
		double w_cross_pt[3] = { w[1] * pt[2] - w[2] * pt[1],
								  w[2] * pt[0] - w[0] * pt[2],
								  w[0] * pt[1] - w[1] * pt[0] };
		double tmp =
			(w[0] * pt[0] + w[1] * pt[1] + w[2] * pt[2]) * (1.0 - costheta);

		result[0] = pt[0] * costheta + w_cross_pt[0] * sintheta + w[0] * tmp;
		result[1] = pt[1] * costheta + w_cross_pt[1] * sintheta + w[1] * tmp;
		result[2] = pt[2] * costheta + w_cross_pt[2] * sintheta + w[2] * tmp;
	}
	else {
		// Near zero, the first order Taylor approximation of the rotation
		// matrix R corresponding to a vector w and angle w is
		//
		//   R = I + hat(w) * sin(theta)
		//
		// But sintheta ~ theta and theta * w = angle_axis, which gives us
		//
		//  R = I + hat(w)
		//
		// and actually performing multiplication with the point pt, gives us
		// R * pt = pt + w x pt.
		//
		// Switching to the Taylor expansion near zero provides meaningful
		// derivatives when evaluated using Jets.
		//
		// Explicitly inlined evaluation of the cross product for
		// performance reasons.
		double w_cross_pt[3] = { angle_axis[1] * pt[2] - angle_axis[2] * pt[1],
								  angle_axis[2] * pt[0] - angle_axis[0] * pt[2],
								  angle_axis[0] * pt[1] - angle_axis[1] * pt[0] };

		result[0] = pt[0] + w_cross_pt[0];
		result[1] = pt[1] + w_cross_pt[1];
		result[2] = pt[2] + w_cross_pt[2];
	}
}
int mycommon::findNcameras(FILE* fp)
{
	int lineno, ncams, ch;

	lineno = ncams = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#') { /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);

		SKIP_LINE(fp);
		++lineno;
		if (ferror(fp))
		{
			fprintf(stderr, "findNcameras(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		++ncams;
	}
	return ncams;
}
void mycommon::PLY()
{
	//????ply???
	FILE* fp1 = NULL, * fp2 = NULL;
	std::string originalPath(dt);
	size_t pos = originalPath.find_last_of("/\\");
	std::string parentPath = (pos != std::string::npos) ? originalPath.substr(0, pos) : "";
	std::string newPath1 = parentPath + "/" + "Cam.txt";
	std::string newPath2 = parentPath + "/" + "XYZ1.txt";
	std::string newPath3 = parentPath + "/" + "3D.ply";
	printf("DataPreprocess: Version 1.0\n");
	FILE* fp;
	char* szCam = const_cast<char*>(newPath1.c_str());
	char* szXYZ = const_cast<char*>(newPath2.c_str());
	char* sz3D = const_cast<char*>(newPath3.c_str());
	fp1 = fopen(szCam, "r");
	int ncams = findNcameras(fp1);
	fp2 = fopen(szXYZ, "r");
	int npts = findNcameras(fp2);
	rewind(fp1);
	rewind(fp2);
	fp = fopen(sz3D, "w");
	fprintf(fp, "%s\n%s\n%s%d\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n", "ply", "format ascii 1.0", "element vertex ", ncams + npts, \
		"property float x", "property float y", "property float z", \
		"property uchar red", "property uchar green", "property uchar blue", "end_header");
	double* motstructs = (double*)malloc((ncams * 3 + npts * 3) * sizeof(double));
	double* ptr1 = motstructs;
	double kappa = 0, phi = 0, omega = 0;
	while (!feof(fp1))
	{
		fscanf(fp1, "%lf %lf %lf %lf %lf %lf", &kappa, &phi, &omega, &ptr1[0], &ptr1[1], &ptr1[2]);
		fprintf(fp, "%lf %lf %lf %d %d %d\n", ptr1[0], ptr1[1], ptr1[2], 0, 255, 0);
		ptr1 += 3;
	}
	while (!feof(fp2))
	{
		fscanf(fp2, "%lf %lf %lf", &ptr1[0], &ptr1[1], &ptr1[2]);
		fprintf(fp, "%lf %lf %lf %d %d %d\n", ptr1[0], ptr1[1], ptr1[2], 255, 255, 255);
		ptr1 += 3;
	}

	int n = readNDoubles(fp, ptr1, 1);
	fclose(fp);
	free(motstructs);
}
int mycommon::countNDoubles(FILE* fp)
{
	int lineno, ch, np, i;
	char buf[MAXSTRLEN], * s;
	double dummy;

	lineno = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) return 0;

		ungetc(ch, fp);
		++lineno;
		if (!fgets(buf, MAXSTRLEN - 1, fp)) { /* read the line found... */
			fprintf(stderr, "countNDoubles(): error reading input file, line %d\n", lineno);
			exit(1);
		}
		/* ...and count the number of doubles it has */
		for (np = i = 0, s = buf; 1; ++np, s += i) {
			ch = sscanf(s, "%lf%n", &dummy, &i);
			if (ch == 0 || ch == EOF) break;
		}

		rewind(fp);
		return np;
	}
	return 0; // should not reach this point
}
void mycommon::readNpointsAndNprojections(FILE* fp, int* n3Dpts, int pnp, int* nprojs, int mnp)
{
	int nfirst, lineno, npts, nframes, ch, n;

	/* #parameters for the first line */
	nfirst = countNDoubles(fp);

	*n3Dpts = *nprojs = lineno = npts = 0;
	while (!feof(fp))
	{
		if ((ch = fgetc(fp)) == '#')
		{ /* skip comments */
			SKIP_LINE(fp);
			++lineno;
			continue;
		}

		if (feof(fp)) break;

		ungetc(ch, fp);
		++lineno;
		//skipNDoubles(fp, pnp);
		n = readNInts(fp, &nframes, 1);

		if (n != 1)
			exit(1);

		//printf("%d ", nframes);

		SKIP_LINE(fp);
		*nprojs += nframes;
		++npts;
	}

	*n3Dpts = npts;
}

Matrix3d mycommon::ExpSO3(const Vector3d& phi) {
	double theta = phi.norm();
	if (theta < 1e-6) return Matrix3d::Identity();
	Vector3d u = phi / theta;
	Matrix3d u_hat;
	u_hat << 0, -u.z(), u.y(),
		u.z(), 0, -u.x(),
		-u.y(), u.x(), 0;
	return Matrix3d::Identity() + sin(theta) * u_hat + (1 - cos(theta)) * u_hat * u_hat;
}


Vector3d mycommon::LogSO3(const Matrix3d& R) {
	double theta = acos((R.trace() - 1) / 2);
	if (theta < 1e-6) return Vector3d::Zero();
	Matrix3d logR = (R - R.transpose()) / (2 * sin(theta));
	return Vector3d(logR(2, 1), logR(0, 2), logR(1, 0)) * theta;
}
double * mycommon::qt2euc(double a, double b, double c, double d, double t1, double t2, double t3)
{
	double R[9];
	R[0] = 1 - 2 * c * c - 2 * d * d;
	R[1] = 2 * b * c - 2 * a * d;
	R[2] = 2 * b * d + 2 * a * c;
	R[3] = 2 * b * c + 2 * a * d;
	R[4] = 1 - 2 * b * b - 2 * d * d;
	R[5] = 2 * c * d - 2 * a * b;
	R[6] = 2 * b * d - 2 * a * c;
	R[7] = 2 * c * d + 2 * a * b;
	R[8] = 1 - 2 * b * b - 2 * c * c;

	double eu[3];
	rotationMatrixToEulerAngles(R, eu);//kappa phi omega 
	// rotationMatrixToEulerAngles_phi_omega_kappa(R, eu);//phi omega kappa
	// double eu1[3], R1[9];
	// rotationMatrixToEulerAngles_phi_omega_kappa(R, eu1);
	// eulerAnglesToRotationMatrix_phi_omega_kappa(eu1, R1);
	// printf("%f %f %f\n", R[0], R[1], R[2]);
	// printf("%f %f %f\n", R[3], R[4], R[5]);
	// printf("%f %f %f\n", R[6], R[7], R[8]);
	// printf("%f %f %f\n", eu[0], eu[1], eu[2]);
	// printf("%f %f %f\n", R1[0], R1[1], R1[2]);
	// printf("%f %f %f\n", R1[3], R1[4], R1[5]);
	// printf("%f %f %f\n", R1[6], R1[7], R1[8]);
	// printf("%f %f %f\n", eu1[0], eu1[1], eu1[2]);

	double c1, c2, c3;
	c1 = -(R[0] * t1 + R[3] * t2 + R[6] * t3);//t??c
	c2 = -(R[1] * t1 + R[4] * t2 + R[7] * t3);
	c3 = -(R[2] * t1 + R[5] * t2 + R[8] * t3);

	//printf("%f %f %f\n", c1, c2, c3);

	static double euc[6];
	euc[0] = eu[0];
	euc[1] = eu[1];
	euc[2] = eu[2];
	euc[3] = c1;
	euc[4] = c2;
	euc[5] = c3;
	//printf("%f %f %f %f %f %f\n", euc[0], euc[1], euc[2], euc[3], euc[4], euc[5]);
	return euc;
}

void mycommon::add_gano2xyz(vector<pt4>& points, double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist(mean, stddev);

	for (auto& val : points) {
		val.xyz[0] += dist(gen); 
		val.xyz[1] += dist(gen);
		val.xyz[2] += dist(gen);
	}
}
bool mycommon::isRotationMatrix(const Eigen::Matrix3d& R) {
	// 1. ??? R^T * R ?? I
	Eigen::Matrix3d shouldBeIdentity = R.transpose() * R;
	Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
	bool isOrthogonal = shouldBeIdentity.isApprox(I, 1e-6);  // 

	// 2. ???????? ?? 1
	bool isDetOne = std::abs(R.determinant() - 1.0) < 1e-6;

	return isOrthogonal && isDetOne;
}
void mycommon::add_gano2eul(vector<eop8_>& eop, double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist(mean, stddev);
	
	for (auto& val : eop) {
		//??kappa??phi??omega???
		//??phi??omega, kappa???
		double disturb[3], disturbR[9];
		disturb[0] = dist(gen); disturb[1] = dist(gen); disturb[2] = dist(gen);
		// eulerAnglesToRotationMatrix_phi_omega_kappa(disturb, disturbR);//
		eulerAnglesToRotationMatrix(disturb, disturbR);
		Matrix3d R_noise;
		R_noise << disturbR[0], disturbR[1], disturbR[2],
			disturbR[3], disturbR[4], disturbR[5],
			disturbR[6], disturbR[7], disturbR[8];
		//const Vector3d phi(dist(gen), dist(gen), dist(gen));//
		//printf("%f %f %f\n", phi(0), phi(1), phi(2));
		//Matrix3d R_noise = ExpSO3(phi);//??????????
		//printf("%f %f %f\n", R_noise(0, 0), R_noise(0, 1), R_noise(0, 2));
		//printf("%f %f %f\n", R_noise(1, 0), R_noise(1, 1), R_noise(1, 2));
		//printf("%f %f %f\n", R_noise(2, 0), R_noise(2, 1), R_noise(2, 2));
		if (isRotationMatrix(R_noise)) {
			double R[9], eulerAngles[3];
			eulerAngles[0] = val.eop[0];//???????????
			eulerAngles[1] = val.eop[1];
			eulerAngles[2] = val.eop[2];
			// eulerAnglesToRotationMatrix_phi_omega_kappa(eulerAngles, R);
			eulerAnglesToRotationMatrix(eulerAngles, R);
			//printf("%f %f %f\n", eulerAngles[0], eulerAngles[1], eulerAngles[2]);
			//printf("%f %f %f\n", R[0], R[1], R[2]);
			//printf("%f %f %f\n", R[3], R[4], R[5]);
			//printf("%f %f %f\n", R[6], R[7], R[8]);
			Matrix3d R_, R1;
			R_ << R[0], R[1], R[2],
				R[3], R[4], R[5],
				R[6], R[7], R[8];
			//printf("%f %f %f\n", R_(0, 0), R_(0, 1), R_(0, 2));
			//printf("%f %f %f\n", R_(1, 0), R_(1, 1), R_(1, 2));
			//printf("%f %f %f\n", R_(2, 0), R_(2, 1), R_(2, 2));
			R1 = R_noise * R_;
			//printf("%f %f %f\n", R1(0, 0), R1(0, 1), R1(0, 2));
			//printf("%f %f %f\n", R1(1, 0), R1(1, 1), R1(1, 2));
			//printf("%f %f %f\n", R1(2, 0), R1(2, 1), R1(2, 2));
			if (isRotationMatrix(R1)) {
				double R2[9], eua[3];
				R2[0] = R1(0, 0);    R2[1] = R1(0, 1);   R2[2] = R1(0, 2);
				R2[3] = R1(1, 0);    R2[4] = R1(1, 1);   R2[5] = R1(1, 2);
				R2[6] = R1(2, 0);    R2[7] = R1(2, 1);   R2[8] = R1(2, 2);
				//printf("%f %f %f\n", R2[0], R2[1], R2[2]);
				//printf("%f %f %f\n", R2[3], R2[4], R2[5]);
				//printf("%f %f %f\n", R2[6], R2[7], R2[8]);
				// rotationMatrixToEulerAngles_phi_omega_kappa(R2, eua);//R???????kappa??phi??omega??
				rotationMatrixToEulerAngles(R2, eua);
				//printf("%f %f %f\n", eua[0], eua[1], eua[2]);
				val.eop[0] = eua[0];
				val.eop[1] = eua[1];
				val.eop[2] = eua[2];
			}
			else {
				printf("error!");
				val.eop[0] += disturb[0];
				val.eop[1] += disturb[1];
				val.eop[2] += disturb[2];
			}
		}else{
			printf("error!");
			val.eop[0] += disturb[0];
			val.eop[1] += disturb[1];
			val.eop[2] += disturb[2];
		}

		//JacobiSVD<Matrix3d> svd(R1, ComputeFullU | ComputeFullV);
		//R1 = svd.matrixU() * svd.matrixV().transpose();
	}
}
void mycommon::add_gano2c(vector<eop8_>& eop, double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist(mean, stddev);

	for (auto& val : eop) {
		//printf("%f %f %f  ", val.eop[3], val.eop[4], val.eop[5]);
		val.eop[3] += dist(gen);
		val.eop[4] += dist(gen);
		val.eop[5] += dist(gen);
		//printf("%f %f %f\n", val.eop[3], val.eop[4], val.eop[5]);
	}
}
void mycommon::add_gano2obs(vector<obs5>& obs, double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist(mean, stddev);

	for (auto& val : obs) {
		val.uv[0] += dist(gen);
		val.uv[1] += dist(gen);
	}
}
void mycommon::add_gano2cal(vector<cal3>& cal, double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist(mean, stddev);

	for (auto& val : cal) {
		val.f += dist(gen);
		val.cx += dist(gen);
		val.cy += dist(gen);
	}
}
void mycommon::add_gano2cal8(vector<cal8>& cal, double mean, double stddev) {
	random_device rd;
	mt19937 gen(rd());
	normal_distribution<> dist(mean, stddev);

	for (auto& val : cal) {
		val.fx += dist(gen);
		val.fy += dist(gen);
		val.cx += dist(gen);
		val.cy += dist(gen);
		val.k1 += dist(gen);
		val.k2 += dist(gen);
		val.p1 += dist(gen);
		val.p2 += dist(gen);
	}
}
double mycommon::deg2rad(double deg) {
	return deg * M_PI / 180.0;
}
double mycommon::rad2deg(double rad) {
	return rad * 180.0 / M_PI;
}

std::ptrdiff_t mycommon::findViewIndex(const std::vector<eop8_>& views, int targetViewId) {
	// ???????lambda??????viewid
	auto comp = [](const eop8_& item, int viewId) {
		return item.viewid < viewId;
	};


	auto it = std::lower_bound(views.begin(), views.end(), targetViewId, comp);

	if (it != views.end() && it->viewid == targetViewId) {
		std::ptrdiff_t index = std::distance(views.begin(), it);
		return index;
	}
	return -1; 
}

void mycommon::groupAndSort(std::vector<track5>& track) {

	std::sort(track.begin(), track.end(),
		[](const track5& a, const track5& b) {
		if (a.ptid != b.ptid) return a.ptid < b.ptid;
		return a.viewid < b.viewid;
	});
}

void mycommon::groupAndSort(vector<pt4>& pt) {

	sort(pt.begin(), pt.end(),
		[](const pt4& a, const pt4& b) { 
		if (a.ptid != b.ptid) return a.ptid < b.ptid;
		return false;
	});
}