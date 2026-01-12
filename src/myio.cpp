#include "myio.h"
//myio::myio(int das)
//{
//
//}
//
//myio::~myio(void)
//{
//}

void myio::ba_readCameraPose(FILE* fp, FILE* fpw, vector<int> deleteid, double* motcam)
{
	int n;
	int lineno = 0;
	double* tofilter = (double*)malloc(6 * sizeof(double));

	while (!feof(fp))
	{
		n = readNDoubles(fp, tofilter, 6);
		if (n == -1)
			break;

		for (int i = 0; i < 6; i++)
			motcam[6 * lineno + i] = tofilter[i];

		auto it = std::find(deleteid.begin(), deleteid.end(), lineno);
		if (it != deleteid.end())
		{
			++lineno;
			continue;
		}
		fprintf(fpw, "%f %f %f %f %f %f\n", tofilter[0], tofilter[1], tofilter[2], tofilter[3], tofilter[4], tofilter[5]);
		++lineno;
	}
	fclose(fp);
	fclose(fpw);
}

vector<int> myio::readf(const char* fname, const char* wname, double* ical, int& ncams)
{
	vector<int> deleteid;
	FILE* fp;
	int  ch = EOF;

	if ((fp = fopen(fname, "r")) == NULL)
	{
		fprintf(stderr, "DataPreprocess: Cannot open calbration file %s, exiting\n", fname);
		return deleteid;
	}

	FILE* fpf_new = fopen(wname, "w");

	double* ptr;
	for (int i = 0; i < ncams; i++)
	{
		ptr = ical + 1 * i;
		int num = fscanf(fp, "%lf", &ptr[0]);
		if (num != 1)
		{
			fprintf(stderr, "DataPreprocess error: Format of Calibration file is wrong");
			return deleteid;
		}
		else
		{
			if (abs(ptr[0]) > 1e+5)//�쳣ֵ
			{
				deleteid.push_back(i);
				continue;
			}
			fprintf(fpf_new, "%f\n", ptr[0]);
		}
	}

	fclose(fp);
	fclose(fpf_new);
	return deleteid;
}
vector<int> myio::readProjection(FILE* fp, FILE* fpw, vector<int> deleteid, double* motz, double* motcam)
{
	vector<int> deleteptid, frn;
	vector<double> pxs;
	int n, nframes, ptno = 0, frameno;
	double ptr[2];
	while (!feof(fp))
	{
		n = readNInts(fp, &nframes, 1);
		if (n != 1)
			break;

		for (int i = 0; i < nframes; ++i)
		{
			n = readNInts(fp, &frameno, 1);
			n += readNDoubles(fp, ptr, 2);

			auto it = std::find(deleteid.begin(), deleteid.end(), frameno);
			if (it != deleteid.end())
				continue;

			//printf("%f %f\n", motz[3 * ptno + 2], motcam[6 * frameno + 5]);
			double dz = motz[3 * ptno + 2] - motcam[6 * frameno + 5];
			if (abs(dz) < 1e-6)
				continue;

			pxs.push_back(ptr[0]);
			pxs.push_back(ptr[1]);


			int maga = 0;
			if (!deleteid.empty())
			{
				
				if (frameno > deleteid[0])
					++maga;
				for (int k = 1; k < deleteid.size(); k++)
				{
					if (frameno > deleteid[k])
						++maga;
					else
						break;
				}
				//frameno -= maga;
			}
			frn.push_back(frameno - maga);
		}
		size_t remdeg = frn.size();
		if (remdeg < 2)
			deleteptid.push_back(ptno);
		else
		{
			fprintf(fpw, "%zd ", remdeg);
			for (int i = 0; i < remdeg; i++)
			{
				fprintf(fpw, "%d %f %f ", frn[i], pxs[i * 2], pxs[i * 2 + 1]);
			}
			fprintf(fpw, "\n");
			//printf("%d ", remdeg);
			//for (int i = 0; i < remdeg; i++)
			//{
			//	printf("%d %f %f ", frn[i], pxs[i * 2], pxs[i * 2 + 1]);
			//}
			//printf("\n");
		}
		frn.clear();
		pxs.clear();
		ptno++;
	}
	fclose(fp);
	fclose(fpw);
	return deleteptid;
}