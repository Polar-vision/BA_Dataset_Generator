#include "myconvert.h"
#include <stdexcept>
#include "triangulation.h"
#include <filesystem> // 需要 C++17 支持

namespace fs = std::filesystem;
// #pragma  comment( lib, "E:/zuo/projects/backup/colmap2pba/opencv/build/x64/vc16/lib/opencv_world4110.lib")

// 返回动态分配的数组，调用者需负责 delete[]
double* argvsparse(const string& str) {
    const double default_values[] = { 0.0, 0.0 }; // 默认值
    double* result = new double[2];
    result[0] = default_values[0];
    result[1] = default_values[1];

    try {
        size_t start = str.find('[');
        size_t end = str.find(']');

        if (start == string::npos || end == string::npos || start >= end) {
            return result; // 格式错误，返回默认值
        }

        string values_str = str.substr(start + 1, end - start - 1);
        size_t comma_pos = values_str.find(',');

        if (comma_pos == string::npos) {
            return result; // 缺少逗号，返回默认值
        }

        string val1 = values_str.substr(0, comma_pos);
        string val2 = values_str.substr(comma_pos + 1);

        // 去除首尾空格
        auto trim = [](string& s) {
            s.erase(0, s.find_first_not_of(" \t"));
            s.erase(s.find_last_not_of(" \t") + 1);
        };

        trim(val1);
        trim(val2);

        // 转换为 double
        result[0] = stod(val1);
        result[1] = stod(val2);

    }
    catch (const exception& e) {
        // 捕获 stod 可能的异常（无效数字或越界）
        // 保持默认值
        std::cerr << "Exception: " << e.what() << std::endl;
    }

    return result;
}
typedef enum Operate_Type
{
    ini = 1,//解析colmap文件+添加白噪声+重新三角化
    abseva = 2,//解析colmap文件
}OperateType;
int main(int argc, char* argv[])
{
    printf("%s\n", "你好colmap2pba!");
        // 方法 1：直接转换为绝对路径（基于当前工作目录）
    fs::path abs_path = fs::absolute(dt);
    std::cout << "绝对路径: " << abs_path << std::endl;
    OperateType otype = ini;
    // OperateType otype = abseva;
    string pP, obsstd, eulstd, posstd, xyzstd, calstd;
    string cNa, xNa;

    string originalPath(dt);
    size_t pos = originalPath.find_last_of("/\\");
    string parentPath = (pos != string::npos) ? originalPath.substr(0, pos) : "";
    string P = parentPath + "/";

    argv[1] = const_cast<char*>(P.data());
    pP = argv[1];

    for(int i=0;i<2;i++){
        if(i==0){otype=ini;}
        if(i==1){otype=abseva;}
        if (otype == ini) {
            argv[2] = const_cast<char*>("--obs-std=[0.0, 0.0]");
            argv[3] = const_cast<char*>("--eul-std=[0.0, 2.0]");//以度为单位
            argv[4] = const_cast<char*>("--pos-std=[0.0, 0.015]");//
            argv[5] = const_cast<char*>("--xyz-std=[0.0, 0.0]");
            argv[6] = const_cast<char*>("--cal-std=[0.0, 0.0]");
            
            obsstd = argv[2];
            eulstd = argv[3];
            posstd = argv[4];
            xyzstd = argv[5];
            calstd = argv[6];

            cNa = "Cam.txt";
            xNa = "XYZ.txt";
        }
        else if (otype == abseva) {
            /*绝对精度评价专用*/
            // pP = const_cast<char*>(argv[1]);
            // obsstd = const_cast<char*>(argv[2]);
            // eulstd = const_cast<char*>(argv[3]);
            // posstd = const_cast<char*>(argv[4]);
            // xyzstd = const_cast<char*>(argv[5]);
            // calstd = const_cast<char*>(argv[6]);
            argv[2] = const_cast<char*>("--obs-std=[0.0, 0.0]");
            argv[3] = const_cast<char*>("--eul-std=[0.0, 0.0]");//以度为单位
            argv[4] = const_cast<char*>("--pos-std=[0.0, 0.0]");//
            argv[5] = const_cast<char*>("--xyz-std=[0.0, 0.0]");
            argv[6] = const_cast<char*>("--cal-std=[0.0, 0.0]");
            // pP = argv[1];
            obsstd = argv[2];
            eulstd = argv[3];
            posstd = argv[4];
            xyzstd = argv[5];
            calstd = argv[6];
            cNa = "G-Cam.txt";
            xNa = "G-XYZ.txt";
            /*绝对精度评价专用*/
        }


        string p1r = pP + "points3D.txt";
        string p2r = pP + "images.txt";
        string p3r = pP + "cameras.txt";
        string p1w = pP + xNa;// "G-XYZ.txt";//G-XYZ.txt
        string p2w = pP + cNa;// "G-Cam.txt";//G-Cam.txt
        string p3w = pP + "cal.txt";
        string p4w = pP + "Feature.txt";
        //string p2w = pP + "Cam_new.txt";
        //string p2w = pP + "SBA-LM-FinalPose.txt";
        //string p3w = pP + "cal_new.txt";
        //string p4w = pP + "Feature_new.txt";

        double* sigobs, * sigeul, * sigpos, * sigxyz, * sigcal;
        sigobs = argvsparse(obsstd);
        sigeul = argvsparse(eulstd);
        sigpos = argvsparse(posstd);
        sigxyz = argvsparse(xyzstd);
        sigcal = argvsparse(calstd);

        myconvert maga(1);
        //maga.bal2pba();
        //maga.outlierdet();

        maga.colmap2pba(p1r, p2r, p3r, p1w, p2w, p3w, p4w, \
            sigobs, sigeul, sigpos, sigxyz, sigcal);
        
        delete[]sigobs; sigobs = nullptr;
        delete[]sigeul; sigeul = nullptr;
        delete[]sigpos; sigpos = nullptr;
        delete[]sigxyz; sigxyz = nullptr;
        delete[]sigcal; sigcal = nullptr;

        if (otype == ini) {
            triangulation tr(1);
            tr.nc = maga.nc;
            tr.ba_initialize(const_cast<char*>(p2w.c_str()), \
                const_cast<char*>(p4w.c_str()), \
                const_cast<char*>(p3w.c_str()), \
                NULL);
        }
        else if (otype == abseva) {

        }
    }
}
