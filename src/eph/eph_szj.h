#ifndef EPH_SZJ_H
#define EPH_SZJ_H

#include <string>
#include <vector>

namespace sxwnl {
struct SJ {
    double z;
    double x;
    double s;
    double j;
    double c;
    double h;
    double c2;
    double h2;
    double c3;
    double h3;
    double H0;
    double H;
    double H1;
    double H2;
    double H3;
    double H4;
    std::string sm;
};

struct SJ_S {
    std::string s;
    std::string z;
    std::string j;
    std::string c;
    std::string h;
    std::string ch;
    std::string sj;
    std::string Ms;
    std::string Mz;
    std::string Mj;
};

/**
 * @brief 日月的升中天降,不考虑气温和气压的影响
 */
class SZJ {
public:
    static double L;               //站点地理经度,向东测量为正
    static double fa;              //站点地理纬度
    static std::vector<SJ_S> rts;  //多天的升中降

    static double getH(double h, double w);
    static void Mcoord(double jd, double H0, SJ &r);
    static void Scoord(double jd, int xm, SJ &r);
    static SJ Mt(double jd);
    static SJ Qt(double jd);
    static SJ St(double jd);
    static void calcRTS(double jd, int n, double Jdl, double Wdl, double sq);

private:
    static double E;   //黄赤交角
    static double dt;  //TD-UT
};
}  // namespace sxwnl

#endif