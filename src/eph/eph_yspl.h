#ifndef EPH_YSPL_H
#define EPH_YSPL_H

#include <array>
#include <string>

namespace sxwnl {
struct RE0 {
    double e_mRad;
    double eShadow;
    double eShadow2;
    double x;
    double y;
    double mr;
    double er;
    double Er;
    double t;
};

class YS_PL {  //月食快速计算器
public:
    static std::array<double, 7> lT;
    static std::string LX;
    static double sf;
    static void lecXY(double jd, RE0 &re);
    static double lineT(const RE0 &G, double v, double u, double r, bool n);
    static void lecMax(double jd);
};
}  // namespace sxwnl

#endif