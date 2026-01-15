// 计算 CFD 或燃烧化学动力学
// 处理化学反应速率、热力学属性及化学源项的数值求解



#ifndef REACTION_H
#define REACTION_H





# include <fstream>
# include <iostream>
# include <string>
# include <vector>





# include "Array.hpp"
using namespace ARRAY;

class Reaction
{
private:

    int ni;                                                                                                                                                    // X 方向的网格节点数
    Array<double,1> xnode;                                                                                                                                     // 存储网格节点坐标的 1D 数组

    int nj;
    Array<double,1> ynode;

    int nk;
    Array<double,1> znode;

    int bc;                                                                                                                                                    // 边界条件的宽度 Ghost Cells

    int NS;                                                                                                                                                    // 化学物种的总数
    int NR;                                                                                                                                                    // 化学反应步的总数
    int TB;

    Array<double,2> Stoi_F;                                                                                                                                    // 正向反应的化学计量系数矩阵
    Array<double,2> Stoi_B;                                                                                                                                    // 逆向反应的化学计量系数矩阵

    // 阿伦尼乌斯公式参数
    Array<double,1> Af;                // 指前因子
    Array<double,1> Bf;                // 温度指数
    Array<double,1> Eaf;               // 活化能

    Array<double,2> React_TB;
    // 存储 NASA 多项式系数
    Array<double,2> Coeff0;
    Array<double,2> Coeff1;

    const double R = 8.31434;
    const double Ru = 1.987;
    const double P0 = 101325;

    Array<double,1> Mw;                // 各物种的分子量
    Array<double,1> Ri;                // 各物种的气体常数

    Array<double,4> Mc;                // 摩尔浓度
    Array<double,4> Mr;                // 摩尔比
    Array<double,1> Mr_temp;
    Array<double,4> Mi;                // 摩尔分数
    Array<double,4> Yi;                // 质量分数
    Array<double,4> Di;                // 密度与质量分数的乘积，即偏密度

    double T_local;

    Array<double,1> Hi;                // 各物种的焓
    Array<double,1> Si;                // 各物种的熵
    Array<double,1> Gi;                // 各物种的吉布斯自由能

    Array<double,1> KF;                // 正向反应速率常数
    Array<double,1> KB;                // 逆向反应速率常数
    Array<double,1> Kp;                // 吉布斯自由能计算平衡常数
    Array<double,1> Kc;                // 吉布斯自由能计算平衡常数

    Array<double,1> RR_F;              // 正向反应速率
    Array<double,1> RR_B;              // 逆向反应速率
    Array<double,1> R_TB;              // 三体效应，考虑某些反应需要第三种媒介分子参与
    Array<double,1> RR;                // 净反应速率

    Array<double,1> Wi;                // 每个物种的净生成率
    Array<double,4> CMS;               // 化学反应产生的质量源项
    Array<double,2> WJH1,WJH2;         // 对浓度的导数项和对温度的导数项
    Array<double,5> MD;
    
    Array<double,1> P,Q;

    Array<int,3> Nchem;                // 存储每个网格点为了保证计算稳定所需的“化学子步”数 因为化学反应通常比流体流动快得多（刚性问题），需要更小的时间步长
    int NchemNow = 1;
    std::vector<int> NchemMax_Rank;    // 当前进程在每一时间步中，网格内最大的亚步值
    Array<int,1> NchemMax_Total;       // 全局汇总的缓冲区
    std::vector<int> NchemMax;         // 全局在每一个时间步中，绝对最大的亚步值
    int NchemSum;                      // 用于 MPI 多进程间的负载统计和均衡
    Array<int,1> NchemTotal;
Array<int,1> NchemTotalBalance;;
    struct intermidatePara
    {
        double P[16];
        double Q[16];
        double Hi[16];
        double Si[16];
        double Gi[16];
        double Wi[16];
        double RR_F[24];
        double RR_B[24];
        double R_TB[24];
        double KF[24];
        double Kp[24];
        double Kc[24];
        double KB[24];
        double RR[24];
    };
    struct interParaDiag
    {
        double P[16];
        double Q[16];
        double Hi[16];
        double Si[16];
        double Gi[16];
        double Wi[16];
        double RR_F[24];
        double RR_B[24];
        double R_TB[24];
        double KF[24];
        double Kp[24];
        double Kc[24];
        double KB[24];
        double RR[24];
        double WJH1[12][24] = {};
        double WJH2[12][24] = {};
    };

public:
    friend class Euler;


    Reaction() = default;

    // 从外部文件读取化学反应机理和热力学数据库 解析物种数、反应步及 NASA 系数
    void ReactionRead(char *reaction_model,char *thermofile,Array<double,1> &xnode,Array<double,1> &ynode,Array<double,1> &znode,int bc);

    // 根据读取到的 NS 和 NR 动态初始化（分配空间）所有的 Array 数组
    void ReactionConstruction();

    void ReConstruction(int meshnum);

    // 设置初始流场状态
    void ReactionInitial();

    // 求解化学刚性常微分方程 ODE
    Array<double,4> Trapezoid(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 3> &T, double dt);
    void Trapezoid(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 4> &Yi_temp, Array<double, 3> &T, double dtm, int step, int i, int j, int k);
    void TrapezoidPrection(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 4> &Yi_temp, Array<double, 3> &T, double dtm);

    // 创建化学雅可比矩阵的对角线部分
    void Diagonalized(Array<double,4> &Mc_temp, Array<double,4> &Di_temp,Array<double,3> &T,Array<double,4> &Partial_T);




    ~Reaction(){ ; };

    // 热力学属性计算：输入温度 T，根据内部存储的 Coeff 系数，利用多项式计算定压比热、焓和熵
    double GetCpi(double T, double R, int SP, Array<double,2> &, Array<double,2> &);


    double GetHi(double T, double R, int SP, Array<double,2> &, Array<double,2> &);


    double GetSi(double T, double R, int SP, Array<double,2> &, Array<double,2> &);


    Array<double,4> GetWi();           // 获得化学源项


    Array<double,5> GetMD();           // 获得雅可比矩阵


    void PushNchemMax();

    void GetNchemMax();
    // 并行数据处理
    void PackagePrev(std::vector<int> &, int i, int j, int k);
    void UnpackagePrev(std::vector<int> &, int, std::vector<int> &, std::vector<int> &);
    void UnpackagePrev(std::vector<int> &, std::vector<int> &, std::vector<int> &);
    void PackageUpdate(std::vector<int> &, int, std::vector<int> &, std::vector<int> &);
    void UnpackageUpdate(std::vector<int> &, std::vector<int> &, int, std::vector<int> &, std::vector<int> &);
};
#endif