// 三维多组分可压缩化学反应流求解器
// 通过 Reaction React 对象及各参数，实现了对多组分化学反应流的支持
// 不仅能计算单一气体的流动，还能计算涉及多种化学物质的混合物特性



#ifndef FlowField
#define FlowField





# include <iostream>





# include"Array.hpp"
# include"TimeAdv.hpp"
# include"Reaction.hpp"
# include"Function.hpp"

using namespace std;
using namespace ARRAY;

class Flowfield
{
private:

    int ni;
    Array<double,1> xnode;
    double dx;

    int nj;
    Array<double,1> ynode;
    double dy;

    int nk;
    Array<double,1> znode;
    double dz;

    int bc;

    double dt;
    double time;

    double ft1 = 0, ft2 = 0, ft3 = 0, ft4 = 0, ft5 = 0, ft6 = 0;

    int NS;
    int NR;
    Reaction React;
    const double R = 8.31434;
    Array<double,1> Mw;                          // 摩尔质量
    Array<double,1> Ri;
    Array<double,2> Coeff0;
    Array<double,2> Coeff1;

    TimeAdv Time;


    double M_wave;                               // 激波马赫数
    double U_post;                               // 激波前的速度
    double P_pre;                                // 激波前的压力
    double P_post;
    double T_pre;                                // 激波前的温度
    double T_post;

    double cita;                                 // 坐标旋转角度或激波角度
    double V3;
    Array<double,3> U;                           // x 方向上速度
    Array<double,3> V;                           // y 方向上速度
    Array<double,3> W;                           // z 方向上速度
    Array<double,3> P;                           // 静压
    double P_bound;
    Array<double,3> D;                           // 密度
    Array<double,3> T;                           // 温度
    double T_bound;
    Array<double,3> C;                           // 当地声速
    Array<double,3> Ma;                          // 马赫数
    Array<double,3> Wav;
    Array<double,3> Rgas;                        // 混合气体常数
    Array<double,3> Cp;                          // 定压比热
    Array<double,3> H;                           // 单位质量总焓
    Array<double,3> E;                           // 单位质量内能
    Array<double,3> Gamma;                       // 比热比 cp/cv
    Array<double,4> Mr;
    Array<double,4> Mc;
    Array<double,4> Mi;
    Array<double,4> Yi;                          // 组分质量分数
    Array<double,4> Di;                          // 各单一分组的偏密度
    Array<double,1> Cpi;                         // 各单一组分的比热
    Array<double,1> Hi;                          // 各单一组分的焓
    Array<double,1> Ei;                          // 各单一组分的内能

    Array<double,4> F;
    Array<double,4> G;
    Array<double,4> Q;
    Array<double,4> CS;                          // 化学反应源项，代表单位时间内化学反应产生的物质变化率

    Function Fun;

    // 入口边界条件，二维：给定一个方向 X 设置边界条件，是一个 Y-Z 平面
    Array<double,2> Uint, Vint, Wint, Pint, Dint, Tint, Hint, Eint, Gint;
    Array<double,3> Yint;


    Array<double,1> Yi_temp0, YL, YR;
    Array<double,4> PLR, DLR, ULR, VLR, WLR, HLR, GLR, YLR, YL_temp, YR_temp;
    Array<double,3> Yi_temp;
    Array<double,4> Partion_T;
    Array<double,4> RHS;                         // 右端项，存储空间离散后的残差，用于时间迭代


    Array<double,1> m_U_s, m_V_s, m_P_s, m_D_s, m_C_s, m_Gamma_s, m_H_s, m_T_s,m_Yi_s;
    Array<double,1> m_U_r, m_V_r, m_P_r, m_D_r, m_C_r, m_Gamma_r, m_H_r, m_T_r,m_Yi_r;
    Array<double,1> send_data_1, recv_data_1, send_data_2, recv_data_2, send_data_3,recv_data_3;

    struct SendRecv_Data                         // 在 MPI 发送和接收时打包数据
    {
        Array<double,1> D, U, V, P, T, H, Gamma, C, Yi;
        void Initial(int size, int _NS)
        {
            D.Initial(size);
            U.Initial(size);
            V.Initial(size);
            P.Initial(size);
            T.Initial(size);
            H.Initial(size);
            Gamma.Initial(size);
            C.Initial(size);
            Yi.Initial(size * _NS);
        }
    };
    SendRecv_Data send_data, recv_data;

public:
    friend class Euler;


    Flowfield() = default;


    void InputRead(char *initialization, Array<double, 1> &xnode,Array<double, 1> &ynode,Array<double, 1> &znode, int bc);


    void Construction();


    void ReConstruction(int meshnum);


    void FieldInitial(Array<double, 1> &Ri, Array<double, 1> &Mw, Array<double, 4> &Mi_temp, Array<double, 4> &Yi_temp, Array<double, 2> &coeff0, Array<double, 2> &coeff1);
    
    
    void CFLcondition(double cfl, double Final_Time);

    // 多种物理边界的处理方式
    void FieldBoundary_3dODM();
    void FieldBoundary_1dZND();                  // 用于爆轰波模拟的 ZND 模型
    void FieldBoundary_3dRSBI();


    void Advection(int, int);


    void Update_after_Adv();


    void Explicit();
    void Explicit(int step, int i, int j, int k);


    void Update_IMEX(Array<double, 4> &Wi, Array<double, 5> &MD);


    void Update_after_CS();
    void Update_after_CS(int step, int i, int j, int k);


    ~Flowfield(){ ; };


    void AUSM(int, Array<double,4> &, void (*Diff)(int, Array<double,1> &, Array<double,1> &, Array<double,1> &, Array<double,4> &, Array<double,3> &, int ,int));


    double Get_temp(double, int, int, int);
    double Get_temp(double, int, int, int, Array<double, 1> &);


    void GetPartial_T();


    void Mpi_Boundary();



    // 在 MPI 通信前，将多维数组中 Ghost Cells 的数据转换成连续的一维向量，发送到相邻进程后再还原
    void PackagePrev(std::vector<double> &, int i, int j, int k);
    void PackageUpdate(std::vector<double> &, int meshnum);
    void PackageUpdate(std::vector<double> &, int, std::vector<int> &, std::vector<int> &);
    void PackageUpdate(std::vector<double> &, std::vector<int> &, std::vector<int> &);


    void UnpackagePrev(std::vector<double> &, int meshnum);
    void UnpackagePrev(std::vector<double> &, int, std::vector<int> &, std::vector<int> &);
    void UnpackagePrev(std::vector<double> &, std::vector<int> &, std::vector<int> &);
    void UnpackagePrev(std::vector<double> &, Array<int, 1> &, Array<int,1> &);
    void UnpackageUpdate(std::vector<double> &, std::vector<int> &, int, std::vector<int> &, std::vector<int> &);
    void UnpackageUpdate(std::vector<double> &, std::vector<int> &, std::vector<int> &, std::vector<int> &);
};
#endif