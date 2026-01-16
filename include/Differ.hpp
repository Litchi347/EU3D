// 实现空间离散与重构
// 根据网格中心点的物理量，预测出网格界面（左侧和右侧）上的物理量数值









# include <iostream>
# include <algorithm>





# include "Array.hpp"
# include "Global.hpp"

using namespace std;
using namespace ARRAY;





// 在处理激波（Shock）或物理量剧烈变化的地方，高阶插值会产生非物理的数值振荡（吉布斯现象）
// 限制器的作用是：在平滑区域保持高精度，在剧烈变化区域自动降阶为一阶，以保证数值稳定性
// 限制器 基于 TVD 理论
// 绝不会产生震荡，非常稳定。但会过度平滑流场，导致激波被抹得很宽，分辨率低
double Minmod(double R)
{
    double Minmod = 0.0;
    // 对比前后的斜率，只取其中绝对值最小的一个；如果两个斜率方向相反（一正一负），它直接把斜率设为 0（降阶为一阶）。
    if(R > 0)
        Minmod = min(R, 1.0);
    else
        Minmod = 0.0;
    
    return Minmod;
}
// 通过一个非线性公式在不同斜率间进行调和
double Van_Leer(double R)
{
    double Van_Leer = 0.0;

    Van_Leer = ( R + abs(R)) / (1.0 + abs(R));

    return Van_Leer;
}
// 在满足 TVD 条件的范围内，尽可能选取最大的斜率，以保持界面的尖锐度
double Superbee(double R)
{
    double Superbee = 0.0;

    Superbee = max({0.0,min(2 * R,1.0),min(R,2.0)});

    return Superbee;
}
// 在高梯度区域能够比 Van Leer 提供更好的解析度。在处理包含接触断层的复杂流动时，它往往能提供比 Minmod 更清晰的边界
double Van_Albada(double R)
{
    double Van_Albada = 0.0;

    Van_Albada = (pow(R,2) + R) / (1.0 + pow(R,2));

    return Van_Albada;
}
// 在二阶中心差分和 Minmod 之间寻找平衡
double Double_Minmod(double R)
{
    double Double_Minmod = 0.0;

    if(R > 0)
        Double_Minmod = min({2 * R,1.0,(1 + R) / 2});
    else
        Double_Minmod = 0.0;

    return Double_Minmod; 
}












void Diff_Initial(int direction, Array<double,1> &xnode, Array<double,1> &ynode, Array<double,1> &znode, Array<double,4> &FLR, Array<double,3> &fi, int bc, int limit)
{
}








// 一阶重构 直接取相邻网格中心的值作为界面值
// 左界面值等于当前网格中心值 右界面值等于相邻网格中心值
// fi 是网格中心的原始物理量场
// FLR 是存储结果的数组，FLR(..., 0)是界面左侧，FLR(..., 1)是界面右侧
void MUSCL_1(int direction, Array<double, 1> &xnode, Array<double, 1> &ynode, Array<double, 1> &znode, Array<double, 4> &FLR, Array<double, 3> &fi, int bc, int limit)
{
    int ni = xnode.GetSize() - 2 * bc;
    int nj = ynode.GetSize() - 2 * bc;
    int nk = znode.GetSize() - 2 * bc;

    for(int i = bc - 1;i < ni + bc; i++)
        for(int j = bc - 1;j < nj + bc; j++)
            for(int k = bc - 1;k < nk + bc; k++)
            {
                if(direction == 1)                         // X 方向
                {
                    FLR(i, j, k, 0) = fi(i, j, k);
                    FLR(i, j, k, 1) = fi(i + 1, j, k);
                }
                else if (direction == 2)                   // Y 方向
                {
                    FLR(i, j, k, 0) = fi(i, j, k);
                    FLR(i, j, k, 1) = fi(i, j + 1, k);
                }
                else if (direction == 3)                   // Z 方向
                {
                    FLR(i, j, k, 0) = fi(i, j, k);
                    FLR(i, j, k, 1) = fi(i, j, k + 1);
                }
            }
}










// 高阶重构
// 利用泰勒展开，结合周围多个网格点的梯度预测界面值
void MUSCL_2(int direction, Array<double, 1> &xnode, Array<double, 1> &ynode, Array<double, 1> &znode, Array<double, 4> &FLR, Array<double, 3> &fi, int bc, int limit)
{
    int ni = xnode.GetSize() - 2 * bc;
    int nj = ynode.GetSize() - 2 * bc;
    int nk = znode.GetSize() - 2 * bc;

    double kk = 1.0 / 3.0;
    double A0 = 0.0,A1 = 0.0,A2 = 0.0,L = 0.0,R = 0.0;
    double L1 = 0.0,L2 = 0.0,R1 = 0.0,R2 = 0.0;

    // 数据部分
    for(int i = bc - 1;i < ni + bc;i ++)
        for(int j = bc - 1;j < nj + bc;j ++)
            for(int k = bc - 1;k < nk + bc;k ++)
            {


                if(direction == 1)
                {   // 计算差分
                    A0 = fi(i, j, k) - fi(i - 1, j, k);
                    A1 = fi(i + 1, j, k) - fi(i, j, k);
                    A2 = fi(i + 2, j, k) - fi(i + 1, j, k);
                }
                else if(direction == 2)
                {
                    A0 = fi(i, j, k) - fi(i, j - 1, k);
                    A1 = fi(i, j + 1, k) - fi(i, j, k);
                    A2 = fi(i, j + 2, k) - fi(i, j + 1, k);
                }
                else if(direction == 3)
                {
                    A0 = fi(i, j, k) - fi(i, j, k - 1);
                    A1 = fi(i, j, k + 1) - fi(i, j, k);
                    A2 = fi(i, j, k + 2) - fi(i, j, k + 1);
                }
                // 计算比率
                L = A1 / (A0 + 1e-12);
                R = A1 / (A2 + 1e-12);
                // 选择并应用限制器
                switch(limit)
                {
                case 0:
                    L1 = Minmod(L);
                    L2 = Minmod(1.0 / L);
                    R1 = Minmod(R);
                    R2 = Minmod(1.0 / R);
                    break;
                case 1:
                    L1 = Van_Leer(L);
                    L2 = Van_Leer(1.0 / L);
                    R1 = Van_Leer(R);
                    R2 = Van_Leer(1.0 / R);
                    break;
                case 2:
                    L1 = Van_Albada(L);
                    L2 = Van_Albada(1.0 / L);
                    R1 = Van_Albada(R);
                    R2 = Van_Albada(1.0 / R);
                    break;
                case 3:
                    L1 = Superbee(L);
                    L2 = Superbee(1.0 / L);
                    R1 = Superbee(R);
                    R2 = Superbee(1.0 / R);
                    break;
                case 4:
                    L1 = Double_Minmod(L);
                    L2 = Double_Minmod(1.0 / L);
                    R1 = Double_Minmod(R);
                    R2 = Double_Minmod(1.0 / R);
                    break;
                }
                // 合成界面值
                if(direction == 1)
                {
                    FLR(i, j, k, 0) = fi(i, j, k) + 0.25 * ((1.0 - kk) * L1 + (1 + kk) * L2 * L) * A0;
                    FLR(i, j, k, 1) = fi(i + 1, j, k) - 0.25 * ((1.0 - kk) * R1 + (1 + kk) * R2 * R) * A2;
                }
                else if(direction == 2)
                {
                    FLR(i, j, k, 0) = fi(i, j, k) + 0.25 * ((1.0 - kk) * L1 + (1 + kk) * L2 * L) * A0;
                    FLR(i, j, k, 1) = fi(i, j + 1, k) - 0.25 * ((1.0 - kk) * R1 + (1 + kk) * R2 * R) * A2;
                }
                else if(direction == 3)
                {
                    FLR(i, j, k, 0) = fi(i, j, k) + 0.25 * ((1.0 - kk) * L1 + (1 + kk) * L2 * L) * A0;
                    FLR(i, j, k, 1) = fi(i, j, k + 1) - 0.25 * ((1.0 - kk) * R1 + (1 + kk) * R2 * R) * A2;
                }
            }
}












// Array<double> WENO(int flag,Array<double> U,Array<double> V,Array<double> xnode,Array<double> ynode,Array<double> fi){
//     int ni = xnode.G
// }