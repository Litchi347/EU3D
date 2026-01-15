







#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <iomanip>
#include <mpi.h>
#include <omp.h>
#include <cmath>
#include <numeric>





#include "Reaction.hpp"
#include "FileReader.hpp"
#include "Function.hpp"
#include "Global.hpp"










using namespace std;










// 初始化 -> 分配空间
void Reaction::ReactionConstruction()                                     // 在程序初始化阶段调用，分配整个计算区域所需的全部数组。全网格
{
    Stoi_F.Initial(NS, NR);
    Stoi_B.Initial(NS, NR);

    Af.Initial(NR);
    Bf.Initial(NR);
    Eaf.Initial(NR);

    React_TB.Initial(NS, NR);

    Mr_temp.Initial(NS);
    Mr.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Mi.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Yi.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);

    Coeff0.Initial(9, NS);
    Coeff1.Initial(9, NS);
    Mw.Initial(NS);
    Ri.Initial(NS);

    Hi.Initial(NS);
    Si.Initial(NS);
    Gi.Initial(NS);

    KF.Initial(NR);
    KB.Initial(NR);
    Kp.Initial(NR);
    Kc.Initial(NR);

    RR_F.Initial(NR);
    RR_B.Initial(NR);
    R_TB.Initial(NR);
    RR.Initial(NR);
    Mc.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Di.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Wi.Initial(NS);
    CMS.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);

    WJH1.Initial(NS, NR);
    WJH2.Initial(NS, NR);
    MD.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 3, NS + 3);    // NS + 3 为了存储雅可比矩阵

    P.Initial(NS);
    Q.Initial(NS);

    Nchem.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    Nchem.Fill(1);
    NchemTotal.Initial(numprocs);
    NchemTotalBalance.Initial(numprocs);
}
void Reaction::ReConstruction(int meshnum)                                // 当需要将一部分网格点抽理出来单独计算，或在 MPI 时，将多维数据“拍扁”成一维序列进行快速操作
{
    Mr.Initial(meshnum, 1, 1, NS);
    Mi.Initial(meshnum, 1, 1, NS);
    Yi.Initial(meshnum, 1, 1, NS);
    Mc.Initial(meshnum, 1, 1, NS);
    Di.Initial(meshnum, 1, 1, NS);
    CMS.Initial(meshnum, 1, 1, NS);
    MD.Initial(meshnum, 1, 1, NS + 3, NS + 3);
    Nchem.Initial(meshnum, 1, 1);
    Nchem.Fill(1);
}









// 将文件输入流的读取指针 in 定位到指定的行号 line
ifstream &seek_to_line(ifstream &in, int line)
{
    int i;
    char buf[1024];

    // 文件指针重置到文件的最开头
    in.seekg(0, ios::beg);

    // 逐行读取当前行并存入临时缓冲区 buf，不被处理，为了跳过这些行
    for (i = 0; i < line; i++)
    {
        in.getline(buf, sizeof(buf));
    }
    // 循环结束后停在第 line 行
    return in;
}










// 整个类的数据入口：从外部文本文件中读取并解析“化学动力学机理”和“热力学参数”，并据此初始化程序内存
// 这次模拟涉及哪些化学物质、它们如何反应、以及在不同温度下它们的物理性质如何
void Reaction::ReactionRead(char *reaction_model, char *thermofile, Array<double, 1> &xnode, Array<double, 1> &ynode, Array<double, 1> &znode, int bc)
{
    this->xnode = xnode;
    this->ynode = ynode;
    this->znode = znode;
    this->bc = bc;

    ni = xnode.GetSize() - 2 * bc;
    nj = ynode.GetSize() - 2 * bc;
    nk = znode.GetSize() - 2 * bc;

    ifstream fin;
    fin.open(reaction_model);                                             // 跳转到指定行，读取规模参数、内存分配、化学计量参数、阿伦尼乌斯参数、三体效率等信息
    int a[100] = {0}, i = 0, line = 2;
    string s, item;
    seek_to_line(fin, line);                                              // 打开 reaction_model 的第2行，为每个物种读取两段温度范围内的系数，用于计算比热、焓和熵
    getline(fin, s, '\n');

    // 将读取到的一行字符串 s 放进 text_stream 中，
    stringstream text_stream(s);


    while (std::getline(text_stream, item, '\t'))
    {
        a[i] = stoi(item);
        i++;
    }
    NS = a[0];
    NR = a[1];
    TB = a[2];




    ReactionConstruction();


    line = line + 6;
    seek_to_line(fin, line);

    for (int j = 0; j < NR; j++)
    {
        i = 0;
        getline(fin, s, '\n');
        stringstream reaction_forward(s);
        while (std::getline(reaction_forward, item, '\t'))
        {
            Stoi_F(i, j) = stoi(item);

            i++;
        }

        i = 0;
        getline(fin, s, '\n');
        stringstream reaction_backward(s);
        while (std::getline(reaction_backward, item, '\t'))
        {
            Stoi_B(i, j) = stoi(item);

            i++;
        }

    }



    line = line + 2 * NR + 2;
    seek_to_line(fin, line);

    for (int j = 0; j < NR; j++)
    {
        i = 0;
        double b[4] = {0};
        while (fin >> s)
        {
            stringstream geek(s);
            geek >> b[i];
            i++;
            if (i == 4)
                break;
        }
        Af(j) = b[1];
        Bf(j) = b[2];
        Eaf(j) = b[3];


    }



    line = line + NR + 3;
    seek_to_line(fin, line);
    React_TB.Fill(0.0);
    for (int j = 0; j < TB; j++)
    {
        i = 0;
        while (fin >> s)
        {
            stringstream geek(s);
            geek >> a[i];
            i++;
            if (i == NS + 1)
                break;
        }

        for (int k = 0; k < NS; k++)
        {
            React_TB(k, a[0] - 1) = a[k + 1];

        }

    }



    line = line + TB + 2;
    seek_to_line(fin, line);
    i = 0;
    double c[100] = {0};
    while (fin >> s)
    {
        stringstream geek(s);
        geek >> c[i];
        i++;
        if (i == NS)
            break;
    }
    for (int i = 0; i < NS; i++)
    {
        Mr_temp(i) = c[i];
        cout << Mr_temp(i) << '\t';
    }
    cout << endl;
    fin.close();

    ifstream in;
    in.open(thermofile);


    line = 15;
    double aa[100] = {0};

    for (int j = 0; j < NS; j++)
    {
        seek_to_line(in, line);
        i = 0;
        while (in >> s)
        {
            stringstream geek(s);
            geek >> aa[i];

            i++;
            if (i == 19)
                break;
        }
        for (int k = 0; k < 9; k++)
        {
            Coeff0(k, j) = aa[k];
        
        }

        for (int k = 0; k < 9; k++)
        {
            Coeff1(k, j) = aa[k + 9];
        
        }

        Mw(j) = aa[18];

        line = line + 6;
    }
    fin.close();

    ReactionInitial();

    return;
}











// 在三维空间中初始化一个气泡及其周围介质的化学组分分布
void Reaction::ReactionInitial()
{
    Function Fun;
    Yi.Fill(0.0);
    Mi.Fill(0.0);
    Mr.Fill(0.0);



    double Bubble_CenetrX = 0.038;                                        // 气泡位置
    double Bubble_CenetrY = 0.0;
    double Bubble_Ir = 0.022;                                             // 气泡内部，填充了另一种物质
    double Bubble_Or = 1.2 * Bubble_Ir;                                   // 外部是空气
    double alpha = 5;                                                     // 控制界面模糊程度的系数（值越大，界面越陡峭）
    for (int s = 0; s < NS; s++)
        Ri(s) = R / Mw(s);
    for (int i = 0; i < ni + 2 * bc; i++)
        for(int j = 0; j < nj + 2 * bc; j++)
            for(int k = 0; k < nk + 2 * bc; k++)
            {
                double Mr_SUM = 0.0;
                double Total_Mw = 0.0;
                if (pow(xnode(i) - Bubble_CenetrX, 2) + pow(ynode(j) - Bubble_CenetrY, 2) > pow(Bubble_Or, 2))
                {
                    Mr(i, j, k, 4) = 1.0;                                 // 索引4（氧气） 占一份
                    Mr(i, j, k, 8) = 3.76;                                // 索引8 （氮气）占一份，空气比例
                    for(int s = 0; s < NS; s++)
                    {
                        Mr_SUM += Mr(i, j, k, s);                         // 总摩尔数
                        Total_Mw += Mr(i, j, k, s) * Mw(s);               // 总质量
                    }
                    for (int s = 0; s < NS; s++)
                    {
                        Mi(i, j, k, s) = Mr(i, j, k, s) / Mr_SUM;
                        Yi(i, j, k, s) = Mr(i, j, k, s) * Mw(s) / Total_Mw;
                    }
                }
                else if (pow(xnode(i) - Bubble_CenetrX, 2) + pow(ynode(j) - Bubble_CenetrY, 2) < pow(Bubble_Ir, 2))
                {
                    Mr(i, j, k, 5) = 1.0;                                 // 内部被某燃料气体占据
                    for (int s = 0; s < NS; s++)
                    {
                        Mr_SUM += Mr(i, j, k, s);
                        Total_Mw += Mr(i, j, k, s) * Mw(s);
                    }
                    for (int s = 0; s < NS; s++)
                    {
                        Mi(i, j, k, s) = Mr(i, j, k, s) / Mr_SUM;
                        Yi(i, j, k, s) = Mr(i, j, k, s) * Mw(s) / Total_Mw;
                    }
                }
                else                                                      // 界面过渡区：利用高斯函数创建了一个平滑的过渡
                {
                    double r = sqrt(pow(xnode(i) - Bubble_CenetrX, 2) + pow(ynode(j) - Bubble_CenetrY, 2));
                    Yi(i, j, k, 5) = exp(-alpha * pow((r - Bubble_Ir) / (Bubble_Or - Bubble_Ir), 2));
                    Yi(i, j, k, 4) = (1 - Yi(i, j, k, 5)) * Mw(4) / (Mw(4) + 3.76 * Mw(8));
                    Yi(i, j, k, 8) = (1 - Yi(i, j, k, 5)) * 3.76 * Mw(8) / (Mw(4) + 3.76 * Mw(8));
                }
            }
}












// 采用半隐式梯形法求解化学动力学常微分方程，计算一个时间步内各物种浓度的变化
Array<double, 4> Reaction::Trapezoid(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 3> &T, double dt)
{

    // 遍历每一个网格点
    for (int i = bc; i < ni + bc; i++)
        for(int j = bc; j < nj + bc; j++)
            for(int k = bc; k < nk + bc; k++)
            {
                T_local = T(i, j, k);
                for (int s = 0; s < NS; s++)
                {
                    Di(i, j, k, s) = Di_temp(i, j, k, s) * 1e-3;          // 转化为标准单位
                    Mc(i, j, k, s) = Mc_temp(i, j, k, s) * 1e-6;
                    Hi(s) = GetHi(T_local, R, s, Coeff0, Coeff1);         // 计算焓和熵，得到吉布斯自由能
                    Si(s) = GetSi(T_local, R, s, Coeff0, Coeff1);
                    Gi(s) = Hi(s) - Si(s) * T_local;
                }

                // 化学反应速率的计算
                RR_F.Fill(1.0);                                           // 可能同时发生 NR 个反应，综合考虑所有包含该物种的反应
                RR_B.Fill(1.0);
                R_TB.Fill(0);
                for (int r = 0; r < NR; r++)
                {
                    double Xr = 0.0, delta_g = 0.0;
                    KF(r) = Af(r) * pow(T_local, Bf(r)) * exp(-Eaf(r) / (Ru * T_local));
                    for (int s = 0; s < NS; s++)
                    {
                        Xr += (Stoi_B(s, r) - Stoi_F(s, r));
                        delta_g += (Stoi_B(s, r) - Stoi_F(s, r)) * Gi(s);
                        RR_F(r) *= pow(Mc(i, j, k, s), Stoi_F(s, r));
                        RR_B(r) *= pow(Mc(i, j, k, s), Stoi_B(s, r));
                        R_TB(r) += Mc(i, j, k, s) * React_TB(s, r);
                    }
                    Kp(r) = exp(-delta_g / (R * T_local)) * pow((P0 * 1e-6), Xr);
                    Kc(r) = Kp(r) * pow(R * T_local, -Xr);
                    KB(r) = KF(r) / Kc(r);
                    if (R_TB(r) == 0.0)
                        R_TB(r) == 1.0;
                    RR(r) = KF(r) * RR_F(r) - KB(r) * RR_B(r);
                }

                // 构建化学源项
                for (int s = 0; s < NS; s++)
                {
                    double temp = 0.0;
                    for (int r = 0; r < NR; r++)
                        temp += (Stoi_B(s, r) - Stoi_F(s, r)) * R_TB(r) * RR(r);
                    Wi(s) = Mw(s) * temp;
                    CMS(i, j, k, s) = Wi(s) * 1e3;
                }

                P.Fill(0), Q.Fill(0);                                     // P 消耗项， Q 生产项
                for (int s = 0; s < NS; s++)                              // 点隐式处理强刚性问题
                {
                    for (int r = 0; r < NR; r++)
                    {
                        Q(s) += (Stoi_B(s, r) * KF(r) * RR_F(r) + Stoi_F(s, r) * KB(r) * RR_B(r) * R_TB(r));
                        P(s) += (Stoi_B(s, r) * KB(r) * RR_B(r) + Stoi_F(s, r) * KF(r) * RR_F(r) * R_TB(r));
                    }
                    Q(s) = Mw(s) * Q(s);
                    if (Mc(i, j, k, s) == 0)
                    {
                        P(s) = 0.0;
                    }
                    else
                    {
                        P(s) = P(s) / Mc(i, j, k, s);
                    }
                    Di(i, j, k, s) = ((1.0 - dt / 2.0 * P(s)) * Di(i, j, k, s) + dt * Q(s)) / (1.0 + dt / 2.0 * P(s)) * 1e3;
                }
            }
    return Di;
}

// 针对单个网格点进行化学刚性方程的亚布迭代求解 -> 传入特定的 i, j, k 坐标
// 每个线程负责一个点，独立完成该点的热力学性质计算、速率计算和梯形格式更新
void Reaction::Trapezoid(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 4> &Yi_temp, Array<double, 3> &T, double dtm, int step, int i, int j, int k)
{
    double dt = 0.0, dt_temp = 0.0;
    double T_local;
    intermidatePara __attribute__((aligned(64))) ip;                      // 为临时变量结构体 ip 分配64字节对齐的内存


    // 支持化学自适应步长，防止数值爆炸
    dt = dtm / Nchem(i, j, k);

    T_local = T(i, j, k);
    for (int s = 0; s < NS; s++)
    {
        Di(i, j, k, s) = Di_temp(i, j, k, s) * 1e-3;
        Mc(i, j, k, s) = Mc_temp(i, j, k, s) * 1e-6;
        ip.Hi[s] = GetHi(T_local, R, s, Coeff0, Coeff1);
        ip.Si[s] = GetSi(T_local, R, s, Coeff0, Coeff1);
        ip.Gi[s] = ip.Hi[s] - ip.Si[s] * T_local;
    }

    // 只存储当前这一个点涉及到的所有反应中间值
    std::fill(ip.RR_F, ip.RR_F + NR, 1.0);                                // 指针的操作，与 Fill 作用相同
    std::fill(ip.RR_B, ip.RR_B + NR, 1.0);
    std::fill(ip.R_TB, ip.R_TB + NR, 1.0);
    for (int r = 0; r < NR; r++)
    {
        double Xr = 0.0, delta_g = 0.0;
        ip.KF[r] = Af(r) * pow(T_local, Bf(r)) * exp(-Eaf(r) / (Ru * T_local));
        for (int s = 0; s < NS; s++)
        {
            Xr += (Stoi_B(s, r) - Stoi_F(s, r));
            delta_g += (Stoi_B(s, r) - Stoi_F(s, r));
            ip.RR_F[r] *= pow(Mc(i, j, k, s), Stoi_F(s, r));
            ip.RR_B[r] *= pow(Mc(i, j, k, s), Stoi_B(s, r));
            ip.R_TB[r] += Mc(i, j, k, s) * React_TB(s, r);
        }
        ip.Kp[r] = exp(-delta_g / (R * T_local)) * pow((P0 * 1e-6), Xr);
        ip.Kc[r] = ip.Kp[r] * pow(R * T_local, -Xr);
        ip.KB[r] = ip.KF[r] / ip.Kc[r];
        if (ip.R_TB[r] = 0.0)
            ip.R_TB[r] = 1.0;
        ip.RR[r] = ip.KF[r] * ip.RR_F[r] - ip.KB[r] * ip.RR_B[r];
    }

    // 构建化学源项
    for (int s = 0; s < NS; s++)
    {
        double temp = 0.0;
        for (int r = 0; r < NR; r++)
            temp += (Stoi_B(s, r) - Stoi_F(s, r)) * ip.R_TB[r] * ip.RR[r];
        ip.Wi[s] = Mw(s) * temp;
        CMS(i, j, k, s) = ip.Wi[s] * 1e3;
    }


























    
    std::fill(ip.P, ip.P + NS, 0.0);
    std::fill(ip.Q, ip.Q + NS, 0.0);
    for (int s = 0; s < NS; s++)
    {
        for (int r = 0; r < NR; r++)
        {
            ip.Q[s] += (Stoi_B(s, r) * ip.KF[r] * ip.RR_F[r] + Stoi_F(s, r) * ip.KB[r] * ip.RR_B[r] * ip.R_TB[r]);
            ip.P[s] += (Stoi_B(s, r) * ip.KB[r] * ip.RR_B[r] + Stoi_F(s, r) * ip.KF[r] * ip.RR_F[r] * ip.R_TB[r]);
        }
        ip.Q[s] = Mw(s) * ip.Q[s];
        if (Mc(i, j, k, s) == 0)
        {
            ip.P[s] = 0.0;
        }
        else
        {
            ip.P[s] = ip.P[s] / Mc(i, j, k, s);
        }
        Di_temp(i, j, k, s) = ((1.0 - dt / 2.0 * ip.P[s]) * Di(i, j, k, s) + dt * ip.Q[s]) / (1.0 + dt / 2.0 * ip.P[s]) * 1e3;
    }

    // return Di;
}
// 预估化学反应的时间尺度，并据此计算每个网格点所需的亚步迭代次数
void Reaction::TrapezoidPrection(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 4> &Yi_temp, Array<double, 3> &T, double dtm)
{

    for (int i = bc; i < ni + bc; i++)
        for (int j = bc; j < nj + bc; j++)
            for (int k = bc; k < nk + bc; k++)
            {
                double dt = dtm, dt_temp = 0.0;
                double T_local = T(i, j, k);
                intermidatePara __attribute__((aligned(64))) ip;

                for (int s = 0; s < NS; s++)
                {
                    Di(i, j, k, s) = Di_temp(i, j, k, s) * 1e-3;
                    Mc(i, j, k, s) = Mc_temp(i, j, k, s) * 1e-6;
                    ip.Hi[s] = GetHi(T_local, R, s, Coeff0, Coeff1);
                    ip.Si[s] = GetSi(T_local, R, s, Coeff0, Coeff1);
                    ip.Gi[s] = ip.Hi[s] - ip.Si[s] * T_local;
                }

                // 计算瞬时反应速率
                std::fill(ip.RR_F, ip.RR_F + NR, 1.0);
                std::fill(ip.RR_B, ip.RR_B + NR, 1.0);
                std::fill(ip.R_TB, ip.R_TB + NR, 1.0);
                for (int r = 0; r < NR; r++)
                {
                    double Xr = 0.0, delta_g = 0.0;
                    ip.KF[r] = Af(r) * pow(T_local, Bf(r)) * exp(-Eaf(r) / (Ru * T_local));
                    for (int s = 0; s < NS; s++)
                    {
                        Xr += (Stoi_B(s, r) - Stoi_F(s, r));
                        delta_g += (Stoi_B(s, r) - Stoi_F(s, r));
                        ip.RR_F[r] *= pow(Mc(i, j, k, s), Stoi_F(s, r));
                        ip.RR_B[r] *= pow(Mc(i, j, k, s), Stoi_B(s, r));
                        ip.R_TB[r] += Mc(i, j, k, s) * React_TB(s, r);
                        // cout << React_TB(s, r) << " ";
                    }
                    // cout << endl;
                    ip.Kp[r] = exp(-delta_g / (R * T_local)) * pow((P0 * 1e-6), Xr);
                    ip.Kc[r] = ip.Kp[r] * pow(R * T_local, -Xr);
                    ip.KB[r] = ip.KF[r] / ip.Kc[r];
                    if (ip.R_TB[r] = 0.0)
                        ip.R_TB[r] = 1.0;
                    ip.RR[r] = ip.KF[r] * ip.RR_F[r] - ip.KB[r] * ip.RR_B[r];
                }


                for (int s = 0; s < NS; s++)
                {
                    double temp = 0.0;
                    for (int r = 0; r < NR; r++)
                        temp += (Stoi_B(s, r) - Stoi_F(s, r)) * ip.R_TB[r] * ip.RR[r];
                    ip.Wi[s] = Mw(s) * temp;
                    // CMS(i, j, k, s) = ip.Wi[s] * 1e3;
                }


                for (int s = 0; s < NS; s++)
                {
                    if (Yi_temp(i, j, k, s) >= 1e-6)
                    {
                        if (ip.Wi[s] != 0)
                        {
                            dt_temp = abs(-Di(i, j, k, s) / ip.Wi[s]);    // 时间 = 当前密度 / 反应速率
                            if (dt_temp < dt)                             // 寻找最小时间步长
                                dt = dt_temp;
                        }
                    }
                }
                Nchem(i, j, k) = int(ceil(dtm / dt));                     // 确定循环次数

                if (Nchem(i, j, k) <= 0)
                {
                    cout << "wrong chem " << '\t' << dtm << '\t' << dt << endl;
                    abort();
                }
            }
}








// 在当前时间步下，跨进程寻找全网格中最大的化学亚步迭代次数
// 便于分析化学反应计算是否成为了整个程序的性能瓶颈
void Reaction::PushNchemMax()
{
    NchemNow = Nchem.MaxValue();
    NchemMax_Rank.push_back(NchemNow);
}








// 通过 mpi 通信汇集所有计算节点的数据，计算全域在每个时间步最大的化学迭代步数，并将其保存到文件中用于负载平衡分析
void Reaction::GetNchemMax()
{
    int size = NchemMax_Rank.size();
    NchemMax_Total.Initial(size * numprocs);                              // 主进程的数组，用来存储所有进程中每个时间步自己区域内最大的亚步数
    MPI_Gather(&NchemMax_Rank[0], size, MPI_INT, &NchemMax_Total(0), size, MPI_INT, 0, MPI_COMM_WORLD);
    Array<int, 1> NchemMax_temp;
    NchemMax_temp.Initial(numprocs);
    if (myid == 0)                                                        // 寻找全局最值
    {
        for (int j = 0; j < size; j++)
        {
            for (int i = 0; i < numprocs; i++)
            {
                NchemMax_temp(i) = NchemTotal(j + i * size);
            }
            NchemMax.push_back(NchemMax_temp.MaxValue());
        }
        ofstream outfile("./output_" + to_string(numprocs) + "/load/Load_Loop.dat");
        for (int i = 1; i <= NchemMax.size(); i++)
            outfile << i << '\t' << NchemMax[i - 1] << '\n';
    }
}
















// 计算化学反应源项对物种浓度的偏导数
// 构建化学雅可比矩阵的对角线部分
// 为隐式数值求解器提供化学源项的对角雅可比矩阵，MD 矩阵会被用在解线性方程组 Ax=b 的过程中
void Reaction::Diagonalized(Array<double, 4> &Mc_temp, Array<double, 4> &Di_temp, Array<double, 3> &T, Array<double, 4> &Partial_T)
{

    for (int i = bc; i < ni + bc; i++)
        for (int j = bc; j < nj + bc; j++)
            for (int k = bc; k < nk + bc; k++)
            {
                interParaDiag __attribute__((aligned(64))) ip;
                const double T_local = T(i, j, k);
                for (int s = 0; s < NS; s++)
                {
                    Di(i, j, k, s) = Di_temp(i, j, k, s) * 1e-3;
                    Mc(i, j, k, s) = Mc_temp(i, j, k, s) * 1e-6;
                    ip.Hi[s] = GetHi(T_local, R, s, Coeff0, Coeff1);
                    ip.Si[s] = GetSi(T_local, R, s, Coeff0, Coeff1);
                    ip.Gi[s] = ip.Hi[s] - ip.Si[s] * T_local;
                }


                std::fill(ip.RR_F, ip.RR_F + NR, 1.0);
                std::fill(ip.RR_B, ip.RR_B + NR, 1.0);
                std::fill(ip.R_TB, ip.R_TB + NR, 1.0);
                for (int r = 0; r < NR; r++)
                {
                    double Xr = 0.0, delta_g = 0.0;
                    ip.KF[r] = Af(r) * pow(T_local, Bf(r)) * exp(-Eaf(r) / (Ru * T_local));
                    for (int s = 0; s < NS; s++)
                    {
                        Xr += (Stoi_B(s, r) - Stoi_F(s, r));
                        delta_g += (Stoi_B(s, r) - Stoi_F(s, r));
                        ip.RR_F[r] *= pow(Mc(i, j, k, s), Stoi_F(s, r));
                        ip.RR_B[r] *= pow(Mc(i, j, k, s), Stoi_B(s, r));
                        ip.R_TB[r] += Mc(i, j, k, s) * React_TB(s, r);
                    }
                    ip.Kp[r] = exp(-delta_g / (R * T_local)) * pow((P0 * 1e-6), Xr);
                    ip.Kc[r] = ip.Kp[r] * pow(R * T_local, -Xr);
                    ip.KB[r] = ip.KF[r] / ip.Kc[r];
                    if (ip.R_TB[r] = 0.0)
                        ip.R_TB[r] = 1.0;
                    ip.RR[r] = ip.KF[r] * ip.RR_F[r] - ip.KB[r] * ip.RR_B[r];
                }


                for (int s = 0; s < NS; s++)
                {
                    double temp = 0.0;
                    for (int r = 0; r < NR; r++)
                        temp += (Stoi_B(s, r) - Stoi_F(s, r)) * ip.R_TB[r] * ip.RR[r];
                    ip.Wi[s] = Mw(s) * temp;
                    CMS(i, j, k, s) = ip.Wi[s] * 1e3;
                }

                // 计算偏导数中间项
                for (int r = 0; r < NR; r++)
                    for (int s = 0; s < NS; s++)
                    {
                        if (Di(i, j, k ,s) == 0)
                        {
                            ip.WJH1[s][r] = 0.0;
                        }
                        else
                        {
                            ip.WJH1[s][r] = ip.KF[r] * Stoi_F(s, r) * ip.RR_F[r] / Di(i, j, k, s) - ip.KB[r] * Stoi_B(s, r) * ip.RR_B[r] / Di(i, j, k, s);
                        }
                        ip.WJH2[s][r] = ip.RR[r] * (Bf[r] / T_local + Eaf(r) / (Ru * pow(T_local, 2)));
                    }
                
                for (int s = 0; s < NS; s++)
                {
                    double temp = 0.0;
                    for (int r = 0; r < NR; r++)
                        temp += (Stoi_B(s, r) - Stoi_F(s, r)) * ((ip.WJH1[s][r] + ip.WJH2[s][r] * Partial_T(i, j, k, s)) * ip.R_TB[r] + ip.RR[r] * React_TB(s, r) / Mw(s));
                    MD(i, j, k, s, s) = Mw(s) * temp;                     // 构建对角矩阵 MD
                }
            }
}
















































































// 根据给定的温度 T，利用 NASA 7项多项式拟合公式，计算特定化学物种 SP 的定压比热容 Cp
double Reaction::GetCpi(double T, double R, int SP, Array<double, 2> &Coeff0, Array<double, 2> &Coeff1)
{
    double Cpi = 0.0;

    if (T < 1000)                                                         // 低温度段使用 Coeff0 数组中的系数
        Cpi = R * (Coeff0(0, SP) * pow(T, -2) + Coeff0(1, SP) * pow(T, -1) + Coeff0(2, SP) + Coeff0(3, SP) * T + Coeff0(4, SP) * pow(T, 2) + Coeff0(5, SP) * pow(T, 3) + Coeff0(6, SP) * pow(T, 4));
    else                                                                  // 高温度段使用 Coeff1 数组中的系数
        Cpi = R * (Coeff1(0, SP) * pow(T, -2) + Coeff1(1, SP) * pow(T, -1) + Coeff1(2, SP) + Coeff1(3, SP) * T + Coeff1(4, SP) * pow(T, 2) + Coeff1(5, SP) * pow(T, 3) + Coeff1(6, SP) * pow(T, 4));
    if (SP == 0)
        // cout << T << " " << R << " " << Coeff0(0, 0) << " " << Coeff0(1, 0) << " " << Coeff0(2, 0) << " " << Coeff0(3, 0) << endl;
    return Cpi;
}









// 根据给定的温度，计算特定化学物种 SP 的绝对焓 Hi
// NASA 7 项多项式的焓值积分公式
double Reaction::GetHi(double T, double R, int SP, Array<double, 2> &Coeff0, Array<double, 2> &Coeff1)
{
    double Hi = 0.0;

    if (T < 1000)
        Hi = R * (-Coeff0(0, SP) * pow(T, -1) + Coeff0(1, SP) * log(T) + Coeff0(2, SP) + Coeff0(3, SP) * pow(T, 2) / 2 + Coeff0(4, SP) * pow(T, 3) / 3 + Coeff0(5, SP) * pow(T, 4) / 4 + Coeff0(6, SP) * pow(T, 5) / 5 + Coeff0(7, SP));
    else
        Hi = R * (-Coeff1(0, SP) * pow(T, -1) + Coeff1(1, SP) * log(T) + Coeff1(2, SP) + Coeff1(3, SP) * pow(T, 2) / 2 + Coeff1(4, SP) * pow(T, 3) / 3 + Coeff1(5, SP) * pow(T, 4) / 4 + Coeff1(6, SP) * pow(T, 5) / 5 + Coeff0(7, SP));
    
    return Hi;
}









// 根据给定的温度，计算特定化学物种 SP 的标准状态熵 Si
// NASA 熵值公式
double Reaction::GetSi(double T, double R, int SP, Array<double, 2> &Coeff0, Array<double, 2> &Coeff1)
{
    double Si = 0.0;

    if (T < 1000)
        Si = R * (-Coeff0(0, SP) * pow(T, -2) / 2 - Coeff0(1, SP) * pow(T, -1) + Coeff0(2, SP) * log(T) + Coeff0(3, SP) * T + Coeff0(4, SP) * pow(T, 2) / 2 + Coeff0(5, SP) * pow(T, 3) / 3 + Coeff0(6, SP) * pow(T, 4) / 4 + Coeff0(8, SP));
    else
        Si = R * (-Coeff1(0, SP) * pow(T, -2) / 2 - Coeff1(1, SP) * pow(T, -1) + Coeff1(2, SP) * log(T) + Coeff1(3, SP) * T + Coeff1(4, SP) * pow(T, 2) / 2 + Coeff1(5, SP) * pow(T, 3) / 3 + Coeff1(6, SP) * pow(T, 4) / 4 + Coeff0(8, SP));
    
    return Si;
}











Array<double, 4> Reaction::GetWi()
{
    return CMS;
}











Array<double, 5> Reaction::GetMD()
{
    return MD;
}
// 在动态负载平衡或跨进程通信之前，将当前进程计算出的 Nchem 打包进一个连续的发送缓冲区 data 中
void Reaction::PackageUpdate(std::vector<int> &data, int meshnum, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    for (int n = 0; n < recvFrom.size(); n++)
    {   // 计算前 n-1 进程的累加和 -> 第 n 进程偏移量
        start = std::accumulate(transferNmesh.begin(), transferNmesh.begin() + recvFrom[n], 0);
        size = transferNmesh[recvFrom[n]];                                // 确认数据长度
        for (int i = 0; i < size; index++, i++)
            data[start++] = Nchem(index, 0, 0);
    }
}
// 把 data 输入缓冲区中一维连续的数据，根据坐标索引重新填回到本地的三维 Nchem(x, y, z) 数组中
void Reaction::UnpackageUpdate(std::vector<int> &data,std::vector<int> &transferIndex, int meshnum, std::vector<int> &sendTo, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    int x = 0, y = 0, z = 0;
    for (int n = 0; n < sendTo.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin(), transferNmesh.begin() + sendTo[n], 0);
        size = transferNmesh[sendTo[n]];
        for (int i = 0; i < size; index++, i++)
        {
            x = transferIndex[3 * index];                                 // transferIndex 是一维数组，存储着连续的 x y z 信息
            y = transferIndex[3 * index + 1];                             // 取出 index 点的 y
            z = transferIndex[3 * index + 2];
            Nchem(x, y, z) = data[start++];
        }
    }
}
// 将特定网格点 (i, j, k) 的 Nchem 压入 data 中 -> 三维映射
void Reaction::PackagePrev(std::vector<int> &data, int i, int j, int k)
{
    data.push_back(Nchem(i, j, k));
}
void Reaction::UnpackagePrev(std::vector<int> &data, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    for (int n = 0; n < recvFrom.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin() + recvFrom[0], transferNmesh.begin() + recvFrom[n], 0);
        size = transferNmesh[recvFrom[n]];
        for (int i = 0; i < size; index++, i++)
            Nchem(index, 0, 0) = data[start++];
    }
}
void Reaction::UnpackagePrev(std::vector<int> &data, int meshnum, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{
    // cout << "unpackage start\n";
    int index = 0, start = 0, size = 0;
    for (int n = 0; n < recvFrom.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin(), transferNmesh.begin() + recvFrom[n], 0);
        size = transferNmesh[recvFrom[n]];
        for (int i = 0; i < size; index++, i++)
            Nchem(index, 0, 0) = data[start++];
    }
}
