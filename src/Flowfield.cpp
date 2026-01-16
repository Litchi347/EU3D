







# include <iomanip>
# include <mpi.h>
# include <omp.h>
# include <numeric>





# include "Flowfield.hpp"
# include "FileReader.hpp"
# include "Differ.hpp"
# include "Global.hpp"











void Flowfield::InputRead(char *initialization, Array<double, 1> &xnode, Array<double, 1> &ynode, Array<double, 1> &znode, int bc)
{
    this->xnode = xnode;
    this->ynode = ynode;
    this->znode = znode;
    this->bc = bc;

    ni = xnode.GetSize() - 2 * bc;
    nj = ynode.GetSize() - 2 * bc;
    nk = znode.GetSize() - 2 * bc;

    FileReader Input;
    Input.readFile(initialization);

    NS = Input.getIntParameter("NS");
    NR = Input.getIntParameter("NR");


    Construction();



    P_pre = Input.getdoubleParameter("P_pre");
    T_pre = Input.getdoubleParameter("T_pre");
    M_wave = Input.getdoubleParameter("M_wave");

    return;
}











void Flowfield::Construction()
{
    int row = xnode.GetSize();
    int col = ynode.GetSize();
    int flo = znode.GetSize();

    U.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    V.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    W.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    P.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    D.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    T.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    C.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    Ma.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    Wav.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    Rgas.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    Cp.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    H.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    E.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    Gamma.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);

    Mr.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    
    Mc.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Mi.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Yi.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    Di.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);

    Cpi.Initial(NS);
    Hi.Initial(NS);
    Ei.Initial(NS);

    Uint.Initial(nj, nk);
    Vint.Initial(nj, nk);
    Wint.Initial(nj, nk);
    Pint.Initial(nj, nk);
    Dint.Initial(nj, nk);
    Tint.Initial(nj, nk);
    Hint.Initial(nj, nk);
    Eint.Initial(nj, nk);
    Gint.Initial(nj, nk);
    Yint.Initial(nj, nk, NS);

    CS.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 4);
    F.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 4);
    G.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 4);
    Q.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 4);

    Yi_temp0.Initial(NS);
    PLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    DLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    ULR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    VLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    WLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    HLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    GLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    YLR.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, 2);
    Yi_temp.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc);
    YL_temp.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    YR_temp.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS);
    YL.Initial(NS);
    YR.Initial(NS);

    Partion_T.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 4);
    RHS.Initial(ni + 2 * bc, nj + 2 * bc, nk + 2 * bc, NS + 4);

    send_data_1.Initial(9 * flo * col * bc + NS * flo * col * bc);
    recv_data_1.Initial(9 * flo * col * bc + NS * flo * col * bc);
    send_data_2.Initial(9 * flo * col * bc + NS * flo * col * bc);
    recv_data_2.Initial(9 * flo * col * bc + NS * flo * col * bc);
    send_data_3.Initial(9 * flo * col * bc + NS * flo * col * bc);
    recv_data_3.Initial(9 * flo * col * bc + NS * flo * col * bc);
    return;
}









// 初始化流体的物理状态：在 x=0.01 处放了一个瞬间移动的激波，并算出每一个点在 0 时刻应有的压力、温度、组分分布和总能量，为后续时间推进（Time Advancing）做好数据准备
// 设置了基础的流速、压力和温度
// 计算了多组分化学反应流所需的全部热力学参数和守恒变量
// 模拟一个激波管或激波诱导燃烧的初始环境
void Flowfield::FieldInitial(Array<double,1> &Ri_temp, Array<double,1> &Mw_temp, Array<double,4> &Mi_temp, Array<double,4> &Yi_temp, Array<double,2> &coeff0_temp, Array<double,2> &coeff1_temp)
{
    dt = 0.0;
    time = 0.0;

    this->Ri = Ri_temp;
    this->Mw = Mw_temp;
    this->Coeff0 = coeff0_temp;
    this->Coeff1 = coeff1_temp;

    W.Fill(0.0);

    // 利用空气动力学公式计算了波后的理论值（假设波前是理想气体）
    U_post = M_wave * (1 - (2 + 0.40 * pow(M_wave, 2)) / (2.40 * pow(M_wave, 2))) * sqrt(1.40 * 287 * T_pre);
    P_post = P_pre * (1 + 2 * 1.4 / 2.4 * (pow(M_wave, 2) - 1.0));
    T_post = T_pre * (2 + 0.40 * pow(M_wave, 2)) / (2.40 * pow(M_wave, 2)) * P_post /  P_pre;
    V.Fill(0.0);
    U.Fill(0.0);
    for(int i = 0; i < ni + 2 * bc; i++)
        for(int j = 0; j < nj + 2 * bc; j ++)
            for(int k = 0; k < nk + 2 * bc; k++)
                if(xnode(i) >= 0.01)             // 根据坐标 x = 0.01 作为分界点，将流场划分为两个区域：波前区域和波后区域
                {
                    P(i, j, k) = P_pre;          // 波前区域，压力、温度设置为起始状态
                    T(i, j, k) = T_pre;
                    U(i, j, k) = 0.0;
                }
                else
                {
                    P(i, j, k) = P_post;         // 波后区域，压力、温度等设为计算出来的激波后状态
                    T(i, j, k) = T_post;
                    U(i, j, k) = U_post;
                }
    // 对每个网格点进行复杂的物理量计算
    for(int i = 0; i < ni + bc; i++)
        for(int j = 0; j < nj + bc; j ++)
            for(int k = 0; k < nk + bc; k++)
            {
                Array<double, 1> Yi_temp1;
                Yi_temp1.Initial(NS);
                for(int s = 0;s < NS; s++)
                    Yi_temp1[s] = Yi_temp(i, j, k, s);
                Wav(i, j, k) = 1.0 / Fun.sum(3, Yi_temp1, Mw);
                Rgas(i, j, k) = R * 1000 / Wav(i, j, k);
                D(i, j, k) = P(i, j, k) / (Rgas(i, j, k) * T(i, j, k));
                for(int s = 0; s < NS;s++)
                {
                    Mi(i, j, k, s) = Mi_temp(i, j, k, s);
                    Yi(i, j, k, s) = Yi_temp(i, j, k, s);
                    Di(i, j, k, s) = Yi(i, j, k, s) * D(i, j, k);
                    Mc(i, j, k, s) = Di(i, j, k, s) / Mw(s) * 1000;
                    Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                    Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                    Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
                }
                Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp1) * 1000;
                H(i, j, k) = Fun.sum(2, Yi_temp1, Hi) * 1000;
                E(i, j, k) = Fun.sum(2, Yi_temp1, Ei) * 1000;
                Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - R * Fun.sum(3, Yi_temp1, Mw) * 1000);
                C(i, j, k) = sqrt(Gamma(i, j, k) * P(i, j, k) / D(i, j, k));
                // 守恒变量赋值
                for(int s = 0; s < NS; s++)
                    CS(i, j, k, s) = D(i, j, k) * Yi(i, j, k, s);
                CS(i, j, k, NS) = D(i, j, k) * U(i, j, k);
                CS(i, j, k, NS + 1) = D(i, j, k) * V(i, j, k);
                CS(i, j, k, NS + 2) = D(i, j, k) * W(i, j, k);
                CS(i, j, k, NS + 3) = D(i, j, k) * E(i, j, k) + 0.5 * (pow(U(i, j, k), 2) + pow(V(i, j, k),2) + pow(W(i, j, k),2)) * D(i, j, k);
            }


            GetPartial_T();
}










// 计算并同步整个计算域的自适应时间步长
void Flowfield::CFLcondition(double cfl, double Final_Time)
{
    dx = xnode(1) - xnode(0);
    dy = ynode(1) - ynode(0);
    dz = znode(1) - znode(0);

    double dt_conv = 1e10, temp = 0.0;
    // 计算局部最小时间步长
    for(int i = bc; i < ni + bc; i++)
        for(int j = bc; j < nj + bc; j++)
            for(int k = bc; k < nk + bc; k++)
            {
                temp = min(dx / (abs(U(i, j, k)) + C(i, j, k)), dy / (abs(V(i, j, k)) + C(i, j, k)));
                if(temp < dt_conv)
                    dt_conv = temp;
            }
    dt = cfl * dt_conv;                          // 将找到的时间步长乘以一个安全系数 cfl 确保即使流场发生剧烈变化，计算依然稳定

    if(time + dt > Final_Time)                   // 当当前步长加上去后会超过设定的总计算时间 Final_Time，则强行缩小步长，使程序精准停在结束时刻
        dt = Final_Time - time;
    
    // 全局并行同步，取最小值
    double dtc = 0;
    MPI_Reduce(&dt, &dtc, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    if(myid == 0)
        dt = dtc;
    MPI_Bcast(&dt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);                 // 强制所有进程在此同步，使用完全相同的时间步长进入下一步计算

    time += dt;                                  // 累加时间步长
}










void Flowfield::FieldBoundary_1dZND()
{




























































































}
void Flowfield::FieldBoundary_3dODM()            // 设置三维流场的物理边界条件，x 方向左边是进气口，右边是自由流出，y、z 方向底面和侧面被设定为反射壁面
{
    // 入口边界，模拟超音速入口或恒定的气流输入
    if(myid_x == 0)
        // 将 2D 的入口剖面数据平铺到 3D 数组最左侧的 Ghost Layers上
        for(int i = 0; i < bc; i++)              // Ghost Cells
            for(int j = bc; j < nj + bc; j++)    // 计算核心网格
                for(int k = bc; k < nk + bc; k++)
                {
                    U(i, j, k) = Uint(j - bc, k - bc);
                    V(i, j, k) = Vint(j - bc, k - bc);
                    W(i, j, k) = Wint(j - bc, k - bc);
                    P(i, j, k) = Pint(j - bc, k - bc);
                    D(i, j, k) = Dint(j - bc, k - bc);
                    T(i, j, k) = Tint(j - bc, k - bc);
                    H(i, j, k) = Hint(j - bc, k - bc);
                    Gamma(i, j, k) = Gint(j - bc, k - bc);
                    E(i, j, k) = Eint(j - bc, k - bc);
                    for(int s = 0; s < NS; s++)  // 已知量 预设好的
                        Yi(i, j, k, s) = Yint(j - bc, k - bc, s);
                }
    
    // 出口边界：使用镜像对称的方式在高速流模拟中，允许波扰动离开计算域而不会再边界产生明显的数值反射
    if(myid_x == (m_block_x - 1))
        for(int i = ni + bc; i < ni + 2 * bc; i++)
            for(int j = bc; j < nj + bc; j++)
                for(int k = bc; k < nk + bc; k++)
                {
                    U(i, j, k) = U(2 * (ni + bc) - 1 - i, j, k);
                    V(i, j, k) = V(2 * (ni + bc) - 1 - i, j, k);
                    W(i, j, k) = W(2 * (ni + bc) - 1 - i, j, k);
                    P(i, j, k) = P(2 * (ni + bc) - 1 - i, j, k);
                    T(i, j, k) = T(2 * (ni + bc) - 1 - i, j, k);

                    for(int s = 0; s < NS; s++)
                    {
                        Yi(i, j, k, s) = Yi(2 * (ni + bc) - 1 - i, j, k, s);
                        Yi_temp0(s) = Yi(i, j, k, s);
                        Mc(i, j, k, s) = Mc(2* (ni + bc) -1 -i, j, k, s);
                        Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
                    }
                
                    Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
                    D(i, j, k) = P(i, j, k) / (T(i, j, k) * Rgas(i, j, k));
                    Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
                    H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
                    E(i, j, k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
                    Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - Rgas(i, j, k));
                }
    
    // 固壁边界，激波反射会导致压力和温度剧烈跳跃，对于热力学变量需重新计算确保壁面反射符合能量守恒
    if(myid_y == 0)
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = 0; j < bc; j++)
                for(int k = bc; k < nk + bc; k++)
                {
                    U(i, j, k) = U(i, 2 * bc - 1 - j, k);
                    V(i, j, k) = -V(i, 2 * bc - 1 - j, k);
                    W(i, j, k) = W(i, 2 * bc - 1 - j, k);
                    P(i, j, k) = P(i, 2 * bc - 1 - j, k);
                    T(i, j, k) = T(i, 2 * bc - 1 - j, k);

                    for(int s = 0; s < NS; s++)
                    {
                        Yi(i, j, k, s) = Yi(i, 2 * bc - 1 -j, k, s);
                        Yi_temp0(s) = Yi(i, j, k, s);
                        Mc(i, j, k, s) = Mc(i, 2 * bc - 1 - j, k, s);
                        Cpi(s) = React.GetCpi(T(i, j, k),Ri(s), s, Coeff0, Coeff1);
                        Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
                    }

                    Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
                    D(i, j, k) = P(i ,j, k) / (T(i, j, k)* Rgas(i, j, k));
                    Cp(i, j ,k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
                    H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
                    E(i, j, k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
                    Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - Rgas(i, j, k));
                }
    // 另一个出口边界，对于热力学变量直接赋值：Y方向流场相对平稳，节省计算量
    if(myid_y == (m_block_y - 1))
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = nj + bc; j < nj + 2 * bc; j++)
                for(int k = bc; k < nk + bc ; k++)
                {
                    U(i, j, k) = U(i, 2 * (nj + bc) - 1 - j ,k);
                    V(i, j, k) = V(i, 2 * (nj + bc) - 1 - j ,k);
                    W(i, j, k) = W(i, 2 * (nj + bc) - 1 - j ,k);
                    P(i, j, k) = P(i, 2 * (nj + bc) - 1 - j ,k);
                    D(i, j, k) = D(i, 2 * (nj + bc) - 1 - j ,k);
                    T(i, j, k) = T(i, 2 * (nj + bc) - 1 - j ,k);
                    H(i, j, k) = H(i, 2 * (nj + bc) - 1 - j ,k);
                    Gamma(i, j, k) = Gamma(i, 2 * (nj + bc) - 1 - j ,k);
                    E(i, j, k) = E(i, 2 * (nj + bc) - 1 - j ,k);
                    for(int s = 0; s < NS; s++)
                        Yi(i, j, k, s) = Yi(i, 2 * (nj + bc) - 1 - j, k, s);
                }
    if(myid_z == 0)                              // 固壁边界，模拟一个封闭的、滑移的物理容器壁面
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = 0; j < nj + 2 * bc; j++)
                for(int k = 0; k < bc; k++)
                {













                    U(i, j, k) = U(i, j, 2 * bc - 1 - k);
                    V(i, j, k) = V(i, j, 2 * bc - 1 - k);
                    W(i, j, k) = -W(i, j, 2 * bc - 1 - k);
                    P(i, j, k) = P(i, j, 2 * bc - 1 - k);
                    T(i, j, k) = T(i, j, 2 * bc - 1 - k);

                    for(int s = 0; s < NS; s++)
                    {
                        Yi(i, j, k, s) = Yi(i, j, 2 * bc - 1 - k, s);
                        Yi_temp0(s) = Yi(i, j, k, s);
                        Mc(i, j, k, s) = Mc(i, j, 2 * bc - 1 - k, s);
                        Cpi(s) = React.GetCpi(T(i, j, k),Ri(s), s, Coeff0, Coeff1);
                        Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
                    }

                    Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
                    D(i, j, k) = P(i, j, k) / (T(i, j, k) * Rgas(i, j, k));
                    Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
                    H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
                    E(i, j, k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
                    Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - Rgas(i, j ,k));
                }


    if(myid_z == (m_block_z - 1))
        for(int i = 0;  i < ni + 2 * bc; i++)
            for(int j = 0; j < nj + 2 * bc; j++)
                for(int k = nk + bc; k < nk + 2 * bc; k++)
                {










                    U(i, j, k) = U(i, j, 2 * (nk + bc) - 1 - k);
                    V(i, j, k) = V(i, j, 2 * (nk + bc) - 1 - k);
                    W(i, j, k) = -W(i, j, 2 * (nk + bc) - 1 - k);
                    P(i, j, k) = P(i, j, 2 * (nk + bc) - 1 - k);
                    T(i, j, k) = T(i, j, 2 * (nk + bc) - 1 - k);

                    for(int s = 0; s < NS; s++)
                    {
                        Yi(i, j, k, s) = Yi(i, j, 2 * (nk + bc) - 1 - k, s);
                        Yi_temp0(s) = Yi(i, j, k, s);
                        Mc(i, j, k, s) = Mc(i, j, 2 * (nk + bc) - 1 - k, s);
                        Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
                    }

                    Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
                    D(i, j, k) = P(i, j, k) / (T(i, j, k) * Rgas(i, j, k));
                    Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
                    H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
                    E(i, j, k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
                    Gamma(i, j, k) = Cp(i, j ,k) / (Cp(i, j, k) - Rgas(i, j, k));
                }
}
void Flowfield::FieldBoundary_3dRSBI()           // 设置六个面的对称边界条件，并包含一个性能计时
{   // 全对称边界，激波发生全反射
    double t1, t2, t3, t4, t5, t6, t7;
    t1 = MPI_Wtime();

    if(myid_x == 0)
    {
        const int xl = 2 * bc - 1;
        // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
        for(int i = 0; i < bc; i++)
            for(int j = bc; j < nj + bc; j++)
                for(int k = bc; k < nk + bc; k++)
                {
                    U(i, j, k) = U(xl - i, j, k);
                    V(i, j, k) = V(xl - i, j, k);
                    W(i, j, k) = W(xl - i, j, k);
                    P(i, j, k) = P(xl - i, j, k);
                    D(i, j, k) = D(xl - i, j, k);
                    T(i, j, k) = T(xl - i, j, k);
                    H(i, j, k) = H(xl - i, j, k);
                    Gamma(i, j, k) = Gamma(xl - i, j, k);
                    E(i, j, k) = E(xl - i, j, k);
                    for(int s = 0; s < NS; s++)
                        Yi(i, j, k, s) = Yi(xl - i, j, k, s);
                }
    }

    t2 = MPI_Wtime();
    ft1 += t2 - t1;


    if(myid_x == (m_block_x - 1))
    {
        const int xr = 2 * (ni + bc) - 1;
        // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
        for(int i = ni + bc; i < ni + 2 * bc; i++)
            for(int j = bc; j < nj + bc; j++)
                for(int k = bc; k < nk + bc; k++)
                {
                    U(i, j, k) = U(xr - i, j, k);
                    V(i, j, k) = V(xr - i, j, k);
                    W(i, j, k) = W(xr - i, j, k);
                    P(i, j, k) = P(xr - i, j, k);
                    D(i, j, k) = D(xr - i, j, k);
                    T(i, j, k) = T(xr - i, j, k);
                    H(i, j, k) = H(xr - i, j, k);
                    Gamma(i, j, k) = Gamma(xr - i, j, k);
                    E(i, j, k) = E(xr - i, j, k);
                    for(int s = 0; s < NS; s++)
                        Yi(i, j, k, s) = Yi(xr - i, j, k, s);
                }
    }

    t3 = MPI_Wtime();
    ft2 += t3 - t2;


    if(myid_y == 0)
    {
        const int yf = 2 * bc - 1;
        // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
        for(int i = 0; i < ni + 2 * bc ;i++)
            for(int j = 0; j < bc; j++)
                for(int k = bc; k < nk + bc; k++)
                {
                    U(i, j, k) = U(i, yf - j, k);
                    V(i, j, k) = V(i, yf - j, k);
                    W(i, j, k) = W(i, yf - j, k);
                    P(i, j, k) = P(i, yf - j, k);
                    T(i, j, k) = T(i, yf - j, k);
                    
                    for(int s = 0; s < NS; s++)
                    {
                        Yi(i, j, k, s) = Yi(i, yf - j, k, s);
                        Yi_temp0(s) = Yi(i, j, k, s);
                        Mc(i, j, k, s) = Mc(i, yf - j, k ,s);
                        Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                        Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
                    }

                    Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
                    D(i, j, k) = P(i, j, k) / (T(i, j, k) * Rgas(i, j, k));
                    Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
                    H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
                    E(i, j ,k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
                    Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - Rgas(i, j, k));
                }
    }

    t4 = MPI_Wtime();
    ft3 += t4 - t3;


    if(myid_y == (m_block_y - 1))
    {
        const int yb = 2 * (nj + bc) - 1;
        // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = nj + bc; j < nj + 2 * bc; j++)
                for(int k =  bc; k < nk + bc; k++)
                {
                    U(i, j, k) = U(i, yb - j, k);
                    V(i, j, k) = V(i, yb - j, k);
                    W(i, j, k) = W(i, yb - j, k);
                    P(i, j, k) = P(i, yb - j, k);
                    D(i, j, k) = D(i, yb - j, k);
                    T(i, j, k) = T(i, yb - j, k);
                    H(i, j, k) = H(i, yb - j, k);
                    Gamma(i, j, k) = Gamma(i, yb - j, k);
                    E(i, j, k) = E(i, yb - j, k);
                    for(int s = 0; s < NS; s++)
                        Yi(i, j, k, s) = Yi(i, yb - j, k, s);
                }
    }

    t5 = MPI_Wtime();
    ft4 += t5 - t4;


    if(myid_z == 0)
    {
        const int zd = 2 * bc - 1;
        // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = 0; j < nj + 2 * bc; j++)
                for(int k = 0; k < bc; k++)
                {
                    U(i, j, k) = U(i, j, zd - k);
                    V(i, j, k) = V(i, j, zd - k);
                    W(i, j, k) = W(i, j, zd - k);
                    P(i, j, k) = P(i, j, zd - k);
                    D(i, j, k) = D(i, j, zd - k);
                    T(i, j, k) = T(i, j, zd - k);
                    H(i, j, k) = H(i, j, zd - k);
                    Gamma(i, j, k) = Gamma(i, j, zd - k);
                    E(i, j, k) = E(i, j, zd - k);
                    for(int s = 0; s < NS; s++)
                        Yi(i, j, k, s) = Yi(i, j, zd - k, s);
					// U(i, j, k) = U(i, j, zd - k);
					// V(i, j, k) = V(i, j, zd - k);
					// W(i, j, k) = -W(i, j, zd - k);
					// P(i, j, k) = P(i, j, zd - k);
					// T(i, j, k) = T(i, j, zd - k);

					// for (int s = 0; s < NS; s++)
					// {
					// 	Yi(i, j, k, s) = Yi(i, j, zd - k, s);
					// 	Yi_temp0(s) = Yi(i, j, k, s);
					// 	Mc(i, j, k, s) = Mc(i, j, zd - k, s);
					// 	Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
					// 	Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
					// 	Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
					// }

					// Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
					// D(i, j, k) = P(i, j, k) / (T(i, j, k) * Rgas(i, j, k));
					// Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
					// H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
					// E(i, j, k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
					// Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - Rgas(i, j, k));
                }
    }

    t6 = MPI_Wtime();
    ft5 += t6 - t5;


    if(myid_z == (m_block_z - 1))
    {
        const int zu = 2 * (nk + bc) - 1;
        // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = 0; j < nj + 2 * bc; j++)
                for(int k = nk + bc; k < nk + 2 * bc; k++)
                {
                    U(i, j, k) = U(i, j, zu - k);
                    V(i, j, k) = V(i, j, zu - k);
                    W(i, j, k) = W(i, j, zu - k);
                    P(i, j, k) = P(i, j, zu - k);
                    D(i, j, k) = D(i, j, zu - k);
                    T(i, j, k) = T(i, j, zu - k);
                    H(i, j, k) = H(i, j, zu - k);
                    Gamma(i, j, k) = Gamma(i, j, zu - k);
                    E(i, j, k) = E(i, j, zu - k);
                    for(int s = 0; s < NS; s++)
                        Yi(i, j, k, s) = Yi(i, j, zu - k, s);
					// U(i, j, k) = U(i, j, zu - k);
					// V(i, j, k) = V(i, j, zu - k);
					// W(i, j, k) = -W(i, j, zu - k);
					// P(i, j, k) = P(i, j, zu - k);
					// T(i, j, k) = T(i, j, zu - k);

					// for (int s = 0; s < NS; s++)
					// {
					// 	Yi(i, j, k, s) = Yi(i, j, zu - k, s);
					// 	Yi_temp0(s) = Yi(i, j, k, s);
					// 	Mc(i, j, k, s) = Mc(i, j, zu - k, s);
					// 	Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
					// 	Hi(s) = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
					// 	Ei(s) = Hi(s) - Ri(s) * T(i, j, k);
					// }

					// Rgas(i, j, k) = Fun.sum(2, Yi_temp0, Ri) * 1000;
					// D(i, j, k) = P(i, j, k) / (T(i, j, k) * Rgas(i, j, k));
					// Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp0) * 1000;
					// H(i, j, k) = Fun.sum(2, Yi_temp0, Hi) * 1000;
					// E(i, j, k) = Fun.sum(2, Yi_temp0, Ei) * 1000;
					// Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - Rgas(i, j, k));
                }
    }

    t7 = MPI_Wtime();
    ft6 += t7 - t6;
}

void Flowfield::Mpi_Boundary()                   // 实现相邻进程之间的数据块交换
{
    const int row = xnode.GetSize();
    const int col = ynode.GetSize();
    const int flo = znode.GetSize();
    const int xsize = flo * col * bc;
    const int ysize = flo * row * bc;
    const int zsize = col * row * bc;


    int tag1 = 0, tag2 = 100, tag3 = 200, tag4 = 300, tag5 = 400, tag6 = 500;
    int idx = 0;

    if(m_block_x > 1)
    {
        // 打包紧邻右边界的数据块，发给右邻居去填充它们的左 Ghost Cells
        idx = 0;
        for(int i = 0; i < bc; i++)
            for(int j = 0; j < col; j++)
                for(int k = 0; k < flo; k++)
                {
                    send_data_1(idx) = D(ni - i, j ,k);
                    send_data_1(idx + xsize) = U(ni - i, j ,k);
                    send_data_1(idx + 2 * xsize) = V(ni - i, j, k);
                    send_data_1(idx + 3 * xsize) = W(ni - i, j, k);
                    send_data_1(idx + 4 * xsize) = P(ni - i, j, k);
                    send_data_1(idx + 5 * xsize) = T(ni - i, j, k);
                    send_data_1(idx + 6 * xsize) = H(ni - i, j, k);
                    send_data_1(idx + 7 * xsize) = Gamma(ni - i, j, k);
                    send_data_1(idx + 8 * xsize) = C(ni - i, j, k);
                    for(int s = 0; s < NS; s++)
                        send_data_1(9 * xsize + idx * NS + s) = Yi(ni - i, j, k, s);
                    idx++;
                }

        MPI_Request request1[2];
        MPI_Status status1[2];
        MPI_Isend(&send_data_1(0), send_data_1.GetSize(), MPI_DOUBLE, m_right, tag1, MPI_COMM_WORLD, &request1[0]);
        MPI_Irecv(&recv_data_1(0), recv_data_1.GetSize(), MPI_DOUBLE, m_left, tag1, MPI_COMM_WORLD, &request1[1]);
        MPI_Waitall(2, request1, status1);
        // 从左邻居的右边界数据块获得数据解包到左 Ghost Cells
        idx = 0;
        if(myid_x != 0)
            for(int i = 0; i < bc; i++)
                for(int j = 0; j < col; j++)
                    for(int k = 0; k < flo; k++)
                    {
                        D(bc - i - 1, j, k) = recv_data_1(idx);
                        U(bc - i - 1, j, k) = recv_data_1(idx + xsize);
                        V(bc - i - 1, j, k) = recv_data_1(idx + 2 * xsize);
                        W(bc - i - 1, j, k) = recv_data_1(idx + 3 * xsize);
                        P(bc - i - 1, j, k) = recv_data_1(idx + 4 * xsize);
                        T(bc - i - 1, j, k) = recv_data_1(idx + 5 * xsize);
                        H(bc - i - 1, j, k) = recv_data_1(idx + 6 * xsize);
                        Gamma(bc - i - 1, j, k) = recv_data_1(idx + 7 * xsize);
                        C(bc - i - 1, j, k) = recv_data_1(idx + 8 * xsize);
                        for(int s = 0; s < NS; s++)
                            Yi(bc - i - 1, j, k, s) = recv_data_1(9 * xsize + idx * NS + s);
                        idx++;
                    }


        // 打包紧邻左边界的数据块，发给左邻居去填充它们的右 Ghost Cells
        idx = 0;
        for(int i = 0; i < bc; i++)
            for(int j = 0; j < col; j++)
                for(int k = 0; k < flo; k++)
                {
                    send_data_1(idx) = D(bc + 1 + i, j, k);
                    send_data_1(idx + xsize) = U(bc + 1 + i, j, k);
                    send_data_1(idx + 2 * xsize) = V(bc + 1 + i, j, k);
                    send_data_1(idx + 3 * xsize) = W(bc + 1 + i, j, k);
                    send_data_1(idx + 4 * xsize) = P(bc + 1 + i, j, k);
                    send_data_1(idx + 5 * xsize) = T(bc + 1 + i, j, k);
                    send_data_1(idx + 6 * xsize) = H(bc + 1 + i, j, k);
                    send_data_1(idx + 7 * xsize) = Gamma(bc + 1 + i, j, k);
                    send_data_1(idx + 8 * xsize) = C(bc + 1 + i, j, k);
                    for(int s = 0; s < NS; s++)
                        send_data_1(9 * xsize + idx * NS + s) = Yi(bc + 1 + i,j, k,s);
                    idx++;
                }

        MPI_Request request2[2];
        MPI_Status status2[2];
        MPI_Isend(&send_data_1(0), send_data_1.GetSize(), MPI_DOUBLE, m_right, tag2, MPI_COMM_WORLD, &request2[0]);
        MPI_Irecv(&recv_data_1(0), recv_data_1.GetSize(), MPI_DOUBLE, m_left, tag2, MPI_COMM_WORLD, &request2[1]);
        MPI_Waitall(2, request2, status2);
        // 从右邻居的左边界数据块获得数据解包到右 Ghost Cells
        idx = 0;
        if(myid_x != (m_block_x - 1))
            for(int i = 0; i < bc; i++)
                for(int j = 0; j < col; j++)
                    for(int k = 0; k < flo; k++)
                    {
                        D(ni + bc + i, j, k) = recv_data_1(idx);
                        U(ni + bc + i, j, k) = recv_data_1(idx + xsize);
                        V(ni + bc + i, j, k) = recv_data_1(idx + 2 * xsize);
                        W(ni + bc + i, j, k) = recv_data_1(idx + 3 * xsize);
                        P(ni + bc + i, j, k) = recv_data_1(idx + 4 * xsize);
                        T(ni + bc + i, j, k) = recv_data_1(idx + 5 * xsize);
                        H(ni + bc + i, j, k) = recv_data_1(idx + 6 * xsize);
                        Gamma(ni + bc + i, j, k) = recv_data_1(idx + 7 * xsize);
                        C(ni + bc + i, j, k) = recv_data_1(idx + 8 * xsize);
                        for(int s = 0; s < NS; s++)
                            Yi(ni + bc + i, j, k, s) = recv_data_1(9 * xsize + idx * NS + s);
                        idx++;
                    }
    }

    if(m_block_y > 1)
    {

        
        idx = 0;
        for(int i = 0; i < bc; i++)
                for(int j = 0; j < col; j++)
                    for(int k = 0; k < flo; k++)
                    {
                        send_data_2(idx) = D(i, nj - j, k);
                        send_data_2(idx + ysize) = U(i, nj - j, k);
                        send_data_2(idx + 2 * ysize) = V(i, nj - j, k);
                        send_data_2(idx + 3 * ysize) = W(i, nj - j, k);
                        send_data_2(idx + 4 * ysize) = P(i, nj - j, k);
                        send_data_2(idx + 5 * ysize) = T(i, nj - j, k);
                        send_data_2(idx + 6 * ysize) = H(i, nj - j, k);
                        send_data_2(idx + 7 * ysize) = Gamma(i, nj - j, k);
                        send_data_2(idx + 8 * ysize) = C(i, nj - j, k);
                        for(int s = 0; s < NS; s++)
                            send_data_2(9 * ysize + idx * NS + s) = Yi(i, nj - j, k,s);
                        idx++;
                    }

        MPI_Request request3[2];
        MPI_Status status3[2];
        MPI_Isend(&send_data_2(0), send_data_2.GetSize(), MPI_DOUBLE, m_back, tag3, MPI_COMM_WORLD, &request3[0]);
        MPI_Irecv(&recv_data_2(0), recv_data_2.GetSize(), MPI_DOUBLE, m_front, tag3, MPI_COMM_WORLD, &request3[1]);
        MPI_Waitall(2, request3, status3);

        idx = 0;
        if(myid_y != 0)
            for(int i = 0; i < row; i++)
                for(int j = 0; j < bc; j++)
                    for(int k = 0; k < flo; k++)
                    {
                        D(i, bc - j - 1, k) = recv_data_2(idx);
                        U(i, bc - j - 1, k) = recv_data_2(idx + ysize);
                        V(i, bc - j - 1, k) = recv_data_2(idx + 2 * ysize);
                        W(i, bc - j - 1, k) = recv_data_2(idx + 3 * ysize);
                        P(i, bc - j - 1, k) = recv_data_2(idx + 4 * ysize);
                        T(i, bc - j - 1, k) = recv_data_2(idx + 5 * ysize);
                        H(i, bc - j - 1, k) = recv_data_2(idx + 6 * ysize);
                        Gamma(i, bc - j - 1, k) = recv_data_2(idx + 7 * ysize);
                        C(i, bc - j - 1, k) = recv_data_2(idx + 8 * ysize);
                        for(int s = 0; s < NS; s++)
                            Yi(i, bc - j - 1, k, s) = recv_data_2(9 * ysize + idx * NS + s);
                        idx++;
                    }



        idx = 0;
        for(int i = 0; i < row; i++)
                for(int j = 0; j < bc; j++)
                    for(int k = 0; k < flo; k++)
                    {
                        send_data_2(idx) = D(i, bc + j + 1, k);
                        send_data_2(idx + ysize) = U(i, bc + j + 1, k);
                        send_data_2(idx + 2 * ysize) = V(i, bc + j + 1, k);
                        send_data_2(idx + 3 * ysize) = W(i, bc + j + 1, k);
                        send_data_2(idx + 4 * ysize) = P(i, bc + j + 1, k);
                        send_data_2(idx + 5 * ysize) = T(i, bc + j + 1, k);
                        send_data_2(idx + 6 * ysize) = H(i, bc + j + 1, k);
                        send_data_2(idx + 7 * ysize) = Gamma(i, bc + j + 1, k);
                        send_data_2(idx + 8 * ysize) = C(i, bc + j + 1, k);
                        for(int s = 0; s < NS; s++)
                            send_data_2(9 * ysize + idx * NS + s) = Yi(i, bc + j + 1, k, s);
                        idx++;
                    }

        MPI_Request request4[2];
        MPI_Status status4[2];
        MPI_Isend(&send_data_2(0), send_data_2.GetSize(), MPI_DOUBLE, m_front, tag4, MPI_COMM_WORLD, &request4[0]);
        MPI_Irecv(&recv_data_2(0), recv_data_2.GetSize(), MPI_DOUBLE, m_back, tag4, MPI_COMM_WORLD, &request4[1]);
        MPI_Waitall(2, request4, status4);

        idx = 0;
        if(myid_y != (m_block_y - 1))
            for(int i = 0; i < row; i++)
                for(int j = 0; j < bc; j++)
                    for(int k = 0; k < flo; k++)
                    {
                        D(i, nj + bc + j, k) = recv_data_2(idx);
                        U(i, nj + bc + j, k) = recv_data_2(idx + ysize);
                        V(i, nj + bc + j, k) = recv_data_2(idx + 2 * ysize);
                        W(i, nj + bc + j, k) = recv_data_2(idx + 3 * ysize);
                        P(i, nj + bc + j, k) = recv_data_2(idx + 4 * ysize);
                        T(i, nj + bc + j, k) = recv_data_2(idx + 5 * ysize);
                        H(i, nj + bc + j, k) = recv_data_2(idx + 6 * ysize);
                        Gamma(i, nj + bc + j, k) = recv_data_2(idx + 7 * ysize);
                        C(i, nj + bc + j, k) = recv_data_2(idx + 8 * ysize);
                        for(int s = 0; s < NS; s++)
                            Yi(i, nj + bc + j, k, s) = recv_data_2(9 * ysize + idx * NS + s);
                        idx++;
                    }
    }

    if(m_block_z > 1)
    {


        idx = 0;
        for(int i = 0; i < row; i++)
            for(int j = 0; j < col; j++)
                for(int k = 0; k < bc; k++)
                {
                    send_data_3(idx) = D(i, j, nk - k);
                    send_data_3(idx + zsize) = U(i, j, nk - k);
                    send_data_3(idx + 2 * zsize) = V(i, j, nk - k);
                    send_data_3(idx + 3 * zsize) = W(i, j, nk - k);
                    send_data_3(idx + 4 * zsize) = P(i, j, nk - k);
                    send_data_3(idx + 5 * zsize) = T(i, j, nk - k);
                    send_data_3(idx + 6 * zsize) = H(i, j, nk - k);
                    send_data_3(idx + 7 * zsize) = Gamma(i, j, nk - k);
                    send_data_3(idx + 8 * zsize) = C(i, j, nk - k);
                    for(int s = 0; s < NS; s++)
                        send_data_3(9 * zsize + idx * NS + s) = Yi(i, j, nk - k, s);
                    idx++;
                }
        
        MPI_Request request5[2];
        MPI_Status status5[2];
        MPI_Isend(&send_data_3(0), send_data_3.GetSize(), MPI_DOUBLE, m_up, tag5, MPI_COMM_WORLD, &request5[0]);
        MPI_Irecv(&send_data_3(0), send_data_3.GetSize(), MPI_DOUBLE, m_down, tag5, MPI_COMM_WORLD, &request5[1]);
        MPI_Waitall(2, request5, status5);

        idx = 0;
        if(myid_z != 0)
            for(int i = 0; i < row; i++)
                for(int j = 0; j < col; j++)
                    for(int k = 0; k < bc; k++)
                    {
                        D(i, j, bc - k - 1) = recv_data_3(idx);
                        U(i, j, bc - k - 1) = recv_data_3(idx + zsize);
                        V(i, j, bc - k - 1) = recv_data_3(idx + 2 * zsize);
                        W(i, j, bc - k - 1) = recv_data_3(idx + 3 * zsize);
                        P(i, j, bc - k - 1) = recv_data_3(idx + 4 * zsize);
                        T(i, j, bc - k - 1) = recv_data_3(idx + 5 * zsize);
                        H(i, j, bc - k - 1) = recv_data_3(idx + 6 * zsize);
                        Gamma(i, j, bc - k - 1) = recv_data_3(idx + 7 * zsize);
                        C(i, j, bc - k - 1) = recv_data_3(idx + 8 * zsize);
                        for(int s = 0; s < NS; s++)
                            Yi(i, j, bc - k - 1, s) = recv_data_3(9 * zsize + idx * NS + s);
                        idx++;
                    }



        idx = 0;
        for(int i = 0; i < row; i++)
            for(int j = 0; j < col; j++)
                for(int k = 0; k < bc; k++)
                {
                    send_data_3(idx) = D(i, j, bc + k + 1);
                    send_data_3(idx + zsize) = U(i, j, bc + k + 1);
                    send_data_3(idx + 2 * zsize) = V(i, j, bc + k + 1);
                    send_data_3(idx + 3 * zsize) = W(i, j, bc + k + 1);
                    send_data_3(idx + 4 * zsize) = P(i, j, bc + k + 1);
                    send_data_3(idx + 5 * zsize) = T(i, j, bc + k + 1);
                    send_data_3(idx + 6 * zsize) = H(i, j, bc + k + 1);
                    send_data_3(idx + 7 * zsize) = Gamma(i, j, bc + k + 1);
                    send_data_3(idx + 8 * zsize) = C(i, j, bc + k + 1);
                    for(int s = 0; s < NS; s++)
                        send_data_3(9 * zsize + idx * NS + s) = Yi(i, j, bc + k + 1, s);
                    idx++;
                }
         
        MPI_Request request6[2];
        MPI_Status statu6[2];
        MPI_Isend(&send_data_3(0), send_data_3.GetSize(), MPI_DOUBLE, m_down, tag6, MPI_COMM_WORLD, &request6[0]);
        MPI_Irecv(&send_data_3(0), send_data_3.GetSize(), MPI_DOUBLE, m_up, tag6, MPI_COMM_WORLD, &request6[1]);
        MPI_Waitall(2, request6, statu6);

        idx = 0;
        if(myid_z != (m_block_z - 1))
            for(int i = 0; i < row; i++)
                for(int j = 0; j < col; j++)
                    for(int k = 0; k < bc; k++)
                    {
                        D(i, j, nk + bc + k) = recv_data_3(idx);
                        U(i, j, nk + bc + k) = recv_data_3(idx + zsize);
                        V(i, j, nk + bc + k) = recv_data_3(idx + 2 * zsize);
                        W(i, j, nk + bc + k) = recv_data_3(idx + 3 * zsize);
                        P(i, j, nk + bc + k) = recv_data_3(idx + 4 * zsize);
                        T(i, j, nk + bc + k) = recv_data_3(idx + 5 * zsize);
                        H(i, j, nk + bc + k) = recv_data_3(idx + 6 * zsize);
                        Gamma(i, j, nk + bc + k) = recv_data_3(idx + 7 * zsize);
                        C(i, j, nk + bc + k) = recv_data_3(idx + 8 * zsize);
                        for(int s = 0; s < NS; s++)
                            Yi(i, j, nk + bc + k, s) = recv_data_3(9 * zsize + idx * NS + s);
                        idx++;
                    }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (T.IsNan() || D.IsNan())
    {
        cout << "sth wrong\n";
        abort();
    }
}











// 空间离散 Diff 与时间推进 TimeAdv 的选择
void Flowfield::Advection(int _TimeAdv, int _Diff)
{
    void (*diff)(int, Array<double, 1> &, Array<double, 1> &, Array<double, 1> &, Array<double, 4> &, Array<double, 3> &, int, int);

    // 利用函数指针 duff 选择 MUSCL 格式
    switch(_Diff)
    {
    case 0:                                      // 一阶精度
        diff = MUSCL_1;
        break;
    case 1:                                      // 二阶精度
        diff = MUSCL_2;
        break;
    }
    F.Fill(0.0);
    G.Fill(0.0);
    Q.Fill(0.0);

    AUSM(1, F, diff);                            // 计算 X 方向通量 F
    AUSM(2, G, diff);                            // 计算 Y 方向通量 G
    AUSM(3, Q, diff);                            // 计算 Z 方向通量 Q

    // 时间积分
    switch (_TimeAdv)
    {
    case 0:                                      // 显示欧拉法
        Time.EE(NS, dt, xnode, ynode, znode, bc, F, G, Q, RHS);
        break;
    case 1:                                      // 三阶 TVD 龙塔-库塔法
        Time.TVD_RK3(NS, dt, xnode, ynode, znode, bc, F, G, Q, RHS);
        break;
    }
}










// 执行 CS 的最终更新 显式时间推进步
void Flowfield::Update_after_Adv()
{
    // #pragma omp simd
    for(int i = bc; i < ni + bc; i++)            // 只更新数据块， Ghost Cells 的值必须由 Mpi_Boundary 从邻居进程传过来
        for(int j = bc; j < nj + bc; j++)
            for(int k = bc; k < nk + bc; k++)
                for(int s = 0; s < NS; s++)      // 保守变量向量的长度
                    CS(i, j, k, s) = CS(i, j, k ,s) + RHS(i, j, k, s);
	// if (RHS.IsNan())
	// {
	// 	cout << "RHS is nan\n";
	// }
	// if (CS.IsNan())
	// {
	// 	cout << myid << " CS is nan\n";
	// }
    Update_after_CS();
}







// 在 CFD 求解器中，通常我们会存储两套变量
// 原始变量：如压力 P、温度 T、速度 U, V, W、组分质量分数 Y_i。这些变量直观、便于设置边界条件
// 保守变量：如密度 rho、动量 rho u、总能 rho E。这是控制方程直接求解的对象，满足守恒律
// 将原始变量转换并初始化为 CS
void Flowfield::Explicit()
{
    for(int i = bc; i < ni + bc; i++)
        for(int j = bc; j < nj + bc; j++)
            for(int k = bc; k < nk + bc; k++)
            {
                D(i, j, k) = 0.0;
                for(int s = 0; s < NS; s++)
                {
                    CS(i, j, k, s) = Di(i, j, k, s);
                    D(i, j, k) += CS(i, j, k, s);// 计算总密度
                }
                CS(i, j, k, NS) = U(i, j, k) * D(i, j, k);
                CS(i, j, k, NS + 1) = V(i, j, k) * D(i, j, k);
                CS(i, j, k, NS + 2) = D(i, j, k) * E(i, j, k) + 0.5 * (pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2)) * D(i, j, k);
            }
    Update_after_CS();
}

// 特定的点 
void Flowfield::Explicit(int step, int i, int j, int k)
{
    D(i, j, k) = 0.0;
    for(int s = 0; s < NS; s++)
    {
        CS(i, j, k, s) = Di(i, j, k, s);
        D(i, j, k) += CS(i, j, k, s);
    }
    CS(i, j, k, NS) = U(i, j, k) * D(i, j, k);
    CS(i, j, k, NS + 1) = V(i, j, k) * D(i, j, k);
    CS(i, j, k, NS + 2) = W(i, j, k) * D(i, j, k);
    CS(i, j, k, NS + 3) = D(i, j, k) * E(i, j, k) + 0.5 * (pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2)) * D(i, j, k);

    Update_after_CS(step, i, j, k);
}










// 利用 隐式处理反应源项-显式处理流动项 方法处理化学反应源项并更新流场物理量
void Flowfield::Update_IMEX(Array<double, 4> &Wi, Array<double, 5> &MD)
{


    for(int i = bc; i < ni + bc; i++)
        for(int j = bc; j < nj + bc; j++)
            for(int k = bc; k < nk + bc; k++)
            {   // 累加源项
                for(int s = 0; s < NS; s++)
                    RHS(i, j, k, s) = RHS(i, j, k, s) + Wi(i, j, k, s) * dt;
                // 隐式修正与变量更新
                for(int s = 0; s < NS + 4; s++)
                    CS(i, j, k, s) = CS(i, j, k, s) + RHS(i, j, k, s) / (1.0 - dt * MD(i, j, k, s, s));
            }            

    Update_after_CS();                           // 将更新后的保守变量转换为原始变量(P, U, V, W, E)


    GetPartial_T();                              // 专门计算或更新各组分相关的温度分量/偏导
}









// 解码：将更新后的保守变量转换为原始变量(P, U, V, W, E)
// 控制方程更新的是保守变量，需要原始变量来导出下一步的通量
void Flowfield::Update_after_CS()
{
    // 全局更新
    // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
    for(int i = bc; i < ni + bc; i++)
        for(int j = bc; j < nj + bc; j ++)
            for(int k = bc; k < nk + bc; k++)
            {
                Array<double, 1> Yi_temp1(NS);
                Array<double, 1> Cpi1(NS);
                D(i, j, k) = 0;
                for(int s = 0; s < NS; s++)
                {
                    D(i, j, k) += CS(i, j, k, s);
                    Di(i, j, k, s) = CS(i, j, k, s);
                    Mc(i, j, k, s) = Di(i, j, k, s) / Mw(s) * 1000;
                }
                for(int s = 0; s < NS; s++)
                {
                    Yi(i, j, k, s) = CS(i, j, k, s) / D(i, j, k);
                    Yi_temp1(s) = Yi(i, j, k, s);
                }
                T(i, j ,k) = Get_temp(T(i, j, k), i, j, k, Yi_temp1);
                U(i, j ,k) = CS(i, j, k, NS) / D(i, j, k);
                V(i, j, k) = CS(i, j, k, NS + 1) / D(i, j, k);
                W(i, j, k) = CS(i, j, k, NS + 2) / D(i, j, k);
                E(i, j, k) = (CS(i, j, k, NS + 3) - 0.5 * D(i, j, k) * (pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2))) / D(i, j, k);
                Wav(i, j ,k) = 1.0 / Fun.sum(3, Yi_temp1, Mw);
                Rgas(i, j, k) = R * 1000 / Wav(i, j, k);
                P(i, j, k) = T(i, j, k) * D(i, j, k) * Rgas(i, j, k);
                H(i, j, k) = P(i, j, k) / D(i, j, k) + E(i, j, k);
                for(int s = 0; s < NS;s++)
                    Cpi1(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                Cp(i, j, k) = Fun.sum(2, Cpi1, Yi_temp1) * 1000;
                Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - R * Fun.sum(3, Yi_temp1, Mw) * 1000);
                C(i, j, k) = sqrt(Gamma(i, j, k) * P(i, j, k) / D(i, j, k));
                Ma(i, j, k) = sqrt(pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2)) / C(i, j, k);
            }
}


void Flowfield::Update_after_CS(int step, int i, int j, int k)
{
    // 仅针对单个网格点
    Array<double, 1> Yi_temp1(NS);

    D(i, j, k) = 0;
    for(int s = 0; s < NS; s++)
    {
        D(i, j ,k) += CS(i, j, k,s);
        Di(i, j, k, s) = CS(i, j, k, s);
        Mc(i, j, k, s) = Di(i, j, k, s) / Mw(s) * 1000;
    }
    for(int s = 0; s < NS; s++)
    {
        Yi(i, j, k, s) = CS(i, j, k, s) / D(i, j, k);
        Yi_temp1(s) = Yi(i, j, k, s);
    }
    T(i, j, k) = Get_temp(T(i, j, k), i, j, k, Yi_temp1);
    U(i, j, k) = CS(i, j, k, NS) / D(i, j, k);
    V(i, j, k) = CS(i, j, k, NS + 1) / D(i, j, k);
    W(i, j, k) = CS(i, j, k, NS + 2) / D(i, j, k);
    E(i, j, k) = (CS(i, j, k, NS + 3) - 0.5 * D(i, j, k) * (pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2))) / D(i, j, k);
    Wav(i, j, k) = 1.0 / Fun.sum(3, Yi_temp1, Mw);
    Rgas(i, j, k) = R * 1000 / Wav(i, j, k);
    P(i, j, k) = T(i, j, k) * D(i, j, k) * Rgas(i, j, k);
    H(i, j, k) = P(i, j, k) / D(i, j, k) + E(i, j, k);
    for(int s = 0; s < NS; s++)
        Cpi(s) = React.GetCpi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
    Cp(i, j, k) = Fun.sum(2, Cpi, Yi_temp1) * 1000;
    Gamma(i, j, k) = Cp(i, j, k) / (Cp(i, j, k) - R * Fun.sum(3, Yi_temp1, Mw) * 1000);
    C(i, j, k) = sqrt(Gamma(i, j, k) * P(i, j, k) / D(i, j, k));
    Ma(i, j, k) = sqrt(pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2)) / C(i, j, k);
}









// 计算对流通量，实现了 AUSM 格式：将通量分解成对流部分（由速度主导）和压力部分（由声波主导）
void Flowfield::AUSM(int direction, Array<double, 4> &Fi, void (*Diff)(int, Array<double, 1> &, Array<double, 1> &, Array<double, 1> &, Array<double, 4> &, Array<double, 3> &, int, int))
{
    double PL = 0.0, PR = 0.0, DL = 0.0, DR = 0.0, UL = 0.0, UR = 0.0, VL = 0.0, VR = 0.0, WL = 0.0, WR = 0.0, HL = 0.0, HR = 0.0, GL = 0.0, GR = 0.0;
    double CL = 0.0, CR = 0.0, MaL = 0.0,  MaR = 0.0, B1 = 0.0, B2 = 0.0, G1 = 0.0, G2 = 0.0;
    double CI = 0.0, MaI = 0.0, VI = 0.0;
    // 空间重构
    Diff(direction, xnode, ynode, znode, PLR, P, bc, 0);
    Diff(direction, xnode, ynode, znode, DLR, D, bc, 0);
    Diff(direction, xnode, ynode, znode, ULR, U, bc, 0);
    Diff(direction, xnode, ynode, znode, VLR, V, bc, 0);
    Diff(direction, xnode, ynode, znode, WLR, W, bc, 0);
    Diff(direction, xnode, ynode, znode, HLR, H, bc, 0);
    Diff(direction, xnode, ynode, znode, GLR, Gamma, bc, 0);

    // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
    for(int s = 0; s < NS; s++)
    {
        for(int i = 0; i < ni + 2 * bc; i++)
            for(int j = 0; j < nj + 2 * bc; j++)
                for(int k = 0; k < nk + 2 * bc; k++)
                    Yi_temp(i, j, k) = Yi(i, j, k, s);

        Diff(direction, xnode, ynode, znode, YLR, Yi_temp, bc, 0);
        for(int i = 1; i < ni; i++)
            for(int j= 1; j < nj + bc; j++)
                for(int k = 1; k < nk + bc; k++)
                {
                    YL_temp(i, j, k, s) = YLR(i, j, k, 0);
                    YR_temp(i, j, k, s) = YLR(i, j, k, 1);
                }
    }

    // #pragma omp parallel for num_threads(num_thread) collapse(3) schedule(static)
    // #pragma omp simd
    for(int i = 1; i < ni + bc; i++)
        for(int j = 1; j < nj + bc; j++)
            for(int k = 1; k < nk + bc; k++)
            {
                // double PL = 0.0, PR = 0.0, DL = 0.0, DR = 0.0, UL = 0.0, UR = 0.0, VL = 0.0, VR = 0.0, WL = 0.0, WR = 0.0, HL = 0.0, HR = 0.0, GL = 0.0, GR = 0.0;
                // double CL = 0.0, CR = 0.0, MaL = 0.0,  MaR = 0.0, B1 = 0.0, B2 = 0.0, G1 = 0.0, G2 = 0.0;
                // double CI = 0.0, MaI = 0.0, VI = 0.0;
                // Array<double, 1> YL(NS);
                // Array<double, 1> YR(NS);
                PL = PLR(i, j, k, 0);
                PR = PLR(i, j, k, 1);
                DL = DLR(i, j, k, 0);
                DR = DLR(i, j, k, 1);
                UL = ULR(i, j, k, 0);
                UR = ULR(i, j, k, 1);
                VL = VLR(i, j, k, 0);
                VR = VLR(i, j, k, 1);
                WL = WLR(i, j, k, 0);
                WR = WLR(i, j, k, 1);
                HL = HLR(i, j, k, 0);
                HR = HLR(i, j, k, 1);
                GL = GLR(i, j, k, 0);
                GR = GLR(i, j, k, 1);

                for(int s = 0; s < NS; s++)
                {
                    YL(s) = YL_temp(i, j, k, s);
                    YR(s) = YR_temp(i, j, k, s);
                }

                CL = sqrt(GL * PL / DL);         // 计算左侧声速
                CR = sqrt(GR * PR / DR);         // 计算右侧声速
                CI = 0.5 * (CL + CR);            // 界面参考声速

                if(direction == 1)
                {
                    MaL = UL / CI;               // 左马赫数
                    MaR = UR / CI;               // 右马赫数
                }
                else if (direction == 2)
                {
                    MaL = VL / CI;
                    MaR = VR / CI;
                }
                else if(direction == 3)
                {
                    MaL = WL / CI;
                    MaR = WR / CI;
                }

                if(abs(MaL) < 1.0)               // 使用多项式进行平滑分裂，确保音速点附近没有数值间断
                {
                    B1 = 0.25 * pow((MaL + 1.0), 2) + 0.125 * pow((pow(MaL, 2) - 1.0), 2);
                    G1 = 0.25 * pow((MaL + 1.0), 2) * (2.0 - MaL) + 3.0 / 16.0 * MaL * pow((pow(MaL, 2) - 1.0), 2);
                }
                else                             // 使用全上风处理
                {
                    B1 = 0.5 * (MaL + abs(MaL));
                    G1 = B1 / MaL;
                }
                if(abs(MaR) < 1.0)
                {
                    B2 = -0.25 * pow((MaR - 1.0), 2) - 0.125 * pow((pow(MaR, 2) - 1.0), 2);
                    G2 = 0.25 * pow((MaR - 1.0), 2) * (2.0 + MaR) - 3.0 / 16.0 * MaR * pow((pow(MaR, 2) - 1.0), 2);
                }
                else
                {
                    B2 = 0.5 * (MaR - abs(MaR));
                    G2 = B2 / MaR;
                }
                // 通量组装
                MaI = B1 + B2;                   // 界面马赫数
                VI = MaI * CI;                   // 界面速度

                if (direction == 1)
                {
                    if (VI >= 0.0)                // 根据 VI 的正负决定使用哪一侧的物理量进行对流更新
                    {
                        for (int s = 0; s < NS; s++)
                            Fi(i, j, k, s) = VI * DL * YL(s);
                        Fi(i, j, k, NS) = VI * DL * UL + G1 * PL + G2 * PR;
                        Fi(i, j, k, NS + 1) = VI * DL * VL;
                        Fi(i, j, k, NS + 2) = VI * DL * WL;
                        Fi(i, j, k, NS + 3) = VI * (DL * HL + 0.5 * DL * (pow(UL, 2) + pow(VL, 2) + pow(WL, 2)));
                    }
                    else
                    {
                        for (int s = 0; s < NS; s++)
                            Fi(i, j, k, s) = VI * DR * YR(s);
                        Fi(i, j, k, NS) = VI * DR * UR + G1 * PL + G2 * PR;
                        Fi(i, j, k, NS + 1) = VI * DR * VR;
                        Fi(i, j, k, NS + 2) = VI * DR * WR;
                        Fi(i, j, k, NS + 3) = VI * (DR * HR + 0.5 * DR * (pow(UR, 2) + pow(VR, 2) + pow(WR, 2)));
                    }
                }
                else if (direction == 2)
                {
                    if (VI >= 0.0)
                    {
                        for (int s = 0; s < NS; s++)
                            Fi(i, j, k, s) = VI * DL * YL(s);
                        Fi(i, j, k, NS) = VI * DL * UL;
                        Fi(i, j, k, NS + 1) = VI * DL * VL + G1 * PL + G2 * PR;
                        Fi(i, j, k, NS + 2) = VI * DL * WL;
                        Fi(i, j, k, NS + 3) = VI * (DL * HL + 0.5 * DL * (pow(UL, 2) + pow(VL, 2) + pow(WL, 2)));
                    }
                    else
                    {
                        for(int s = 0; s < NS; s++)
                            Fi(i, j, k ,s) = VI * DR * YR(s);
                        Fi(i, j, k, NS) = VI * DR * UR;
                        Fi(i, j, k, NS + 1) = VI * DR * VR + G1 * PL + G2 * PR;
                        Fi(i, j, k, NS + 2) = VI * DR * WR;
                        Fi(i, j, k, NS + 3) = VI * (DR * HR + 0.5 * DR * (pow(UR, 2) + pow(VR, 2) + pow(WR, 2)));
                    }
                }
                else if(direction == 3)
                {
                    if(VI >= 0.0)
                    {
                        for (int s = 0; s < NS; s++)
                            Fi(i, j, k, s) = VI * DL * YL(s);
                        Fi(i, j, k, NS) = VI * DL * UL;
                        Fi(i, j, k, NS + 1) = VI * DL * VL;
                        Fi(i, j, k, NS + 2) = VI * DL * WL + G1 * PL + G2 * PR;
                        Fi(i, j, k, NS + 3) = VI * (DL * HL + 0.5 * DL * (pow(UL, 2) + pow(VL, 2) + pow(WL, 2)));
                    }
                    else
                    {
                        for(int s = 0; s < NS; s++)
                            Fi(i, j, k, s) = VI * DR * YR(s);
                        Fi(i, j, k, NS) = VI * DR * UR;
                        Fi(i, j, k, NS + 1) = VI * DR * VR;
                        Fi(i, j, k, NS + 2) = VI * DR * WR + G1 * PL + G2 * PR;
                        Fi(i, j, k, NS + 3) = VI * (DR * HR + 0.5 * DR * (pow(UR, 2) + pow(VR, 2) + pow(WR, 2)));
                    }
                }
            }
}











// 具备从内能和组分反算温度的能力
double Flowfield::Get_temp(double T, int i, int j, int k)
{
    double T0 = T, T_temp = T;
    double temp1 = 0.0, temp2 = 0.0;
    int count = 0;
    Array<double, 1> Cpi0, Hi0;
    Cpi0.Initial(NS);
    Hi0.Initial(NS);
    while(count < 10)
    {
        T0 = T_temp;
        for(int s = 0; s < NS; s++)
        {
            Cpi0(s) = React.GetCpi(T0, Ri(s), s, Coeff0, Coeff1);
            Hi0(s) = React.GetHi(T0, Ri(s), s, Coeff0, Coeff1);
        }
        temp1 = (Fun.sum(2, Yi_temp0, Hi) - E(i, j, k) * 1e-3) - Fun.sum(2, Yi_temp0, Ri) * T0;
        temp2 = Fun.sum(2, Yi_temp0, Cpi) - Fun.sum(2, Yi_temp0, Ri);
        T_temp = T0 - temp1 / temp2;
        if(abs(T_temp - T0) < 1e-6)
            break;
        else
            count = count + 1;
    }

    return T_temp;
}
double Flowfield::Get_temp(double T, int i, int j, int k, Array<double, 1> &Yi_temp0)
{
    double T0 = T, T_temp = T;
    double temp1 = 0.0, temp2 = 0.0;
    int count = 0;
    Array<double, 1> Cpi0(NS), Hi0(NS);
    while(count < 10)
    {
        T0 = T_temp;
        for(int s = 0; s < NS; s++)
        {
            Cpi0(s) = React.GetCpi(T0, Ri(s), s, Coeff0, Coeff1);
            Hi0(s) = React.GetHi(T0, Ri(s), s, Coeff0, Coeff1); 
        }
        temp1 = (Fun.sum(2, Yi_temp0, Hi0) - E(i, j, k) * 1e-3) - Fun.sum(2, Yi_temp0, Ri) * T0;
        temp2 = Fun.sum(2, Yi_temp0, Cpi0) - Fun.sum(2, Yi_temp0, Ri);
        T_temp = T0 - temp1 / temp2;
        if(abs(T_temp - T0) < 1e-6)
            break;
        else
            count = count + 1;
    }

    return T_temp;
}











void Flowfield::GetPartial_T()
{
    double CV = 0.0, Hi = 0.0;

    for(int i = bc; i < ni + bc; i++)
        for(int j = bc; j < nj + bc; j++)
            for(int k = bc; k < nk + bc; k++)
            {
                CV = Cp(i, j, k) - Rgas(i, j, k);
                for(int s = 0; s < NS; s++)
                {
                    Hi = React.GetHi(T(i, j, k), Ri(s), s, Coeff0, Coeff1);
                    Partion_T(i, j, k, s) = (0.5 * (pow(U(i, j, k), 2) + pow(V(i, j, k), 2) + pow(W(i, j, k), 2)) - Hi + Ri(s) * 1e3 * T(i, j, k)) / (CV * D(i, j, k));
                }
                Partion_T(i, j, k, NS + 1) = -U(i, j, k) / (CV * D(i, j, k));
                Partion_T(i, j, k, NS + 2) = 1.0 / (CV * D(i, j, k));
            }
}

void Flowfield::PackagePrev(std::vector<double> &data, int i, int j, int k)
{
    for (int s = 0; s < NS; s++)
        data.push_back(Mc(i, j, k, s));
    for (int s = 0; s < NS; s++)
        data.push_back(Di(i, j, k, s));
    for (int s = 0; s < NS; s++)
        data.push_back(Yi(i, j, k ,s));
    data.push_back(T(i, j, k));
    data.push_back(E(i, j, k));
    data.push_back(U(i, j, k));
    data.push_back(V(i, j, k));
    data.push_back(W(i, j, k));
}
void Flowfield::UnpackagePrev(std::vector<double> &data, int meshnum)
{
    int index = 0;
    for (int i = 0; i < meshnum; i++)
        for (int s = 0; s < NS; s++)
        {
            Mc(i, 0, 0, s) = data[index++];
            for (int s = 0; s < NS; s++)
                Di(i, 0, 0, s) = data[index++];
            for (int s = 0; s < NS; s++)
                Yi(i, 0, 0, s) = data[index++];
            T(i, 0, 0) = data[index++];

            E(i, 0, 0) = data[index++];
            U(i, 0, 0) = data[index++];
            V(i, 0, 0) = data[index++];
            W(i, 0, 0) = data[index++];
        }

}
void Flowfield::UnpackagePrev(std::vector<double> &data, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    for (int n = 0; n < recvFrom.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin() + recvFrom[0], transferNmesh.begin() + recvFrom[n], 0) * (5 + 3 * NS);
        size = transferNmesh[recvFrom[n]];
        for (int i = 0; i < size; index++, i++)
        {
            for (int s = 0; s < NS; s++)
                Mc(index, 0, 0, s) = data[start++];
            for (int s = 0; s < NS; s++)
                Di(index, 0, 0, s) = data[start++];
            for (int s = 0; s < NS; s++)
                Yi(index, 0, 0, s) = data[start++];
            T(index, 0, 0) = data[start++];
            E(index, 0, 0) = data[start++];
            U(index, 0, 0) = data[start++];
            V(index, 0, 0) = data[start++];
            W(index, 0, 0) = data[start++];
        }
    }
}
void Flowfield::UnpackagePrev(std::vector<double> &data, int meshnum, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{

    int index = 0, start = 0, size = 0;
    for (int n = 0; n < recvFrom.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin(), transferNmesh.begin() + recvFrom[n], 0) * (5 + 3 * NS);
        size = transferNmesh[recvFrom[n]];

        for (int i = 0; i < size; index++, i++)
        {
            for (int s = 0; s < NS; s++)
                Mc(index, 0, 0, s) = data[start++];
            for (int s = 0; s < NS; s++)
                Di(index, 0, 0, s) = data[start++];
            for (int s = 0; s < NS; s++)
                Yi(index, 0, 0, s) = data[start++];
            T(index, 0, 0) = data[start++];
            E(index, 0, 0) = data[start++];
            U(index, 0, 0) = data[start++];
            V(index, 0, 0) = data[start++];
            W(index, 0, 0) = data[start++];
        }
    }


}
void Flowfield::UnpackagePrev(std::vector<double> &data, Array<int, 1> &mesh, Array<int, 1> &NchemIndex)
{
    int start = 0, index = 0, size = 0, id = 0;
    for (int i = 0; i < numprocs; i++)
        for (int j = 0; j < numprocs; j++)
        {
            id = NchemIndex(numprocs - i - 1) * numprocs + j;
            if (mesh(id) > 0)
            {
                if (id % numprocs == myid)
                {
                    index = start * (5 + 3 * NS);
                    for (int k = 0; k < mesh(id) * (5 + 3 * NS); k++)
                    {
                        for (int s = 0; s < NS; s++)
                            Mc(index, 0, 0, s) = data[start++];
                        for (int s = 0; s < NS; s++)
                            Di(index, 0, 0, s) = data[start++];
                        for (int s = 0; s < NS; s++)
                            Yi(index, 0, 0, s) = data[start++];
                        T(k + size, 0, 0) = data[index++];
                        E(k + size, 0, 0) = data[index++];
                        U(k + size, 0, 0) = data[index++];
                        V(k + size, 0, 0) = data[index++];
                        W(k + size, 0, 0) = data[index++];
                    }
                    size += mesh(id);
                }
                start = mesh(id);
            }
        }


}

void Flowfield::PackageUpdate(std::vector<double> &data, int meshnum)
{
    for(int i = 0; i < meshnum; i++)
    {
        data.push_back(D(i, 0, 0));
        data.push_back(T(i, 0, 0));

        data.push_back(U(i, 0, 0));
        data.push_back(V(i, 0, 0));
        data.push_back(W(i, 0, 0));
        data.push_back(E(i, 0, 0));
        data.push_back(Wav(i, 0, 0));
        data.push_back(Rgas(i, 0, 0));
        data.push_back(P(i, 0, 0));
        data.push_back(H(i, 0, 0));
        data.push_back(Cp(i, 0, 0));
        data.push_back(Gamma(i, 0, 0));
        data.push_back(C(i, 0, 0));
        data.push_back(Ma(i, 0, 0));

        for(int s = 0; s < NS; s++)
            data.push_back(CS(i, 0, 0, s));
        for(int s = 0; s < NS; s++)
            data.push_back(Di(i, 0, 0, s));
        for(int s = 0; s < NS; s++)
            data.push_back(Mc(i, 0, 0, s));
        for(int s = 0; s < NS; s++)
            data.push_back(Yi(i, 0, 0, s));
    }
}
void Flowfield::PackageUpdate(std::vector<double> &data, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    for(int n = 0; n < recvFrom.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin() + recvFrom[0], transferNmesh.begin() + recvFrom[n], 0) * (18 + 4 * NS);
        size = transferNmesh[recvFrom[n]];

        for(int i = 0; i < size; index++, i++)
        {
            data[start++] = D(index, 0, 0);
            data[start++] = T(index, 0, 0);
            data[start++] = U(index, 0, 0);
            data[start++] = V(index, 0, 0);
            data[start++] = W(index, 0, 0);
            data[start++] = E(index, 0, 0);
            data[start++] = Wav(index, 0, 0);
            data[start++] = Rgas(index, 0, 0);
            data[start++] = P(index, 0, 0);
            data[start++] = H(index, 0, 0);
            data[start++] = Cp(index, 0, 0);
            data[start++] = Gamma(index, 0, 0);
            data[start++] = C(index, 0, 0);
            data[start++] = Ma(index, 0, 0);

            for(int s = 0; s < NS + 4; s++)
                data[start++] = CS(index, 0, 0, s);
            for(int s = 0; s < NS; s++)
                data[start++] = Di(index, 0, 0, s);
            for(int s = 0; s < NS; s++)
                data[start++] = Mc(index, 0, 0, s);
            for(int s = 0; s < NS; s++)
                data[start++] = Yi(index, 0, 0, s);
        }
    }
}
void Flowfield::PackageUpdate(std::vector<double> &data, int meshnum, std::vector<int> &recvFrom, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    for(int n = 0; n < recvFrom.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin(), transferNmesh.begin() + recvFrom[n], 0) * (18 + 4 * NS);
        size = transferNmesh[recvFrom[n]];

        for(int i = 0; i < size; index++, i++)
        {
            data[start++] = D(index, 0, 0);
            data[start++] = T(index, 0, 0);
            data[start++] = U(index, 0, 0);
            data[start++] = V(index, 0, 0);
            data[start++] = W(index, 0, 0);
            data[start++] = E(index, 0, 0);
            data[start++] = Wav(index, 0, 0);
            data[start++] = Rgas(index, 0, 0);
            data[start++] = P(index, 0, 0);
            data[start++] = H(index, 0, 0);
            data[start++] = Cp(index, 0, 0);
            data[start++] = Gamma(index, 0, 0);
            data[start++] = C(index, 0, 0);
            data[start++] = Ma(index, 0, 0);

            for(int s = 0; s < NS + 4; s++)
                data[start++] = CS(index, 0, 0, s);
            for(int s = 0; s < NS; s++)
                data[start++] = Di(index, 0, 0, s);
            for(int s = 0; s < NS; s++)
                data[start++] = Mc(index, 0, 0, s);
            for(int s = 0; s < NS; s++)
                data[start++] = Yi(index, 0, 0, s);
        }
    }
}

void Flowfield::UnpackageUpdate(std::vector<double> &data, std::vector<int> &transferIndex, int meshnum, std::vector<int> &sendTo, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    int x = 0, y = 0, z = 0;
    for(int n = 0; n < sendTo.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin(), transferNmesh.begin() + sendTo[n], 0) * (18 + 4 * NS);
        size = transferNmesh[sendTo[n]];
        for(int i = 0; i < size; index++, i++)
        {
            x = transferIndex[3 * index];
            y = transferIndex[3 * index + 1];
            z = transferIndex[3 * index + 2];
            D(x, y, z) = data[start++];
            T(x, y, z) = data[start++];
            U(x, y, z) = data[start++];
            V(x, y, z) = data[start++];
            W(x, y, z) = data[start++];
            E(x, y, z) = data[start++];
            Wav(x, y, z) = data[start++];
            Rgas(x, y, z) = data[start++];
            P(x, y, z) = data[start++];
            H(x, y, z) = data[start++];
            Cp(x, y, z) = data[start++];
            Gamma(x, y, z) = data[start++];
            C(x, y, z) = data[start++];
            Ma(x, y, z) = data[start++];

            for(int s = 0; s < NS + 4; s++)
                CS(x, y, z, s) = data[start++];
            for(int s = 0; s < NS; s++)
                Di(x, y, z, s) = data[start++];
            for(int s = 0; s < NS; s++)
                Mc(x, y, z, s) = data[start++];
            for(int s = 0; s < NS; s++)
                Yi(x, y, z, s) = data[start++];
        }
    }
}

void Flowfield::UnpackageUpdate(std::vector<double> &data, std::vector<int> &transferIndex, std::vector<int> &sendTo, std::vector<int> &transferNmesh)
{
    int index = 0, start = 0, size = 0;
    int x = 0, y = 0, z = 0;
    for(int n = 0; n < sendTo.size(); n++)
    {
        start = std::accumulate(transferNmesh.begin() + sendTo[0], transferNmesh.begin() + sendTo[n], 0) * (18 + 4 * NS);
        size = transferNmesh[sendTo[n]];
        for(int i = 0; i < size; index++, i++)
        {
            x = transferIndex[3 * index];
            y = transferIndex[3 * index + 1];
            z = transferIndex[3 * index + 2];
            D(x, y, z) = data[start++];
            T(x, y, z) = data[start++];
            U(x, y, z) = data[start++];
            V(x, y, z) = data[start++];
            W(x, y, z) = data[start++];
            E(x, y, z) = data[start++];
            Wav(x, y, z) = data[start++];
            Rgas(x, y, z) = data[start++];
            P(x, y, z) = data[start++];
            H(x, y, z) = data[start++];
            Cp(x, y, z) = data[start++];
            Gamma(x, y, z) = data[start++];
            C(x, y, z) = data[start++];
            Ma(x, y, z) = data[start++];


            for(int s = 0; s < NS + 4; s++)
                CS(x, y, z, s) = data[start++];
            for(int s = 0; s < NS; s++)
                Di(x, y, z, s) = data[start++];
            for(int s = 0; s < NS; s++)
                Mc(x, y, z, s) = data[start++];
            for(int s = 0; s < NS; s++)
                Yi(x, y, z, s) = data[start++];
        }
    }
}
// 实现了高阶空间离散，用于提高计算精度
void Flowfield::ReConstruction(int meshnum)
{
    U.Initial(meshnum, 1, 1);
    V.Initial(meshnum, 1, 1);
    W.Initial(meshnum, 1, 1);
    P.Initial(meshnum, 1, 1);
    D.Initial(meshnum, 1, 1);
    T.Initial(meshnum, 1, 1);
    C.Initial(meshnum, 1, 1);
    Ma.Initial(meshnum, 1, 1);
    Wav.Initial(meshnum, 1, 1);
    Rgas.Initial(meshnum, 1, 1);
    Cp.Initial(meshnum, 1, 1);
    H.Initial(meshnum, 1, 1);
    E.Initial(meshnum, 1, 1);
    Gamma.Initial(meshnum, 1, 1);

    Mr.Initial(meshnum, 1, 1, NS);
    Mc.Initial(meshnum, 1, 1, NS);
    Mi.Initial(meshnum, 1, 1, NS);
    Yi.Initial(meshnum, 1, 1, NS);
    Di.Initial(meshnum, 1, 1, NS);
    CS.Initial(meshnum, 1, 1, NS + 4);
}
