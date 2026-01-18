







# include <cstring>
# include <iostream>
# include <fstream>
# include <sstream>
# include <string>
# include <mpi.h>
# include <omp.h>
# include <vector>
# include <numeric>
# include <filesystem>





# include "Euler.hpp"
# include "_string_switch.hpp"
# include "FileReader.hpp"

using namespace std;













// 构造函数
Euler::Euler()
{

    cout << "Input ctrl file is:" << this->ctrlfile << endl;

    // t1 - t2 读取配置文件
    double t1 = MPI_Wtime();
    FileRead();
    double t2 = MPI_Wtime();
    // t2 - t3 并行环境初始化
    Mpi_Initial();
    double t3 = MPI_Wtime();
    // t3 - t4 读取网格文件，建立网格节点的坐标系统，计算网格单元的体积、面积及法向量
    Mymesh.MeshProcess(meshFile);
    double t4 = MPI_Wtime();
    // t4 - t5 读取输入数据
    TwoDim.InputRead(inputFile, Mymesh.xnode, Mymesh.ynode, Mymesh.znode, Mymesh.bc);
    ExtraFlow.InputRead(inputFile, Mymesh.xnode, Mymesh.ynode, Mymesh.znode, Mymesh.bc);
    double t5 = MPI_Wtime();
    // t5 - t6 读取化学反应机理和 NASA 热力学多项式系数
    React.ReactionRead(Reaction_Model, thermoFile, Mymesh.xnode, Mymesh.ynode, Mymesh.znode, Mymesh.bc);
    ExtraReact.ReactionRead(Reaction_Model, thermoFile, Mymesh.xnode, Mymesh.ynode, Mymesh.znode, Mymesh.bc);
    double t6 = MPI_Wtime();
    // t6 - t7 流场初始化
    TwoDim.FieldInitial(React.Ri, React.Mw, React.Mi, React.Yi, React.Coeff0, React.Coeff1);
    ExtraFlow.FieldInitial(React.Ri, React.Mw, React.Mi, React.Yi, React.Coeff0, React.Coeff1);
    double t7 = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);                                     // 进程同步
    double t8 = MPI_Wtime();

    Computing(count);                                                // 核心计算循环
    double t9 = MPI_Wtime();

    Output(count);                                                   // 结果输出与耗时统计
    double t10 = MPI_Wtime();
    Output_Total_Time(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
    return;
}










// 读取配置文件并创建输出文件目录
void Euler::FileRead()
{
    FileReader Ctrl;
    Ctrl.readFile(ctrlfile);

    // 赋值配置文件中设置的参数
    Final_Time = Ctrl.getdoubleParameter("Final-Time");
    cfl = Ctrl.getdoubleParameter("cfl");
    if(cfl == 0)
        dt = Ctrl.getdoubleParameter("dt");
    count = Ctrl.getIntParameter("count");


    strcpy(meshFile, Ctrl.getStringParameter("gridname").c_str());
    strcpy(inputFile, Ctrl.getStringParameter("initialization").c_str());
    strcpy(Reaction_Model, Ctrl.getStringParameter("reaction-model").c_str());
    strcpy(thermoFile, Ctrl.getStringParameter("Thermofile").c_str());


    strcpy(react, Ctrl.getStringParameter("reaction").c_str());
    switch(hash_(react))
    {
    case "Trapezoid"_hash:
        Reaction_Sch = Global::Trapezoid;
        type = 0;
        break;
    case "Trapezoid_Adaptive"_hash:
        Reaction_Sch = Global::Trapezoid;
        type = 1;
        break;
    case "IMEX"_hash:
        Reaction_Sch = Global::IMEX;
        break;
    case "Cantera"_hash:
        Reaction_Sch = Global::DNN;
        type = 0;
        break;
    case "DNN"_hash:
        Reaction_Sch = Global::DNN;
        type = 1;
        break;
    default:
        cout << "Invalid reaction scheme" << endl;
        abort();
    }


    strcpy(diff, Ctrl.getStringParameter("difference").c_str());
    switch(hash_(diff))
    {
    case "MUSCL_1"_hash:
        Diff_Sch = Global::MUSCL_1;
        break;
    case "MUSCL_2"_hash:
        Diff_Sch = Global::MUSCL_2;
        break;
    case "WENO"_hash:
        Diff_Sch = Global::WENO;
        break;
    default:
        cout << "Invalid difference scheme" << endl;
    }


    strcpy(timeadv, Ctrl.getStringParameter("time-adv").c_str());
    switch(hash_(timeadv))
    {
    case "EE"_hash:
        TimeAdv_Sch = Global::EE;
        break;
    case "TVD_RK3"_hash:
        TimeAdv_Sch = Global::TVD_RK3;
        break;
    default:
        cout << "Invalid time-advance scheme" << endl;
    }


    m_block_x = Ctrl.getIntParameter("m_block_x");
    m_block_y = Ctrl.getIntParameter("m_block_y");
    m_block_z = Ctrl.getIntParameter("m_block_z");

    
    num_thread = Ctrl.getIntParameter("num_thread");
    

    DLB_step = Ctrl.getIntParameter("DLB_step");
    DLB_tol = Ctrl.getIntParameter("DLB_tol");


    if (myid == 0)
    {
        std::string outputFolder = "./output_" + std::to_string(numprocs);
        std::filesystem::create_directory(outputFolder);

        std::filesystem::create_directory(outputFolder + "/grid");
        std::filesystem::create_directory(outputFolder + "/load");
        std::filesystem::create_directory(outputFolder + "/time_record");
        std::filesystem::create_directory(outputFolder + "/time_record/all");
        std::filesystem::create_directory(outputFolder + "/time_record/compare");
        std::filesystem::create_directory(outputFolder + "/time_record/compute");
        std::filesystem::create_directory(outputFolder + "/time_record/total");
        std::filesystem::create_directory(outputFolder + "/results");

        std::cout << "目录结构创建成功：" << outputFolder << std::endl;
    }

    return;
}









// 并行环境初始化
void Euler::Mpi_Initial()
{
    myid_x = int(myid % m_block_x);
    myid_y = int(floor(myid % (m_block_x * m_block_y) / m_block_x));
    myid_z = int(floor(myid / (m_block_x * m_block_y)));
    if(myid_x > 0)                                                   // 在 x 方向设置进程的左邻居进程
        m_left = myid - 1;
    else
        m_left = MPI_PROC_NULL;
    if(myid_x < m_block_x - 1)                                       // 在 x 方向设置进程的右邻居进程
        m_right = myid + 1;
    else
        m_right = MPI_PROC_NULL;
    
    if(myid_y > 0)                                                   // 在 y 方向设置进程的前邻居进程
        m_front = myid - m_block_x;
    else
        m_front = MPI_PROC_NULL;
    if (myid_y < m_block_y - 1)                                      // 在 y 方向设置进程的后邻居进程
        m_back = myid + m_block_x;
    else
        m_back = MPI_PROC_NULL;

    if (myid_z > 0)                                                   // 在 z 方向设置进程的下邻居进程
        m_down = myid - m_block_x * m_block_y;
    else
        m_down = MPI_PROC_NULL;
    if(myid_z < m_block_z - 1)                                       // 在 z 方向设置进程的上邻居进程
        m_up = myid + m_block_x * m_block_y;
    else
        m_up = MPI_PROC_NULL;
    // cout << myid << " " << myid_x << " " << myid_y << " " << myid_z << " " << m_left << " " << m_right << " " << m_front << " " << m_back << " " << m_down << " " << m_up << endl;
}












void Euler::Computing(int count)
{
    double time = 0.0;                                               // 当前时间
    int num = 0;
    iteration = 0;                                                   // 最大迭代次数
    double t1 = 0,  t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
    double t12 = 0, t23 = 0, t34 = 0, t45 = 0, t67 = 0;
    double tt1 = 0, tt2 = 0;
    ofstream loadfile("./output_" + to_string(numprocs) + "/load/LoadDegree.dat");
    tt1 = MPI_Wtime();
    while (time < Final_Time)                                        // 时间还没到
    {

        if (cfl == 0)
            TwoDim.dt = dt;                                          // cfl 等于 0 使用用户指定固定步长
        else
            TwoDim.CFLcondition(cfl, Final_Time);                    // cfl 不等于 0 时，为防止计算发散，根据当前声速和流速，动态计算安全的时间步长，

        // 根据总物理时间 Final_Time 和预设的输出文件总数 count，判断当前时刻是否跨越了输出节点。如果是，则调用 Output(num) 保存 .dat 结果文件
        if (time <= num * Final_Time / count && time + TwoDim.dt > num * Final_Time / count)
        {
            t6 = MPI_Wtime();
            Output(num);
            num = num + 1;
            t7 = MPI_Wtime();
            t67 += t7 - t6;                                          // t67 记录输出 dat 文件操作需要的时间
        }
        t1 = MPI_Wtime();

        TwoDim.FieldBoundary_3dRSBI();                               // 设置六个面全对称物理边界，激波发生全反射
        t2 = MPI_Wtime();
        t12 += t2 - t1;
        TwoDim.Mpi_Boundary();                                       // 实现相邻进程之间的数据块交换
        t3 = MPI_Wtime();
        t23 += t3 - t2;

        TwoDim.Advection(TimeAdv_Sch, Diff_Sch);                     // 求解欧拉方程中的通量项（对流部分）


        t4 = MPI_Wtime();
        t34 += t4 - t3;
        
        switch (Reaction_Sch)                                         // 化学反应源项
        {
        case 0:
            Trapezoid(type);                                         // 简单的二阶梯形公式求解
            break;
        case 1:
            IMEX();                                                  // 显隐式混合法，用于处理化学反应中极其“刚性”的时间尺度问题
            break;
        case 2:
            DNN(type);                                               // 神经网络辅助求解反应源项
            break;
        }
        t5 = MPI_Wtime();
        t45 += t5 - t4;

        time += TwoDim.dt;
        iteration += 1;
        if (myid == 0)
            cout << "Process " << myid << ": Iteration: " << iteration << " " << TwoDim.dt << " " << time << endl;
        
        // 数值稳定性检查：如果密度或温度变成了 NaN，意味着由于物理模型不合理或时间步长过大，计算已经在数学上彻底崩溃了。此时继续计算没有任何意义。
        if (TwoDim.D.IsNan())
        {
            cout << "D is nan\n";
            Output(999);
            abort();
        }

        if (TwoDim.T.IsNan())
        {
            cout << "T is nan\n";
            Output(999);
            abort();
        }
        loadfile << iteration << '\t' << DLB.LoadDegree << '\t' << DLB.LoadDegreeBalance << '\n';
    }
    MPI_Barrier(MPI_COMM_WORLD);
    tt2 = MPI_Wtime();
    trackt = tt2 - tt1;
    Output_Compute_Time(t12, t23, t34, t45, t67);


    if (Reaction_Sch == 0 && type == 1)                              // 在开启了梯形法 + 负载均衡时，统计并更新全场最大的化学反应子步数。
        React.GetNchemMax();

    return;
}










// 求解化学反应的常微分方程，集成了动态负载均衡逻辑
// 在处理具有强刚性的化学反应计算时，通过 MPI 进程间的数据迁移，平衡各核心的计算压力，从而显著提升大规模并行计算效率
void Euler::Trapezoid(int type)
{
    int ni = Mymesh.ni;
    int nj = Mymesh.nj;
    int nk = Mymesh.nk;
    int bc = Mymesh.bc;

    double t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0, t6 = 0, t7 = 0;
    double tw1 = 0, tw2 = 0, tw3 = 0, tw4 = 0, tw5 = 0, tw6 = 0, tw7 = 0;

    t1 = MPI_Wtime();

    TwoDim.Update_after_Adv();                                       // 显式时间推进，更新保守变量和原始变量，为反应项计算做准备
    t2 = MPI_Wtime();
    tr1 += t2 - t1;
    if (type == 0)                                                   // 标准模式
    {
        // 采用半隐式梯形法求解化学动力学常微分方程，计算每个网格点的反应源项
        TwoDim.Di = React.Trapezoid(TwoDim.Mc, TwoDim.Di, TwoDim.T, TwoDim.dt);

        // 将原始变量（如温度、压力）转换成求解器计算需要的保守变量（密度、动量、总能、组分密度）
        TwoDim.Explicit();
    }
    else if (type == 1)                                              // 负载均衡模式
    {   // 先预测每个网格点执行化学反应所需的计算量
        double ttp1 = MPI_Wtime();
        React.TrapezoidPrection(TwoDim.Mc, TwoDim.Di, TwoDim.Yi, TwoDim.T, TwoDim.dt);
        double ttp2 = MPI_Wtime();
        dt4 += ttp2 - ttp1;

        int nchemSumBalance = 0;
        // 负载均衡检查
        if (iteration % DLB_step == 0)                               // 每隔 DLB_step 步重新执行一次负载均衡检查
        {
            DLB.NchemSumBalance = 0;
            React.NchemSum = React.Nchem.SumNoBoundary(bc);          // 本进程内部所有数据点的化学计算步数总和
            MPI_Allgather(&React.NchemSum, 1, MPI_INT, &React.NchemTotal(0), 1, MPI_INT, MPI_COMM_WORLD);
            DLB.LoadDegree = double((React.NchemTotal.MaxValue() - React.NchemTotal.AveValue())) / React.NchemTotal.MaxValue();
            if (myid == 0)
                cout << "Load Degree: " << DLB.LoadDegree << endl;
            if(DLB.LoadDegree > DLB_tol)                             // 如果负载不均衡度超过阈值 DLB_tol，则启动 DLBPriorityQueue 来计算哪些网格（Mesh）需要从“忙碌进程”迁移到“空闲进程”
            {
                // if (myid == 0)
                // {
                //     cout << "Before DLB: ";
                //     React.NchemTotal.Print();
                // }
                DLB.transferNchem.clear();                           // 清空迁移记录

                DLB.DLBPriorityQueue(React.NchemTotal);               // 利用优先队列算法，计算哪些网格点该从哪个进程搬到哪个进程
                DLB.transferMeshNum.resize(DLB.transferNchem.size()); // 各种迁移缓冲区（buffer）的 resize 和清空操作，准备接收数据
                DLB.transferNmesh.resize(DLB.transferNchem.size());
                std::fill(DLB.transferMeshNum.begin(), DLB.transferMeshNum.end(), 0);
                std::fill(DLB.transferNmesh.begin(), DLB.transferNmesh.end(), 0);
                DLB.transferPrevDataLocal.clear();
                DLB.transferPrevData.clear();
                DLB.transferPrevNchemLocal.clear();
                DLB.transferPrevNchem.clear();
                DLB.transferIndex.clear();
                DLB.sentNmesh = 0;
                DLB.prevFlowDataLocal.clear();
                DLB.recvPrevReactData.clear();
            }
        }
        MPI_Request sendRequests1[DLB.transferNchem.size()];         // 并行通信准备
        MPI_Request sendRequests2[DLB.transferNchem.size()];
        MPI_Request sendRequests3[DLB.transferNchem.size()];
        MPI_Request recvRequests1[DLB.transferNchem.size()];
        MPI_Request recvRequests2[DLB.transferNchem.size()];
        MPI_Request recvRequests3[DLB.transferNchem.size()];
        MPI_Status statuses[DLB.transferNchem.size()];
        for (int i = 0; i < DLB.transferNchem.size(); i++)
        {
            sendRequests1[i] = MPI_REQUEST_NULL;
            sendRequests2[i] = MPI_REQUEST_NULL;
            sendRequests3[i] = MPI_REQUEST_NULL;
            recvRequests1[i] = MPI_REQUEST_NULL;
            recvRequests2[i] = MPI_REQUEST_NULL;
            recvRequests3[i] = MPI_REQUEST_NULL;
        }
        if (DLB.LoadDegree > DLB_tol)
        {
            double tt1 = MPI_Wtime();

            int sentNchem = DLB.GetSentNchem(myid);                  // 获取本进程需要发出的计算量
            // cout << "\tmyid " << myid << " Nchem: " << sentNchem << endl;
            if(sentNchem != 0)                                       // 本进程是忙碌进程
            {
                DLB.prevFlowDataLocal.clear();
                DLB.prevReactDataLocal.clear();

                std::vector<int> sendToTemp = DLB.sendTo;
                int nchem_temp = 0, id = 0, start = 0, n = 0;


                if (iteration % DLB_step == 0)
                {
                    DLB.GetTransferMesh(React.Nchem, bc);            // 获取待迁移网格的坐标索引
                    DLB.sentNmesh = std::accumulate(DLB.transferMeshNum.begin(), DLB.transferMeshNum.end(), 0);
                }
                int x = 0, y = 0, z = 0;
                for (int i = 0; i < DLB.transferIndex.size() / 3; i++)
                {   // 遍历所有待迁移的网格点索引，打包流场数据和反应数据
                    x = DLB.transferIndex[3 * i];
                    y = DLB.transferIndex[3 * i + 1];
                    x = DLB.transferIndex[3 * i + 2];
                    TwoDim.PackagePrev(DLB.prevFlowDataLocal, x, y, z);
                    React.PackagePrev(DLB.prevReactDataLocal, x, y, z);
                }
            }
            double tt2 = MPI_Wtime();
            dt1 += tt2- tt1;
            if (iteration % DLB_step == 0)                           // 同步迁移规模信息
            {
                MPI_Reduce(DLB.transferMeshNum.data(), DLB.transferNmesh.data(), DLB.transferNmesh.size(), MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
                MPI_Bcast(DLB.transferMeshNum.data(), DLB.transferNmesh.size(), MPI_INT, 0, MPI_COMM_WORLD);
                if(myid == 0)
                {
                    cout << "TransferNmesh: " << std::accumulate(DLB.transferNmesh.begin(), DLB.transferNmesh.end(), 0) << endl;
                }
            }
            double ttt1 = MPI_Wtime();
            if (sentNchem != 0)                                      // 发送数据
                for (int i = 0; i < DLB.sendTo.size(); i++)
                {
                    int start1 = std::accumulate(DLB.transferNmesh.begin() + DLB.sendTo[0], DLB.transferNmesh.begin() + DLB.sendTo[i], 0) * (5 + 3 * TwoDim.NS);
                    int start2 = std::accumulate(DLB.transferNmesh.begin() + DLB.sendTo[0], DLB.transferNmesh.begin() + DLB.sendTo[i], 0);
                    int size1 = DLB.transferNmesh[DLB.sendTo[i]] * (5 + 3 * TwoDim.NS);
                    int size2 = DLB.transferNmesh[DLB.sendTo[i]];
                    if (size2 != 0)
                    {
                        MPI_Isend(&DLB.prevFlowDataLocal[start1], size1, MPI_DOUBLE, std::get<1>(DLB.transferNchem[DLB.sendTo[i]]),
                                  DLB.sendTo[i], MPI_COMM_WORLD, &sendRequests1[DLB.sendTo[i]]);
                        MPI_Isend(&DLB.prevFlowDataLocal[start2], size2, MPI_DOUBLE, std::get<1>(DLB.transferNchem[DLB.sendTo[i]]),
                                  DLB.sendTo[i] + 10000, MPI_COMM_WORLD, &sendRequests2[DLB.sendTo[i]]);
                    }
                }
            double ttt2 = MPI_Wtime();
            int recvNmesh = DLB.GetRecvNmesh(myid);

            if (recvNmesh != 0)                                      // 接收数据
            {
                DLB.recvPrevFlowData.clear();
                DLB.recvPrevReactData.clear();
                DLB.recvPrevFlowData.resize(recvNmesh * (5 + 3 * TwoDim.NS));
                DLB.recvPrevReactData.resize(recvNmesh);
                for (int i = 0; i < DLB.recvFrom.size(); i++)
                {
                    int start1 = std::accumulate(DLB.transferNmesh.begin() + DLB.recvFrom[0], DLB.transferNmesh.begin() + DLB.recvFrom[i], 0) * (5 + 3 * TwoDim.NS);
                    int start2 = std::accumulate(DLB.transferNmesh.begin() + DLB.recvFrom[0], DLB.transferNmesh.begin() + DLB.recvFrom[i], 0);
                    int size1 = DLB.transferNmesh[DLB.recvFrom[i]] * (5 + 3 * TwoDim.NS);
                    int size2 = DLB.transferNmesh[DLB.recvFrom[i]];
                    if (size2 != 0)
                    {




                        MPI_Irecv(&DLB.recvPrevFlowData[start1], size1, MPI_DOUBLE, std::get<1>(DLB.transferNchem[DLB.recvFrom[i]]),
                                  DLB.recvFrom[i], MPI_COMM_WORLD, &recvRequests1[DLB.recvFrom[i]]);
                        MPI_Irecv(&DLB.recvPrevFlowData[start2], size2, MPI_DOUBLE, std::get<1>(DLB.transferNchem[DLB.recvFrom[i]]),
                                  DLB.recvFrom[i] + 10000, MPI_COMM_WORLD, &recvRequests2[DLB.recvFrom[i]]);
                    }
                }
            }
        }
        int tag = 0;                                                 // 本地计算，在等待网络传输的同时，先处理留在本地的网格
        nchemSumBalance = 0;

        for (int i = bc; i < (ni - 1) / 2 + bc; i++)
            for (int j = bc; j < nj + bc; j++)
                for (int k = bc; k < nk + bc; k++)
                {
                    int n = 0;
                    bool flag = false;                               // 检查点是否已被迁出
                    while (n < DLB.sentNmesh && tag < DLB.sentNmesh)
                    {
                        if (i == DLB.transferIndex[3 * n] &&
                            j == DLB.transferIndex[3 * n + 1] &&
                            k == DLB.transferIndex[3 * n + 2])
                        {
                            flag = true;
                            tag++;
                            break;
                        }
                        n++;
                    }
                    if (!flag)
                    {
                        nchemSumBalance += React.Nchem(i, j, k);
                        for (int step = 0; step < React.Nchem(i, j, k); step++)
                        {

                            React.Trapezoid(TwoDim.Mc, TwoDim.Di, TwoDim.Yi, TwoDim.T, TwoDim.dt, step, i, j, k);


                            TwoDim.Explicit(step, i, j, k);
                        }
                    }
                }
            DLB.NchemSumBalance += nchemSumBalance;

            tw1 = MPI_Wtime();
            MPI_Waitall(DLB.transferNchem.size(), sendRequests1, statuses);
            MPI_Waitall(DLB.transferNchem.size(), sendRequests2, statuses);
            MPI_Waitall(DLB.transferNchem.size(), recvRequests1, statuses);
            MPI_Waitall(DLB.transferNchem.size(), recvRequests2, statuses);
            tw2 = MPI_Wtime();
            // 处理收到的来自其他进程的网格点
            if (DLB.LoadDegree > DLB_tol)
            {
                int recvNmesh = DLB.GetRecvNmesh(myid);

                if(recvNmesh != 0)
                {
                    nchemSumBalance = 0;
                    ExtraFlow.ReConstruction(recvNmesh);             // 为客场网格重建临时流场对象
                    ExtraFlow.UnpackagePrev(DLB.recvPrevFlowData, DLB.recvFrom, DLB.transferNmesh);
                    ExtraReact.ReConstruction(recvNmesh);
                    ExtraReact.UnpackagePrev(DLB.recvPrevReactData, DLB.recvFrom, DLB.transferNmesh);
                    // cout << myid << " Start DLB\n";

                    // #pragma omp parallel for num_threads(num_thread) collapse(1) schedule(dynamic) reduction(+ :nchemSumBalance)
                    for (int i = 0; i < recvNmesh; i++)
                    {
                        // DLB.NchemSumBalance += ExtraReact.Nchem(i, 0, 0);
                        nchemSumBalance += ExtraReact.Nchem(i, 0, 0);
                        for (int step = 0; step < ExtraReact.Nchem(i, 0, 0); step ++)
                        {

                            ExtraReact.Trapezoid(ExtraFlow.Mc, ExtraFlow.Di, ExtraFlow.Yi, ExtraFlow.T, TwoDim.dt, step, i, 0, 0);

                            ExtraFlow.Explicit(step, i, 0, 0);
                        }
                    }
                    DLB.NchemSumBalance += nchemSumBalance;
                    double tt4 = MPI_Wtime();
                    // 计算完成后，将结果打包发还给原主人
                    DLB.updateFlowDataLocal.clear();
                    DLB.updateFlowDataLocal.resize(recvNmesh *(18 + 4 * TwoDim.NS));
                    ExtraFlow.PackageUpdate(DLB.updateFlowDataLocal, DLB.recvFrom, DLB.transferNmesh);
                    bool has_nan = std::any_of(DLB.updateFlowDataLocal.begin(), DLB.updateFlowDataLocal.end(), [](double x)
                                                {   return std::isnan(x);});
                    if (has_nan)
                    {
                        std::cout << "DLB contains NaN." << std::endl;
                        abort();
                    }
                    for (int i = 0; i < DLB.recvFrom.size(); i++)
                    {
                        int start = std::accumulate(DLB.transferNmesh.begin() + DLB.recvFrom[0], DLB.transferNmesh.begin() + DLB.recvFrom[i], 0) * (18 + 4 * TwoDim.NS);
                        int size = DLB.transferNmesh[DLB.recvFrom[i]] * (18 + 4 * TwoDim.NS);
                        if (size != 0)
                            MPI_Isend(&DLB.updateFlowDataLocal[start], size, MPI_DOUBLE, std::get<0>(DLB.transferNchem[DLB.recvFrom[i]]),
                                      DLB.recvFrom[i] + 20000, MPI_COMM_WORLD, &sendRequests3[DLB.recvFrom[i]]);
                    }
                }
            }
            t3 = MPI_Wtime();

            t4 = MPI_Wtime();
            tr2 += t4 - t3;
            // 回收那些发送给其他进程代为计算的网格点结果，并将其更新回本地的主流场数组中
            if (DLB.LoadDegree > DLB_tol)
            {
                int sentNchem = DLB.GetSentNchem(myid);              // 发送出去的任务量
                DLB.recvUpdateFlowData.clear();
                DLB.recvUpdateFlowData.resize(DLB.sentNmesh * (18 + 4 * TwoDim.NS));

                if (sentNchem != 0)
                {
                    double tt3 = MPI_Wtime();
                    for (int i = 0; i < DLB.sendTo.size(); i++)
                    {
                        int start = std::accumulate(DLB.transferNmesh.begin() + DLB.sendTo[0], DLB.transferNmesh.begin() + DLB.sendTo[i], 0) * (18 + 4 * TwoDim.NS);
                        int size = DLB.transferNmesh[DLB.sendTo[i]] * (18 + 4 * TwoDim.NS);
                        if (size != 0)
                        {
                            // MPI_Recv(&DLB.recvUpdateFlowData[start], size, MPI_DOUBLE, std::get<1>(DLB.transferNchem[DLB.sendTo[i]]),
                            //         MPI_DOUBLE + 20000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                            MPI_Irecv(&DLB.recvUpdateFlowData[start], size, MPI_DOUBLE, std::get<1>(DLB.transferNchem[DLB.sendTo[i]]),
                                    DLB.sendTo[i] + 20000, MPI_COMM_WORLD, &recvRequests3[DLB.sendTo[i]]);
                        }
                    }
                    double tt4 = MPI_Wtime();
                    dt6 += tt4 - tt3;
                }

                // dt5 += tt2 - tt1;
                // dt6 += tt3 - tt2;
            }

            nchemSumBalance = 0;
            // #pragma omp parallel for num_thread(num_thread) collapse(3) schedule(dynamic, 100) reduction(+ : nchemSumBalance)
            for (int i = (ni - 1) / 2 + bc; i < ni + bc; i++)        // 本地剩余任务计算
                for (int j = bc; j < nj + bc; j++)
                    for (int k = bc; k < nk + bc; k++)
                    {
                        int n = 0;
                        bool flag = false;
                        while (n < DLB.sentNmesh && tag < DLB.sentNmesh)
                        {
                            if (i == DLB.transferIndex[3 * n] &&
                                j == DLB.transferIndex[3 * n + 1] &&
                                k == DLB.transferIndex[3 * n + 2])
                            {
                                flag = true;
                                tag ++;
                                break;
                            }
                            n++;
                        }
                        if(!flag)
                        {
                            nchemSumBalance += React.Nchem(i, j, k);
                            for (int step = 0; step < React.Nchem(i, j, k); step++)
                            {

                                React.Trapezoid(TwoDim.Mc, TwoDim.Di, TwoDim.Yi, TwoDim.T, TwoDim.dt, step, i, j, k);


                                TwoDim.Explicit(step, i, j, k);
                            }
                        }
                    }
            DLB.NchemSumBalance += nchemSumBalance;

            tw3 = MPI_Wtime();
            MPI_Waitall(DLB.transferNchem.size(), sendRequests3, statuses);
            MPI_Waitall(DLB.transferNchem.size(), recvRequests3, statuses);
            tw4 = MPI_Wtime();

            if (DLB.LoadDegree > DLB_tol)
            {
                int sentNchem = DLB.GetSentNchem(myid);
                if (sentNchem != 0)                                  // 将接收缓冲区 (recvUpdateFlowData) 里的数据，根据坐标解包回 TwoDim 主数组，迁出的那部分网格点的值被异地计算的结果覆盖更新
                    TwoDim.UnpackageUpdate(DLB.recvUpdateFlowData, DLB.transferIndex, DLB.sendTo, DLB.transferNmesh);
            }
            t5 = MPI_Wtime();

            // dt3 += t5 - t4;
            React.PushNchemMax();
            t6 = MPI_Wtime();
            // MPI_Barrier(MPI_COMM_WORLD);
            t7 = MPI_Wtime();
            tr3 += t3 - ttp2 + t5 - t4;
            tr4 += tw4 - tw3 + tw2 - tw1;


            if (iteration % DLB_step == 0)                           // 记录并输出均衡后的负载状态
            {
                MPI_Allgather(&DLB.NchemSumBalance, 1, MPI_INT, &React.NchemTotalBalance(0), 1, MPI_INT, MPI_COMM_WORLD);
                DLB.LoadDegreeBalance = double((React.NchemTotalBalance.MaxValue() - React.NchemTotalBalance.AveValue())) / React.NchemTotalBalance.MaxValue();
                if (myid == 0)
                {
                    // cout << "After DLB: ";
                    // React.NchemTotal.Print();
                    cout << "Balance Load Degree: " << DLB.LoadDegreeBalance << endl;
                    // cout << "Balanced Load: " ;
                    // React.NchemTotalBalance.Print();
                    // cout << "Diff: \n";
                    // for (int i = 0; i < React.NchemTotal.GetSize(); i++)
                    //     cout << React.NchemTotal(i) - React.NchemTotalBalance(i) << " ";
                    // cout << endl;
                    // std::cout << "\tTransfer records:" << std::endl;
                    // for (int i = 0; i < DLB.transferNchem.size(); i++)
                    // {
                    //     int sender = std::get<0>(DLB.transferNchem[i]);
                    //     int receiver = std::get<1>(DLB.transferNchem[i]);
                    //     int amount = std::get<1>(DLB.transferNchem[i]);
                    //     std::cout << "\tSender: " << sender << " Receiver: " << receiver << " Transfer Nchem: " << amount << " Transfer Nmesh: " << DLB.transferNmesh[i] << std::endl; 
                    // }
                }
            }
            int max_it = 0;
            MPI_Reduce(&React.NchemNow, &max_it, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
            if (myid == 0)
                cout << "Max Reaction iteration step: " << max_it << endl;
            cout << "Process " << myid << ": Reaction iteration step: " << React.NchemNow << endl; 
    }
}









// 利用 隐式处理反应源项-显式处理流动项 方法处理化学反应源项并更新流场物理量
void Euler::IMEX()
{
    double t1 = 0, t2 = 0, t3 = 0;
    
    t1 = MPI_Wtime();
    React.Diagonalized(TwoDim.Mc, TwoDim.Di, TwoDim.T, TwoDim.Partion_T);
    t2 = MPI_Wtime();
    tr1 += t2 - t1;

    TwoDim.Update_IMEX(React.CMS, React.MD);
    t3 = MPI_Wtime();
    tr2 += t3 - t2;
}












void Euler::DNN(int type)
{
    TwoDim.Update_after_Adv();

    

    TwoDim.Explicit();
}











void Euler::Output(int num)
{
    int ni = Mymesh.ni;
    int nj = Mymesh.nj;
    int nk = Mymesh.nk;
    int bc = Mymesh.bc;


    string rea = react;
    string filename;
    if (numprocs == 1)
        filename = "./output_" + to_string(numprocs) + "/results/Series" + rea + to_string(num) + +"_" + to_string(myid) + ".dat";
    else
        filename = "./output_" + to_string(numprocs) + "/results" + rea + to_string(num) + +"_" + to_string(myid) + ".dat";


        ofstream outfile;
        outfile.open(filename, ios::out);


        outfile << "Variables = X,Y,Z,T,D,P,U,V,Ma,Y1,Y2,Y3,Y4,Y5,Y6,Y7,L,WaitTime\n";
        outfile << "ZONE I=" << ni << '\t' << "J=" << nj << '\t' << "K=" << nk << '\n';
        outfile << "datapacking=block\n";
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << Mymesh.xnode(i) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << Mymesh.ynode(j) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << Mymesh.znode(k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.T(i, j, k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.D(i, j, k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.P(i, j, k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.U(i, j, k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.V(i, j, k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Ma(i, j, k) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 0) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 1) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 2) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 3) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 4) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 5) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << TwoDim.Yi(i, j, k, 6) << '\n';

        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << React.Nchem(i, j, k) << '\n';
        
        for (int k = bc; k < nk + bc; k++)
            for (int j = bc; j < nj + bc; j++)
                for (int i = bc; i < ni + bc; i++)
                    outfile << tr4 << '\n';


































































    outfile.close();

    return;
}










void Euler::Output_Total_Time(double t1, double t2, double t3, double t4, double t5, double t6, double t7, double t8, double t9, double t10)
{
    auto t12 = (t2 - t1);
    auto t23 = (t3 - t2);
    auto t34 = (t4 - t3);
    auto t45 = (t5 - t4);
    auto t56 = (t6 - t5);
    auto t67 = (t7 - t6);
    auto t78 = (t8 - t7);
    auto t89 = (t9 - t8);
    auto t90 = (t10 - t9);

    string filename = "./output_" + to_string(numprocs) + "/time_record/total/Total_Time_" + to_string(myid) + ".txt";
    ofstream outfile(filename);
    outfile << "myid:" << '\t' << myid << '\n';
    outfile << "FileReader:" << '\t' << t12 << '\n';
    outfile << "Mpi_Initial:" << '\t' << t23 << '\n';
    outfile << "MeshProcess:" << '\t' << t34 << '\n';
    outfile << "InputRead:" << '\t' << t45 << '\n';
    outfile << "ReactionRead:" << '\t' << t56 << '\n';
    outfile << "FieldInitial:" << '\t' << t67 << '\n';
    outfile << "MPI_Barrier:" << '\t' << t78 << '\n';
    outfile << "Computing:" << '\t' << t89 << '\n';
    outfile << "Output:" << '\t' << t90 << '\n';
    outfile.close();

    double maxt = 0;
    MPI_Reduce(&t89, &maxt, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (myid == 0)
    {
        string rea = react;
        filename = "./output_" + to_string(numprocs) + "/time_record/compare/" + rea + "_mpi_" + to_string(m_block_x * m_block_y) + "_omp_" + to_string(num_thread) + ".txt";
        ofstream outfile1(filename);
        outfile1 << maxt << '\t' << trackt;
        outfile1.close();
        std::cout << "Loop time: " << maxt << std::endl;
    }
}
void Euler::Output_Compute_Time(double t12, double t23, double t34, double t45, double t67)
{
    string file_compute = "./output_" + to_string(numprocs) + "/time_record/compute/Compute_Time_" + to_string(myid) + ".txt";
    ofstream outfile(file_compute);
    outfile << "myid:" << '\t' << myid << '\n';
    outfile << "FieldBoundary:" << '\t' << t12 << '\n';
    outfile << "Mpi_Boundary:" << '\t' << t23 << '\n';
    outfile << "Advection:" << '\t' << t34 << '\n';
    outfile << "Reaction:" << '\t' << t45 << '\n';
    outfile << "Output:" << '\t' << t67 << '\n';
    switch (Reaction_Sch)
    {
    case 0:
        outfile << "Update_after_Advection:" << '\t' << tr1 << '\n';
        outfile << "Trapezoid_Prediction:" << '\t' << dt4 << '\n';
        outfile << "Trapezoid:" << '\t' << tr2 + dt3 << '\n';
        outfile << "DLB:" << '\t' << tr3 - dt3 << '\n';
        outfile << "Wait:" << '\t' << tr4 << '\n';
        break;
    case 1:
        outfile << "Diagnalized:" << '\t' << tr1 << '\n';
        outfile << "Update_IMEX:" << '\t' << tr2 << '\n';
        outfile << "FileBoundary:" << TwoDim.ft1 << " " << TwoDim.ft2 << " " << TwoDim.ft3 << " " << TwoDim.ft4 << " " << TwoDim.ft5 << " " << TwoDim.ft5 << " " << TwoDim.ft6 << '\n';
        break;
    }

    outfile.close();
    double recv_wait_time[m_block_x * m_block_y];
    double recv_react_time[m_block_x * m_block_y];
    double recv_advect_time[m_block_x * m_block_y];
    double recv_DLB_time[m_block_x * m_block_y];
    double recv_DLB_time1[m_block_x * m_block_y];
    double recv_DLB_time2[m_block_x * m_block_y];
    double recv_DLB_time3[m_block_x * m_block_y];
    double recv_DLB_time4[m_block_x * m_block_y];
    double recv_DLB_time5[m_block_x * m_block_y];
    double recv_DLB_time6[m_block_x * m_block_y];
    double react_time = tr1 + tr2 + dt3 + dt4;
    double dlb_time = tr3 - dt3;
    MPI_Gather(&tr4, 1, MPI_DOUBLE, &recv_wait_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&react_time, 1, MPI_DOUBLE, &recv_react_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&t34, 1, MPI_DOUBLE, &recv_advect_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dlb_time, 1, MPI_DOUBLE, &recv_DLB_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dt1, 1, MPI_DOUBLE, &recv_DLB_time1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dt2, 1, MPI_DOUBLE, &recv_DLB_time2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dt3, 1, MPI_DOUBLE, &recv_DLB_time3, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dt5, 1, MPI_DOUBLE, &recv_DLB_time4, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dt6, 1, MPI_DOUBLE, &recv_DLB_time5, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&dt7, 1, MPI_DOUBLE, &recv_DLB_time6, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (myid == 0)
    {
        string file_total_time = "./output_" + to_string(numprocs) + "/time_record/all/All_wait.txt";
        ofstream outfile1(file_total_time);
        for (int i = 0; i < m_block_x * m_block_y; i++)
            outfile1 << recv_wait_time[i] << endl;
        outfile1.close();

        file_total_time = "./output_" + to_string(numprocs) + "/time_record/all/All_react.txt";
        ofstream outfile2(file_total_time);
        for (int i = 0; i < m_block_x * m_block_y; i++)
            outfile2 << recv_react_time[i] << endl;
        outfile2.close();

        file_total_time = "./output_" + to_string(numprocs) + "/time_record/all/All_advert.txt";
        ofstream outfile3(file_total_time);
        for (int i = 0; i < m_block_x * m_block_y; i++)
            outfile2 << recv_advect_time[i] << endl;
        outfile3.close();

        file_total_time = "./output_" + to_string(numprocs) + "/time_record/all/All_DLB.txt";
        ofstream outfile4(file_total_time);
        for (int i = 0; i < m_block_x * m_block_y; i++)
            outfile4 << recv_DLB_time[i] << '\t' << recv_DLB_time1[i] << '\t' << recv_DLB_time2[i] << '\t'
                     << recv_DLB_time3[i] << '\t' << recv_DLB_time4[i] << '\t' << recv_DLB_time5[i] << '\t' << recv_DLB_time6[i] << endl;
        outfile4.close();
    }
}
