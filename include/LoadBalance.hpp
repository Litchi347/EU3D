#ifndef DLB_H
#define DLB_H
# include <vector>
# include <queue>
# include <tuple>
# include <mpi.h>
# include "Array.hpp"

using namespace ARRAY;
class DynamicLoadBalancer
{
public:

    std::vector<std::tuple<int,int,int>> transferNchem;    // 记录谁给谁发送了多少载荷
    double LoadDegree = 0.0;                               // 负载不均衡度 =  (最大负载 - 平均负载) / 最大负载
    double LoadDegreeBalance = 0.0;                        // 平衡后的负载度
    int NchemSumBalance = 0;                               // 平衡后的全场化学反应总数
    std::vector<int> NchemTotalBalance;                    // 平衡后每一个进程的分得的总载荷量
    std::vector<int> sendTo;                               // 存储当前进程作为 sender 时，在 tranferNchem 中对应的行号
    std::vector<int> recvFrom;                             // 存储当前进程作为 receiver 时，在 tranferNchem 中对应的行号
    std::vector<int> transferIndex;                        // 存储被选中的网格点在原流场中的三维坐标
    std::vector<double> prevFlowDataLocal;
    std::vector<int> prevReactDataLocal;
    std::vector<double> updateFlowDataLocal;
    int sentNmesh = 0;                                     // 当前进程总共要送出的网格数量
    std::vector<int> transferMeshNum;                      // 每一笔迁移交易中包含的具体网格点数
    std::vector<int> transferNmesh;                        // 接收到的每一笔交易中包含的网格点数
    std::vector<double> transferPrevDataLocal;             // 存放待发送的上一步流场数据
    std::vector<double> transferPrevData;
    std::vector<double> transferUpdateDataLocal;           // 存放待发送的当前步流场数据
    std::vector<int> transferPrevNchem;
    std::vector<int> transferPrevNchemLocal;               // 存放待发送网格点的化学载荷量
    std::vector<double> TransferUpdateData;
    std::vector<int> TransferUpdateNchem;
    std::vector<int> TransferUpdateNchemLocal;
    std::vector<double> recvUpdateFlowData;
    std::vector<double> recvPrevFlowData;
    std::vector<int> recvPrevReactData;


    DynamicLoadBalancer() = default;
    void DLBPriorityQueue(Array<int,1> &arr);
    int GetSentNchem(int element);
    int GetRecvNmesh(int element);
    void GetTransferMesh(Array<int,3> &Nchem,int bc);
    ~DynamicLoadBalancer(){ ; };
};
#endif