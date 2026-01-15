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

    std::vector<std::tuple<int,int,int>> transferNchem;
    double LoadDegree = 0.0;
    double LoadDegreeBalance = 0.0;
    int NchemSumBalance = 0;
    std::vector<int> NchemTotalBalance;
    std::vector<int> sendTo;
    std::vector<int> recvFrom;
    std::vector<int> transferIndex;
    std::vector<double> prevFlowDataLocal;
    std::vector<int> prevReactDataLocal;
    std::vector<double> updateFlowDataLocal;
    int sentNmesh = 0;
    std::vector<int> transferMeshNum;
    std::vector<int> transferNmesh;
    std::vector<double> transferPrevDataLocal;
    std::vector<double> transferPrevData;
    std::vector<double> transferUpdateDataLocal;
    std::vector<int> transferPrevNchem;
    std::vector<int> transferPrevNchemLocal;
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