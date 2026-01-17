#include <LoadBalance.hpp>
#include <algorithm>
#include <numeric>
#include <Global.hpp>
#include <limits.h>
#include <tuple>
#include <unordered_map>
void DynamicLoadBalancer::DLBPriorityQueue(Array<int, 1> &arr)                                     // arr 数组存储每个进程当前的计算负荷，通常是该进程拥有的化学反应步数总和
{
    int n = arr.GetSize();
    int average = 0;
    int sum;

    average = arr.AveValue();

    std::priority_queue<std::pair<int, int>> maxHeap;                                              // 放入过载的进程
    std::priority_queue<std::pair<int, int>> minHeap;                                              // 最闲的进程


    for (int i = 0; i < n; i++)
    {
        if (arr(i) >= average)
        {
            maxHeap.push({arr(i), i});
        }
        else
        {
            minHeap.push({arr(i), i});
        }
    }


    while (!maxHeap.empty())
    {
        int sender = maxHeap.top().second;                                                         // 表示当前计算任务过重的进程编号
        int diff = maxHeap.top().first - average;                                                  // 当前 sender 比平均多多少任务
        maxHeap.pop();

        while (!minHeap.empty() && diff > 0)
        {
            int receiver = minHeap.top().second;                                                   // 表示当前计算任务过轻的进程编号
            int surplus = average - minHeap.top().first;                                           // 当前 receiver 比平均少多少任务，剩余空间
            minHeap.pop();

            int transferAmount = std::min(diff, surplus);
            arr(sender) -= transferAmount;
            arr(receiver) += transferAmount;

            transferNchem.push_back(std::make_tuple(sender, receiver, transferAmount));            // <sender, receiver, transferAmount>，迁移了多少计算负荷            diff -= transferAmount;

            diff -= transferAmount;

            if(arr(receiver) < average)
            {
                minHeap.push({arr(receiver), receiver});
            }
        }
    }
}
// 迁移当前进程需要负责发送的所有任务信息
int DynamicLoadBalancer::GetSentNchem(int element)                                                 // element 通常是 MPI Rank
{
    sendTo.clear();
    int totalSentAmount = 0;
    int line = 0;
    for (const auto &transfer : transferNchem)                                                     // 记录发送给每个接收者的网格点总数
    {
        int sender = std::get<0>(transfer);                                                        // 从 transferNchem 中获取信息 <sender, receiver, transferAmount>
        int receiver = std::get<1>(transfer);
        int amount = std::get<2>(transfer);

        if (sender == element)
        {
            sendTo.push_back(line);                                                                // 记录该数据在 transferNchem 中的位置
            totalSentAmount += amount;                                                             // 累计自己需要由其他进程帮忙承担的载荷
        }
        line++;
    }
    return totalSentAmount;
}

int DynamicLoadBalancer::GetRecvNmesh(int element)
{
    recvFrom.clear();
    int totalRecvAmount = 0;
    int line = 0;
    for (const auto &transfer : transferNchem)
    {
        int sender = std::get<0>(transfer);
        int receiver = std::get<1>(transfer);

        if(receiver == element)
        {
            recvFrom.push_back(line);
            totalRecvAmount += transferNmesh[line];
        }
        line++;
    }
    return totalRecvAmount;
}
// 从本进程的三维网格中，具体挑选出哪些网格点被打包带走
void DynamicLoadBalancer::GetTransferMesh(Array<int, 3> &Nchem, int bc)
{
    int ni = Nchem.GetNi() - 2 *bc;
    int nj = Nchem.GetNj() - 2 *bc;
    int nk = Nchem.GetNk() - 2 *bc;
    std::vector<int> nchem(Nchem.Getbuf(), Nchem.Getbuf() + Nchem.GetSize());
    std::vector<int> index(nchem.size());
    std::iota(index.begin(), index.end(), 0);

    // Global::quickSort(nchem, index, 0, nchem.size() - 1);

    std::vector<int> targetSums;
    for (int line : sendTo)
        targetSums.push_back(std::get<2>(transferNchem[line]));

    std::vector<std::vector<int>> packNchem;
    std::vector<std::vector<int>> packIndices;
    packNchem.resize(targetSums.size());
    packIndices.resize(targetSums.size());

    // for (int n = nchem.size() - 1; n >=0; n --)
    for(int n = 0; n < nchem.size(); n++)
    {
        std::tuple<int, int, int> arrayIndex = Nchem.Get3DIndices(index[n]);
        if (std::get<0>(arrayIndex) >= bc && std::get<1>(arrayIndex) >= bc && std::get<2>(arrayIndex) >= bc &&
            std::get<0>(arrayIndex) < ni + bc && std::get<1>(arrayIndex) < nj + bc && std::get<2>(arrayIndex) < nk + bc)
        {
            int num = nchem[n];
            int minDiff = INT_MAX;
            int targetIndex = -1;


            for (int i = 0; i < targetSums.size(); i++)
            {
                int currentSum = std::accumulate(packNchem[i].begin(), packNchem[i].end(), 0);     // 当前已经装载的量
                int currentDiff = std::abs(targetSums[i] - currentSum);

                if (currentDiff < minDiff && currentSum + num <= targetSums[i])                    // 如果还装得下
                {
                    minDiff = currentDiff;
                    targetIndex = i;
                }
            }

            if (targetIndex != -1)
            {
                packNchem[targetIndex].push_back(num);
                packIndices[targetIndex].push_back(std::get<0>(arrayIndex));                       // 记录 X 坐标
                packIndices[targetIndex].push_back(std::get<1>(arrayIndex));                       // 记录 Y 坐标
                packIndices[targetIndex].push_back(std::get<2>(arrayIndex));                       // 记录 Z 坐标
                transferMeshNum[sendTo[targetIndex]]++;                                            // 交易的网格数加 1
            }
        }
    }

    for (const auto &row : packIndices)
        transferIndex.insert(transferIndex.end(), row.begin(), row.end());

        // for (int i = 0; i < packNchem.size(); i++)
        // {
        //     std::cout << "Array " << i + 1 << ": ";
        //     for (int num : packNchem[i])
        //     {
        //         std::cout << num << " ";
        //     }
        //     std::cout << std::endl;

        //     std::cout << "Indices " << i + 1 << ": ";
        //     for (int index : packIndices[i])
        //     {
        //         std::cout << index << " ";
        //     }
        //     std::cout << std::endl;
        // }
}