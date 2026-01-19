



#pragma once

#include "Array.hpp"
#include <vector>
using namespace ARRAY;
#ifndef GLOBAL
#define GLOBAL

namespace Global
{
    // 结构体
    typedef enum
    {
        Trapezoid,
        IMEX,
        DNN
    } REACT;

    // 结构体
    typedef enum
    {
        MUSCL_1,
        MUSCL_2,
        WENO
    } DIFF;

    // 结构体
    typedef enum
    {
        EE,
        TVD_RK3
    } ADV;

    // 判断输入字符行是有效内容、空白行还是注释行
    int SenJud(char *sentence);

    // 得到一个形如 "value 10" 字符行的键值对信息
    void ParaGet(char *sentence, char *parameter, char *value);


    template <typename T>
    void swap(T &a, T &b)
    {
        T temp = a;
        a = b;
        b = temp;
    }
    
    template<typename T>                                                            // 快速排序中的分区操作：所有大于基准值的数在基准值左边，所有小于基准值的数在基准值右边
    int partition(std::vector<T> &arr, std::vector<T> &index, int low, int high)    // index 通常存储原始索引或相关联数据，跟随 arr 变动而变动，从而得知排序后的数据原来在哪
    {
        T pivot = arr[high];                                                        // 基准值
        int i = low - 1;

        for (int j = low; j < high; j ++)
        {
            if (arr[j] > pivot)
            {
                i++;
                swap(arr[i], arr[j]);
                swap(index[i], index[j]);
            }
        }

        swap(arr[i + 1], arr[high]);
        swap(index[i+1], index[high]);
        return i + 1;
    }

    template <typename T>
    void quickSort(std::vector<T> &arr, std::vector<T> &index, int low, int high)   // 快排
    {
        if (low < high)
        {

            int pi = partition(arr, index, low, high);


            quickSort(arr, index, low, pi - 1);                                     // 迭代
            quickSort(arr, index, pi + 1, high);
        }
    }

    template <typename T>
    int findIndex(std::vector<T> &index, int target)                                // 在 index 中寻找值为 target 的项的索引
    {
        for (int i = 0; i < index.GetSize(); i++)
        {
            if (index[i] == target)
            {
                return i;
            }
        }

        return -1;
    }


    template <typename T>
    T Power(T a, size_t b)                                                          // 计算 a 的 b 次方的值
    {
        T ans = 1;
        for (size_t i = 0; i < b; i++)
            ans *= a;
        return ans;
    }
}
#endif

extern int myid;                                                                    // 当前进程编号
extern int numprocs;                                                                // 总进程数
extern int m_block_x, m_block_y, m_block_z;                                         // 分配给当前进程处理的局部网格规模
extern int m_left, m_right, m_up, m_down, m_front, m_back;                          // X 方向上的左邻居和右邻居的 myid，Y 方向上的上邻居和下邻居的 myid，Z 方向上的前邻居和后邻居的 myid
extern int myid_x, myid_y, myid_z;                                                  // 当前进程在 X、Y、Z 三个方向上的坐标索引
extern int num_thread;                                                              // OpenMP 线程数