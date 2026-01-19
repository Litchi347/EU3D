







#include <iostream>





#include "Global.hpp"














// 判断输入字符行的属性，决定它是有效内容、空白行还是注释行
int Global::SenJud(char *sentence)
{
    int i = 0;
    char c = sentence[i];
    if (c == '#')
        return 2;
    if (c == '\n' || c == '\0' || c == ' ')
        return 1;
    else
        return 0;
}












// 提取一个字符行 sentence 的键值对信息，键存入 parameter，值存入 value
void Global::ParaGet(char *sentence, char *parameter, char *value)
{
    int i = 0;
    char c = sentence[i];
    // 提取参数名到 parameter 中，并在结尾处增加一个\0符号
    while(c != ' ')
    {
        parameter[i] = sentence[i];
        i++;
        c = sentence[i];
    }
    parameter[i] = '\0';
    // 跳过空格
    while(c == ' ')
    {
        i++;
        c = sentence[i];
    }
    // 提取数值到 value ，并在结尾加上 \0 符号
    int j = 0;
    while(c != '\0')
    {
        value[j] = sentence[i];
        i++;
        j++;
        c = sentence[i];
    }
    value[j] = '\0';
}

int myid;
int numprocs;
int m_block_x, m_block_y, m_block_z;
int m_left = 0, m_right = 0, m_up = 0, m_down = 0, m_front = 0, m_back = 0;
int myid_x, myid_y, myid_z;
int num_thread;