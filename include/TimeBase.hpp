#ifndef TIMEBASE
#define TIMEBASE
#include <chrono>
using namespace std::chrono;

class TimeBase
{
private:
    system_clock::time_point m_start;

public:
    TimeBase() : m_start(system_clock::time_point::min()) { ; }      // 构造函数用于初始化

    // 重置计时器：将开始时间恢复到最小值
    void Clear()
    {
        m_start = system_clock::time_point::min();
    }

    // 判断计时器是否已经启动
    bool IsStart() const
    {
        return (m_start.time_since_epoch() != system_clock::duration(0));
    }

    // 按下秒表
    void Start()
    {
        m_start = system_clock::now();
    }

    // 计算耗时
    unsigned long GetMs()
    {
        if(IsStart())
        {
            system_clock::duration diff;
            diff = system_clock::now() - m_start;
            return (unsigned)(duration_cast<milliseconds>(diff).count());
        }
        return 0;
    }
};
#endif