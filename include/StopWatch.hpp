#ifndef STOPWATCH
#define STOPWATCH
# include "TimeBase.hpp"
# include <iostream>
# include "Global.hpp"
template <typename T>
class basic_stopwatch;
typedef basic_stopwatch<TimeBase> Stopwatch;

template <typename T>
class basic_stopwatch
{
    T BaseTimer;

public:

    explicit basic_stopwatch(bool start);
    explicit basic_stopwatch(char const *activity = "Stopwatch", bool start = true);



    ~basic_stopwatch();


    unsigned LapGet() const;


    bool IsStarted() const;


    void Show(char const *event_name = "show");


    void Start(char const *event_name = "start");


    void Stop(char const *event_name = "stop");

private:
};

template <typename T>
basic_stopwatch<T>::basic_stopwatch(bool start)                           // 构造函数1：通过布尔值决定是否立即开始计时
{
    if(start)
        Start("start");
}

template <typename T>
basic_stopwatch<T>::basic_stopwatch(char const *activity, bool start)     // 构造函数2：传入活动名称（如"计算化学项"），并默认启动计时
{
    std::cout << "Activity:\t" << activity << std::endl;
    if(start)
        Start();
}

template <typename T>
basic_stopwatch<T>::~basic_stopwatch()                                    // 自动结束
{
    Stop();
}

template <typename T>
unsigned basic_stopwatch<T>::LapGet() const                               // 获取从开始到现在经过的时间（毫秒），不会停止计时
{
    return BaseTimer.GetMs();
}

template <typename T>
bool basic_stopwatch<T>::IsStarted() const                                // 检查计时器当前是否正在运行
{
    return BaseTimer.IsStart()();
}

template <typename T>
void basic_stopwatch<T>::Show(char const *event_name)                     // 在不停止计时的情况下，打印当前已经过去的时间
{
    std::cout << "Activity:\t" << event_name << "\t" << BaseTimer.getMs() << " mS" << std::endl;
}

template <typename T>
void basic_stopwatch<T>::Start(char const *event_name)                    // 重置并开始计时，同时打印进程 ID ，方便在 MPI 并行环境下调试
{
    std::cout << "Process " << myid << ":\t" << event_name << std::endl;
    BaseTimer.Start();
}

template <typename T>
void basic_stopwatch<T>::Stop(char const *event_name)                     // 打印最终耗时，并清除计时器状态
{
    std::cout << "process " << myid << ":\t" << event_name << "\t" << BaseTimer.GetMs() << " mS" << std::endl;
    BaseTimer.Clear();
}
#endif

