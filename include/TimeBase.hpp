#ifndef TIMEBASE
#define TIMEBASE
#include <chrono>
using namespace std::chrono;

class TimeBase
{
private:
    system_clock::time_point m_start;

public:
    TimeBase() : m_start(system_clock::time_point::min()) { ; }


    void Clear()
    {
        m_start = system_clock::time_point::min();
    }


    bool IsStart() const
    {
        return (m_start.time_since_epoch() != system_clock::duration(0));
    }


    void Start()
    {
        m_start = system_clock::now();
    }


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