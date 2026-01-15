





#ifndef TIME
#define TIME

#include "Array.hpp"
using namespace ARRAY;

class TimeAdv
{
public:

    TimeAdv() = default;


    void EE(int, double, Array<double,1>&, Array<double,1>&, Array<double,1>&, int, Array<double,4>&, Array<double,4>&, Array<double,4>&, Array<double,4>&);


    void TVD_RK3(int,double,Array<double,1>&,Array<double,1>&,Array<double,1>&, int, Array<double,4>&,Array<double,4>&,Array<double,4>&,Array<double,4>&);


    ~TimeAdv(){ ; };
};

#endif
