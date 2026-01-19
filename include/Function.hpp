



#ifndef FUNCTION_H
#define FUNCTION_H





# include <iostream>
# include <fstream>
# include <math.h>





# include "Array.hpp"
using namespace ARRAY;
class Function
{
public:

    Function() = default;


    double sum(int flag, Array<double, 1> &a, Array<double, 1> &b);


    ~Function(){ ; };
};
#endif