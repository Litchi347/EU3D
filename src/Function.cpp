#include <iostream>
#include <fstream>
#include <math.h>

#include "Function.hpp"











double Function::sum(int flag, Array<double, 1> &a, Array<double, 1> &b)
{
    int size = a.GetNi();
    double sum = 0;
    switch(flag)
    {
    case 0:
        for(int i = 0; i < size; i++)
            sum += (a(i) + b(i));
        break;
    case 1:
        for(int i = 0; i < size; i++)
            sum += (a(i) - b(i));
        break;
    case 2:
        for(int i = 0; i < size; i++)
            sum += (a(i) * b(i));
        break;
    case 3:
        for(int i = 0; i < size; i++)
            sum += (a(i) / b(i));
        break;
    }

    return sum;
}