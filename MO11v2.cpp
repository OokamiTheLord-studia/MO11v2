// MO11v2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <functional>;
#include"net.h"

int main()
{
    constexpr double t_max{ 1 };
    constexpr double b{ 0.1 };
    constexpr double d{ 1 };

    
    double x_end{6 * sqrt(d*t_max)};
    double x_begin{-x_end};
    constexpr double t_begin{ 0 };
    constexpr double h{ 0.01 };
    constexpr double dt{ 0.001 };

    std::function<double(double)> temporary_function{ [](double temp) {return 0; } };

    MO::Net tempNet(x_begin, x_end, h, t_begin, t_max, dt, temporary_function, temporary_function, temporary_function);
}
