// MO11v2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <cmath>
#include <functional>
#include "net.h"
#include "FillableNet.h"
#include "KMB.h"
#include "CN_Thomas.h"
#include "CN_GaussSiedel.h"

constexpr double t_max{ 1 };
constexpr double b{ 0.1 };
constexpr double d{ 1 };

double analyticSolution(double t, double x)
{
    return 0.5 * exp(((d * t) / (b * b)) - (x / b)) * erfc((((2 * d * t) /b) - x) / (2 * sqrt(d * t)));
}

int main()
{
    double x_end{6 * sqrt(d*t_max)};
    double x_begin{-x_end};
    constexpr double t_begin{ 0 };
    constexpr double h{ 0.01 };
    constexpr double dt{ 0.001 };

    std::function<double(double)> temporary_function{ [](double temp) {return 0; } };

    std::function<double(double)> start_condition{ [&](double x) {return x < 0 ? 0 : exp(-x / b); } };
    std::function<double(double)> corrected_start_condition{ [&](double x) {return x > 0 ? 0 : exp(-x / b); } };
    std::function<double(double)> edge_condition_derivative_parameter{ [](double t) {return 0; } };
    std::function<double(double)> edge_condition_function_parameter{ [](double t) {return 1; } };
    std::function<double(double)> edge_condition_free_function_parameter{ [](double t) {return 0; } };



    //MO::Net tempNet(x_begin, x_end, h, t_begin, t_max, dt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    /*MO::Net tempNet(x_begin, x_end, h, t_begin, t_max, dt, temporary_function, temporary_function, temporary_function);
    std::cout << tempNet.at(t_begin, x_begin) << std::endl;
    std::cout << tempNet.at(t_max, x_end) << std::endl;
    tempNet.at(t_max, x_end) = 0.75;
    std::cout << tempNet.at(t_max, x_end) << std::endl;*/

    //tempNet.dump("testfile2.csv");

    /*MO::FillableNet analyticNet(x_begin, x_end, h, t_begin, t_max, dt, analyticSolution);
    analyticNet.dump("analyticSolution.csv");*/

    constexpr double KMBh{ 0.05 };
    constexpr double KMBlambda{ dt / (KMBh * KMBh) };

    MO::Net KMBSolvedNet(x_begin, x_end, KMBh, t_begin, t_max, dt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    MO::KMB kmbSolver(KMBlambda);
    /*KMBSolvedNet.solve(&kmbSolver);
    KMBSolvedNet.dump("KMB_corrected.csv");*/

    constexpr double CNh{ 0.05 };
    constexpr double CNdt{ 0.0025 };

    MO::Net CNThomasSolvedNet(x_begin, x_end, CNh, t_begin, t_max, CNdt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    MO::Crank_Nicolson::CN_Thomas CNThomasSolver{ CNdt / (CNh * CNh) };
    /*CNThomasSolvedNet.solve(&CNThomasSolver);
    CNThomasSolvedNet.dump("CNThomas.csv");*/

    MO::Net CNGaussSiedelSolvedNet(x_begin, x_end, CNh, t_begin, t_max, CNdt, start_condition, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter, edge_condition_derivative_parameter, edge_condition_function_parameter, edge_condition_free_function_parameter);
    MO::Crank_Nicolson::CN_GaussSiedel CNGaussSiedelSolver{ CNdt / (CNh * CNh) };
    CNGaussSiedelSolvedNet.solve(&CNGaussSiedelSolver);
    CNGaussSiedelSolvedNet.dump("CNGaussSiedel.csv");
}
