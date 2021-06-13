#pragma once

void Zad1(double x_begin, double x_end, const double& t_begin, std::function<double(double)>& start_condition, std::function<double(double)>& edge_condition_derivative_parameter, std::function<double(double)>& edge_condition_function_parameter, std::function<double(double)>& edge_condition_free_function_parameter);
