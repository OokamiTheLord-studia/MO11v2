#pragma once
#include <cstddef>
#include <array>
#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <cmath>
#include <iostream>
#include <string>
#include <map>

namespace MO
{
    class FillableNet
    {

		std::map<double, unsigned int> x_values;
		std::map<double, unsigned int> t_values;
		std::vector<std::vector<double>> matrix;

	public:
		FillableNet(
			const double a
			, const double b
			, const double h
			, const double c
			, const double d
			, const double dt
			, std::function<double(double, double)> filling_function
		);
		
		
		double& at(const double t, const double x);

		void dump(std::string filename);
		void dumpForT(std::string, size_t);
		void dumpForT(std::string, double);

		std::vector<std::vector<double>>* getMatrix();

    };

}
