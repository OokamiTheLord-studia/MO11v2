#pragma once

#include <cstddef>
#include <array>
#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <cmath>
#include <iostream>
#include <map>



namespace MO
{
	//TODO: Przenieœæ definicje do cpp
	//TODO: Napisaæ funkcje do dostêpu
	class Net
	{
		//std::vector<double> x_values;
		//std::vector<double> t_values;
		std::map<double, unsigned int> x_values;
		std::map<double, unsigned int> t_values;
		std::vector<std::vector<double>> matrix;

	public:
		Net(
			const double a
			, const double b
			, const double h
			, const double c
			, const double d
			, const double dt
			//sprawdziæ czy mog¹ byæ const
			, std::function<double(double)> start_condition
			, std::function<double(double)> left_edge_condition
			, std::function<double(double)> right_edge_condition
		);

		double& at(const double t, const double x);
	};
};

