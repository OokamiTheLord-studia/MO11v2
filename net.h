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
#include "SolvingMethod.h"



namespace MO
{
	class SolvingMethod;
	class Net
	{
	protected:
		std::map<double, size_t> x_values;
		std::vector<double> x_positions;
		std::map<double, size_t> t_values;
		std::vector<double> t_positions;
		std::vector<std::vector<double>> matrix;


	public:

		const std::function<double(double)> left_edge_condition_derivative;
		const std::function<double(double)> left_edge_condition_function;
		const std::function<double(double)> left_edge_condition_free_function;
		const std::function<double(double)> right_edge_condition_derivative;
		const std::function<double(double)> right_edge_condition_function;
		const std::function<double(double)> right_edge_condition_free_function;
		const double h;

		Net(
			const double a
			, const double b
			, const double h
			, const double c
			, const double d
			, const double dt
			, std::function<double(double)> start_condition
			, std::function<double(double)> left_edge_condition_derivative
			, std::function<double(double)> left_edge_condition_function
			, std::function<double(double)> left_edge_condition_free_function
			, std::function<double(double)> right_edge_condition_derivative
			, std::function<double(double)> right_edge_condition_function
			, std::function<double(double)> right_edge_condition_free_function
		);



		double& at(const double t, const double x);

		std::vector<std::vector<double>>* getMatrix();
		std::vector<double>* getTPositions();

		void dump(std::string filename);
		void dumpForT(std::string, size_t);
		void dumpForT(std::string, double);

		void solve(SolvingMethod*);

		void solveLeftEdgeCondition(size_t);
		void solveLeftEdgeCondition(double);

		void solveRightEdgeCondition(size_t);
		void solveRightEdgeCondition(double);

		double stepToTime(size_t);
		double positionToCoordinate(size_t);
	};
};

