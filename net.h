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
	//TODO: Przenieœæ definicje do cpp
	//TODO: Napisaæ funkcje do dostêpu
	//	at przyjmuj¹cy integery
	//  getx zwracaj¹cy wartoœæ
	//  gett zwracaj¹cy wartoœæ
	class Net
	{
		//std::vector<double> x_values;
		//std::vector<double> t_values;
		//TODO: Zastanowiæ siê czy korzystaæ z mapy, wektora, czy dwóch map w obie strony (mapowanie czego na co bêdzie czêœciej potrzebne, i dostêp po czym)
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
			//sprawdziæ czy mog¹ byæ const
			, std::function<double(double)> start_condition
			/*, std::function<double(double)> left_edge_condition
			, std::function<double(double)> right_edge_condition*/
			, std::function<double(double)> left_edge_condition_derivative
			, std::function<double(double)> left_edge_condition_function
			, std::function<double(double)> left_edge_condition_free_function
			, std::function<double(double)> right_edge_condition_derivative
			, std::function<double(double)> right_edge_condition_function
			, std::function<double(double)> right_edge_condition_free_function
		);

		

		double& at(const double t, const double x);

		std::vector<std::vector<double>>* getMatrix();

		void dump(std::string filename);

		void solve(SolvingMethod*);

		void solveLeftEdgeCondition(size_t);
		void solveLeftEdgeCondition(double);

		void solveRightEdgeCondition(size_t);
		void solveRightEdgeCondition(double);

		double stepToTime(size_t);
		double positionToCoordinate(size_t);
	};
};

