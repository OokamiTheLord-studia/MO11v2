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
	//TODO: Przenie�� definicje do cpp
	//TODO: Napisa� funkcje do dost�pu
	//	at przyjmuj�cy integery
	//  getx zwracaj�cy warto��
	//  gett zwracaj�cy warto��
	class Net
	{
		//std::vector<double> x_values;
		//std::vector<double> t_values;
		//TODO: Zastanowi� si� czy korzysta� z mapy, wektora, czy dw�ch map w obie strony (mapowanie czego na co b�dzie cz�ciej potrzebne, i dost�p po czym)
	protected:
		std::map<double, unsigned int> x_values;
		std::vector<double> x_positions;
		std::map<double, unsigned int> t_values;
		std::vector<double> t_positions;
		std::vector<std::vector<double>> matrix;
		std::function<double(double)> left_edge_condition_second_derivative;
		std::function<double(double)> left_edge_condition_first_derivative;
		std::function<double(double)> left_edge_condition_function;
		std::function<double(double)> right_edge_condition_second_derivative;
		std::function<double(double)> right_edge_condition_first_derivative;
		std::function<double(double)> right_edge_condition_function;

	public:
		Net(
			const double a
			, const double b
			, const double h
			, const double c
			, const double d
			, const double dt
			//sprawdzi� czy mog� by� const
			, std::function<double(double)> start_condition
			/*, std::function<double(double)> left_edge_condition
			, std::function<double(double)> right_edge_condition*/
			, std::function<double(double)> left_edge_condition_second_derivative
			, std::function<double(double)> left_edge_condition_first_derivative
			, std::function<double(double)> left_edge_condition_function
			, std::function<double(double)> right_edge_condition_second_derivative
			, std::function<double(double)> right_edge_condition_first_derivative
			, std::function<double(double)> right_edge_condition_function
		);

		

		double& at(const double t, const double x);

		std::vector<std::vector<double>>* getMatrix();

		void dump(std::string filename);

		void solve(MO::SolvingMethod&);
	};
};

