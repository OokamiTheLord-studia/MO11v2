#pragma once

#include <cstddef>
#include <array>
#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <cmath>
#include <iostream>
#include "simpleLogger.h"



namespace MO
{
	//TODO: Przenieœæ definicje do cpp
	//TODO: Napisaæ funkcje do dostêpu
	class Net
	{
		std::vector<double> x_values;
		std::vector<double> t_values;
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
		)
		{
			//TODO: Przemyœleæ optymalizacjê
			//TODO: Asercja wartoœci

			LOG("Hello, Net constructor here!")

			

			const unsigned int x_count{ static_cast<unsigned int>(std::floor((b - a) / h)) };
			const unsigned int t_count{ static_cast<unsigned int>(std::floor((d - c) / dt)) };

			LOG("x_count - " << x_count)
			LOG("t_count - " << t_count)
			//reserve
			x_values.resize(x_count);
			t_values.resize(t_count);
			matrix.resize(t_count);
			for (auto& i : matrix)
			{
				i.resize(x_count);
			}

			LOG("Net reserved memory")

			//fill x
			auto& matrix_first_row = matrix.front();
			x_values.front() = a;
			matrix_first_row.front() = start_condition(a);
			//mo¿e da siê tu u¿yæ iteratora?
			{
				double i{ a + h };
				unsigned int idx{ 1 };
				while (idx < x_values.size())
				{
					x_values[idx] = i;
					matrix_first_row[idx] = start_condition(i);

					i += h; 
					idx++;
				}
			}
			x_values.back() = b;
			matrix_first_row.back() = start_condition(b);

			LOG("Net filled first row")
			//fill edges and set middle to 0
			t_values.front() = c;
			{
				double i{ c + dt };
				unsigned int idx{ 1 };
				while (idx < t_values.size())
				{
					t_values[idx] = i;
					matrix[idx].front() = left_edge_condition(i);
					for (auto j{ std::next(matrix[idx].begin()) }; j < std::prev(matrix[idx].end()); j++)
					{
						*j = 0;
					}
					matrix[idx].back() = right_edge_condition(i);

					if (!(idx % 100))
					{
						LOG("Net filled row " << idx)
					}

					i += dt;
					idx++;
				}
			}

			LOG("Net will fill last row now")

			t_values.back() = d;
			auto& matrix_last_row = matrix.back();
			matrix_last_row.front() = left_edge_condition(d);
			for (auto j{ std::next(matrix_last_row.begin()) }; j < std::prev(matrix_last_row.end()); j++)
			{
				*j = 0;
			}
			matrix_last_row.back() = right_edge_condition(d);

			LOG("Constructor is done. Thank you forever")
		};
	};
};

