#include "CN_Thomas.h"
//#include <iterator>
#include <limits>
#include <iostream>
#include "simpleLogger.h"

namespace MO
{
	namespace Crank_Nicolson
	{
		using namespace MO;

		void CN_Thomas::solveNet(Net* net)
		{
			auto matrix{ net->getMatrix() };

			const double half_lambda{ lambda / 2 };
			const double minus_one_lambda{ -1 - lambda };
			const double lambda_minus_one{ 1 - lambda };


			/*for (
				auto past_iterator{ matrix->begin() }, current_iterator{ std::next(past_iterator) };
				current_iterator != matrix->end();
				past_iterator++, current_iterator++
				)*/
			for (
				unsigned int past_position = 0, current_position = 1;
				current_position < matrix->size();
				past_position++, current_position++
				)
			{
				const auto current_time{ net->stepToTime(current_position) };

				// tworzenie macierzy A i wektora b
				const auto x_size{ matrix->at(0).size() };
				std::vector<double> l;
				std::vector<double> d;
				std::vector<double> u;
				std::vector<double> b;

				l.reserve(x_size - 1);
				d.reserve(x_size);
				u.reserve(x_size - 1);
				b.reserve(x_size);

				u.push_back(net->left_edge_condition_derivative(current_time) / net->h);
				d.push_back(-(net->left_edge_condition_derivative(current_time)) / net->h + net->left_edge_condition_function(current_time));
				b.push_back(-(net->left_edge_condition_free_function(current_time)));
				for (size_t i = 1; i <= x_size - 2; i++)
				{
					l.push_back(half_lambda);
					u.push_back(half_lambda);
					d.push_back(minus_one_lambda);
					b.push_back(-((half_lambda * matrix->at(past_position).at(i - 1)) + (lambda_minus_one * matrix->at(past_position).at(i)) + (half_lambda * matrix->at(past_position).at(i + 1))));
				}
				l.push_back(-(net->right_edge_condition_derivative(current_time)) / net->h);
				b.push_back(-(net->right_edge_condition_free_function(current_time)));
				d.push_back(net->right_edge_condition_derivative(current_time) / net->h + net->right_edge_condition_function(current_time));

				//uzyskanie rozwi¹zañ
				auto solutions = solveLinearEquation(u, d, l, b);

				matrix->at(current_position).swap(solutions);

			}
		}

		std::vector<double> CN_Thomas::solveLinearEquation(std::vector<double>& u, std::vector<double>& d, std::vector<double>& l, std::vector<double>& b)
		{
			const size_t element_count{ d.size() };

			std::vector<double> mi;
			std::vector<double> ni;
			mi.reserve(element_count);
			ni.reserve(element_count);

			mi.push_back(d[0]);
			ni.push_back(b[0]);
			for (size_t i = 1; i < element_count; i++)
			{
				mi.push_back(d[i] - l[i - 1] * (1 / mi[i - 1]) * u[i - 1]);
				ni.push_back(b[i] - l[i - 1] * (1 / mi[i - 1]) * ni[i - 1]);
			}

			std::vector<double> solutions;
			solutions.resize(element_count);

			solutions[element_count - 1] = (1 / mi[element_count - 1]) * ni[element_count - 1];
			for (size_t i = element_count - 2; i != std::numeric_limits<size_t>::max(); i--)
			{
				solutions[i] = (1 / mi[i]) * (ni[i] - u[i] * solutions[i + 1]);
			}

			return solutions;
		}

	}
}