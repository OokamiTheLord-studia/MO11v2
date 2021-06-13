#include "net.h"
#include <fstream>

namespace MO
{

	Net::Net(
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
	) :
		left_edge_condition_derivative(left_edge_condition_derivative)
		, left_edge_condition_function(left_edge_condition_function)
		, left_edge_condition_free_function(left_edge_condition_free_function)
		, right_edge_condition_derivative(right_edge_condition_derivative)
		, right_edge_condition_function(right_edge_condition_function)
		, right_edge_condition_free_function(right_edge_condition_free_function)
		, h(h)
	{


		const size_t x_count{
		static_cast<size_t>(std::floor((b - a) / h))
		};
		const size_t t_count{ static_cast<size_t>(std::floor((d - c) / dt)) };


		matrix.resize(t_count);
		for (auto& i : matrix)
		{
			i.resize(x_count);
		}
		x_positions.resize(x_count);
		t_positions.resize(t_count);

		auto& matrix_first_row = matrix.front();
		x_values.insert({ a, 0 });
		x_positions[0] = a;
		matrix_first_row.front() = start_condition(a);
		{
			double i{ a + h };
			unsigned int idx{ 1 };
			while (idx < x_count - 1)
			{
				x_values.insert({ i, idx });
				x_positions[idx] = i;
				matrix_first_row[idx] = start_condition(i);

				i += h;
				idx++;
			}
		}
		x_values.insert({ b, x_count - 1 });
		x_positions[x_count - 1] = b;
		matrix_first_row.back() = start_condition(b);

		t_values.insert({ c, 0 });
		t_positions[0] = c;
		{
			double i{ c + dt };
			unsigned int idx{ 1 };
			while (idx < t_count - 1)
			{
				t_values.insert({ i, idx });
				t_positions[idx] = i;
				for (auto j{ std::next(matrix[idx].begin()) }; j < std::prev(matrix[idx].end()); j++)
				{
					*j = 0;
				}


				i += dt;
				idx++;
			}
		}

		t_values.insert({ d, t_count - 1 });
		t_positions[t_count - 1] = d;
		auto& matrix_last_row = matrix.back();
		for (auto j{ std::next(matrix_last_row.begin()) }; j < std::prev(matrix_last_row.end()); j++)
		{
			*j = 0;
		}
	};



	double& Net::at(const double t, const double x)
	{
		return matrix.at(t_values.at(t)).at(x_values.at(x));
	}

	std::vector<std::vector<double>>* Net::getMatrix()
	{
		return &matrix;
	}

	std::vector<double>* Net::getTPositions()
	{
		return &t_positions;
	}

	void Net::dump(std::string filename)
	{
		std::ofstream file;
		file.open(filename);
		for (auto it = x_values.cbegin(); it != x_values.cend(); it++)
		{
			file << "," << it->first;
		}
		file << std::endl;
		for (auto t = t_values.crbegin(); t != t_values.crend(); t++)
		{
			file << t->first;
			for (auto x = x_values.cbegin(); x != x_values.cend(); x++)
			{
				file << "," << matrix.at(t->second).at(x->second);
			}
			file << std::endl;
		}
		file.close();
	}

	void Net::dumpForT(std::string filename, size_t t)
	{
		std::ofstream file;
		file.open(filename);
		for (auto ix_values{ x_values.cbegin() }; ix_values != x_values.cend(); ix_values++)
		{
			file << ix_values->first << "," << matrix.at(t).at(ix_values->second) << std::endl;
		}
		file.close();
	}

	void Net::dumpForT(std::string filename, double t)
	{
		dumpForT(filename, t_values.at(t));
	}

	void Net::solve(SolvingMethod* method)
	{
		method->solveNet(this);
	}

	void Net::solveLeftEdgeCondition(size_t t_pos)
	{

		auto u0 = matrix.at(t_pos).at(0);
		auto u1 = matrix.at(t_pos).at(1);

		u0 = (1 / (left_edge_condition_function(t_positions[t_pos]) - left_edge_condition_derivative(t_positions[t_pos]))) * (((-left_edge_condition_derivative(t_positions[t_pos]) / h) * u1) - left_edge_condition_free_function(t_positions[t_pos]));
	}

	void Net::solveLeftEdgeCondition(double t)
	{
		this->solveLeftEdgeCondition(t_values.at(t));
	}

	void Net::solveRightEdgeCondition(size_t t_pos)
	{

		auto u0 = matrix.at(t_pos).at(0);
		auto u1 = matrix.at(t_pos).at(1);

		u0 = (1 / (right_edge_condition_function(t_positions[t_pos]) - right_edge_condition_derivative(t_positions[t_pos]))) * (((-right_edge_condition_derivative(t_positions[t_pos]) / h) * u1) - right_edge_condition_free_function(t_positions[t_pos]));
	}

	void Net::solveRightEdgeCondition(double t)
	{
		this->solveRightEdgeCondition(t_values.at(t));
	}

	double Net::stepToTime(size_t t_pos)
	{
		return t_positions.at(t_pos);
	}

	double Net::positionToCoordinate(size_t x_pos)
	{
		return x_positions.at(x_pos);
	}
}

