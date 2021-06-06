#include "net.h"
#include "simpleLogger.h"
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
		//sprawdziæ czy mog¹ byæ const
		, std::function<double(double)> start_condition
		//, std::function<double(double)> left_edge_condition
		//, std::function<double(double)> right_edge_condition
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
		//TODO: Przemyœleæ optymalizacjê
		//TODO: Asercja wartoœci

		LOG("Hello, Net constructor here!")



			const unsigned int x_count{ static_cast<unsigned int>(std::floor((b - a) / h)) };
		const unsigned int t_count{ static_cast<unsigned int>(std::floor((d - c) / dt)) };

		LOG("x_count - " << x_count)
			LOG("t_count - " << t_count)
			//reserve
			//x_values.resize(x_count);
		//t_values.resize(t_count);
		matrix.resize(t_count);
		for (auto& i : matrix)
		{
			i.resize(x_count);
		}
		x_positions.resize(x_count);
		t_positions.resize(t_count);

		LOG("Net reserved memory")

			//fill x
			auto& matrix_first_row = matrix.front();
		//x_values.front() = a;
		x_values.insert({ a, 0 });
		x_positions[0] = a;
		matrix_first_row.front() = start_condition(a);
		//mo¿e da siê tu u¿yæ iteratora?
		{
			double i{ a + h };
			unsigned int idx{ 1 };
			//while (idx < x_values.size())
			while (idx < x_count - 1)
			{
				//x_values[idx] = i;
				x_values.insert({ i, idx });
				x_positions[idx] = i;
				matrix_first_row[idx] = start_condition(i);

				i += h;
				idx++;
			}
		}
		//x_values.back() = b;
		x_values.insert({ b, x_count - 1 });
		x_positions[x_count - 1] = b;
		matrix_first_row.back() = start_condition(b);

		LOG("Net filled first row")
			//set middle to 0
			//t_values.front() = c;
			t_values.insert({ c, 0 });
		t_positions[0] = c;
		{
			double i{ c + dt };
			unsigned int idx{ 1 };
			//while (idx < t_values.size())
			while (idx < t_count - 1)
			{
				//t_values[idx] = i;
				t_values.insert({ i, idx });
				t_positions[idx] = i;
				//matrix[idx].front() = left_edge_condition(i);
				for (auto j{ std::next(matrix[idx].begin()) }; j < std::prev(matrix[idx].end()); j++)
				{
					*j = 0;
				}
				//matrix[idx].back() = right_edge_condition(i);

				if (!(idx % 100))
				{
					LOG("Net filled row " << idx)
				}

				i += dt;
				idx++;
			}
		}

		LOG("Net will fill last row now");

			//t_values.back() = d;
		t_values.insert({ d, t_count - 1 });
		t_positions[t_count - 1] = d;
		auto& matrix_last_row = matrix.back();
		//matrix_last_row.front() = left_edge_condition(d);
		for (auto j{ std::next(matrix_last_row.begin()) }; j < std::prev(matrix_last_row.end()); j++)
		{
			*j = 0;
		}
		//matrix_last_row.back() = right_edge_condition(d);

		LOG("Constructor is done. Thank you forever")
	};

		

	double& Net::at(const double t, const double x)
	{
		//TODO: Dodaæ obs³ugê wyj¹tków
		return matrix.at(t_values.at(t)).at(x_values.at(x));
	}

	std::vector<std::vector<double>>* Net::getMatrix()
	{
		return &matrix;
	}

	void Net::dump(std::string filename)
	{
		std::ofstream file;
		file.open(filename);
		//file << ",";
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

	void Net::solve(SolvingMethod* method)
	{
		method->solveNet(this);
	}

	void Net::solveLeftEdgeCondition(unsigned int t_pos)
	{

		auto u0 = matrix.at(t_pos).at(0);
		auto u1 = matrix.at(t_pos).at(1);

		u0 = (1 / (left_edge_condition_function(t_positions[t_pos]) - left_edge_condition_derivative(t_positions[t_pos]))) * (((-left_edge_condition_derivative(t_positions[t_pos]) / h) * u1) - left_edge_condition_free_function(t_positions[t_pos]));
	}

	void Net::solveLeftEdgeCondition(double t)
	{
		this->solveLeftEdgeCondition(t_values.at(t));
	}

	void Net::solveRightEdgeCondition(unsigned int t_pos)
	{

		auto u0 = matrix.at(t_pos).at(0);
		auto u1 = matrix.at(t_pos).at(1);

		u0 = (1 / (right_edge_condition_function(t_positions[t_pos]) - right_edge_condition_derivative(t_positions[t_pos]))) * (((-right_edge_condition_derivative(t_positions[t_pos]) / h) * u1) - right_edge_condition_free_function(t_positions[t_pos]));
	}

	void Net::solveRightEdgeCondition(double t)
	{
		this->solveRightEdgeCondition(t_values.at(t));
	}
}

