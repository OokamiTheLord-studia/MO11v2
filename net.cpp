#include "net.h"
#include "simpleLogger.h"

namespace MO
{
	Net::Net(
		const double a
		, const double b
		, const double h
		, const double c
		, const double d
		, const double dt
		//sprawdzi� czy mog� by� const
		, std::function<double(double)> start_condition
		, std::function<double(double)> left_edge_condition
		, std::function<double(double)> right_edge_condition
	)
	{
		//TODO: Przemy�le� optymalizacj�
		//TODO: Asercja warto�ci

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

		LOG("Net reserved memory")

			//fill x
			auto& matrix_first_row = matrix.front();
		//x_values.front() = a;
		x_values.insert({ a, 0 });
		matrix_first_row.front() = start_condition(a);
		//mo�e da si� tu u�y� iteratora?
		{
			double i{ a + h };
			unsigned int idx{ 1 };
			//while (idx < x_values.size())
			while (idx < x_count - 1)
			{
				//x_values[idx] = i;
				x_values.insert({ i, idx });
				matrix_first_row[idx] = start_condition(i);

				i += h;
				idx++;
			}
		}
		//x_values.back() = b;
		x_values.insert({ b, x_count - 1 });
		matrix_first_row.back() = start_condition(b);

		LOG("Net filled first row")
			//fill edges and set middle to 0
			//t_values.front() = c;
			t_values.insert({ c, 0 });
		{
			double i{ c + dt };
			unsigned int idx{ 1 };
			//while (idx < t_values.size())
			while (idx < t_count - 1)
			{
				//t_values[idx] = i;
				t_values.insert({ i, idx });
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

		LOG("Net will fill last row now");

			//t_values.back() = d;
		t_values.insert({ d, t_count - 1 });
		auto& matrix_last_row = matrix.back();
		matrix_last_row.front() = left_edge_condition(d);
		for (auto j{ std::next(matrix_last_row.begin()) }; j < std::prev(matrix_last_row.end()); j++)
		{
			*j = 0;
		}
		matrix_last_row.back() = right_edge_condition(d);

		LOG("Constructor is done. Thank you forever")
	};

	double& Net::at(const double t, const double x)
	{
		//TODO: Doda� obs�ug� wyj�tk�w
		return matrix.at(t_values.at(t)).at(x_values.at(x));
	}
}
