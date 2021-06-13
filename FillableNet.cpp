#include "FillableNet.h"
#include "simpleLogger.h"
#include <fstream>

namespace MO
{
	FillableNet::FillableNet(
		const double a
		, const double b
		, const double h
		, const double c
		, const double d
		, const double dt
		//sprawdziæ czy mog¹ byæ const
		, std::function<double(double, double)> filling_function
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
			//x_values.resize(x_count);
		//t_values.resize(t_count);
			matrix.resize(t_count);
		for (auto& i : matrix)
		{
			i.resize(x_count);
		}

		LOG("Net reserved memory")

			//fill x
		//x_values.front() = a;
			x_values.insert({ a, 0 });
		//mo¿e da siê tu u¿yæ iteratora?
		{
			double i{ a + h };
			unsigned int idx{ 1 };
			//while (idx < x_values.size())
			while (idx < x_count - 1)
			{
				//x_values[idx] = i;
				x_values.insert({ i, idx });

				i += h;
				idx++;
			}
		}
		//x_values.back() = b;
		x_values.insert({ b, x_count - 1 });

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


				i += dt;
				idx++;
			}
		}

		//t_values.back() = d;
		t_values.insert({ d, t_count - 1 });

		//fill matrix
		for (auto t = t_values.cbegin(); t != t_values.cend(); t++)
		{
			for (auto x = x_values.cbegin(); x != x_values.cend(); x++)
			{
				matrix.at(t->second).at(x->second) = filling_function(t->first, x->first);
			}
		}

		LOG("Constructor is done. Thank you forever")
	};

	double& FillableNet::at(const double t, const double x)
	{
		//TODO: Dodaæ obs³ugê wyj¹tków
		return matrix.at(t_values.at(t)).at(x_values.at(x));
	}

	void FillableNet::dump(std::string filename)
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

	void FillableNet::dumpForT(std::string filename, size_t t)
	{
		std::ofstream file;
		file.open(filename);
		for (auto ix_values{ x_values.cbegin() }; ix_values != x_values.cend(); ix_values++)
		{
			file << ix_values->first << "," << matrix.at(t).at(ix_values->second) << std::endl;
		}
		file.close();
	}

	void FillableNet::dumpForT(std::string filename, double t)
	{
		dumpForT(filename, static_cast<size_t>(t_values.at(t)) );
	}

	std::vector<std::vector<double>>* FillableNet::getMatrix()
	{
		return &matrix;
	}
}