#include "CN_GaussSiedel.h"
#include <cmath>
#include <iostream>

namespace MO
{
	namespace Crank_Nicolson
	{
		using namespace MO;

		

		void CN_GaussSiedel::solveNet(Net* net)
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

				//DEBUG
				std::vector<double> debugu{ 1, 2, 3, 4 };
				std::vector<double> debugd{ 1, 2, 3, 4, 5 };
				std::vector<double> debugl{ 1, 2, 3, 4 };
				std::vector<double> debugb{ 1, 2, 3, 4, 5 };
				std::vector<double> debugx{ 1, 2, 3, 4, 5 };
				auto debug_solutions = solveLinearEquation(debugu, debugd, debugl, debugb, debugx);



				std::cout << "Rozwi¹zywanie czasu " << current_time << std::endl;
				auto solutions = solveLinearEquation(u, d, l, b, matrix->at(past_position));

				matrix->at(current_position).swap(solutions);

			}
		}

		std::vector<double> CN_GaussSiedel::solveLinearEquation(std::vector<double>& u, std::vector<double>& d, std::vector<double>& l, std::vector<double>& b, std::vector<double> x0)
		{
			using my_matrix_type = std::vector<std::vector<double>>;

			const size_t matrix_size{ d.size() };

			my_matrix_type inv_ud;
			inv_ud.resize(matrix_size);
			//TODO: Po³¹czyæ ze sob¹
			{
				auto i{ inv_ud.begin() };
				size_t element_count{ matrix_size };
				while (i != inv_ud.end() )
				{
					i->reserve(element_count);
					i++;
					element_count--;
				}
			}
			{
				auto i{ inv_ud.begin() };
				auto id{ d.begin() };
				while (i != inv_ud.end())
				{
					i->push_back(1/(*id));
					i++;
					id++;
				}
			}
			for (size_t i{ 0 }; i < matrix_size; i++)
			{
				for (size_t j{ i + 1 }; j < matrix_size; j++)
				{
					auto prev{ -(inv_ud.at(i).back()) };
					inv_ud.at(i).push_back(prev * u.at(j - 1) * (1 / (d.at(j))));
				}
			}

			my_matrix_type M;
			/*M.resize(matrix_size - 1);
			{
				auto i{ M.begin() };
				auto il{ l.begin() };
				auto iinv_ud{ inv_ud.begin() };
				auto element_count{ matrix_size };

				while (i != M.end())
				{
					i->reserve(matrix_size);

					for (auto iinv_ud_element : *iinv_ud)
					{
						i->push_back((-(iinv_ud_element)) * *il);
					}

					i++;
					il++;
					iinv_ud++;
					element_count--;
				}

			}*/
			M.resize(matrix_size);
			{
				for (size_t i{ 0 }; i < matrix_size; i++)
				{
					M.at(i).reserve(i > 0 ? matrix_size - i + 1 : matrix_size);
					for (size_t j{ i > 0 ? i - 1 : 0 }; j < matrix_size - 1; j++)
					{
						M.at(i).push_back((-(inv_ud.at(i).at(i > 0 ? 0 : 1)) * (l.at(j))));
					}
					M.at(i).push_back(0);
				}

			}

			std::vector<double> x;
			x.reserve(matrix_size);

			std::vector<double> C;
			C.reserve(b.size());

			{
				auto iM{ M.begin() };
				//auto ix0{ std::next(x0.begin()) };
				auto iinv_ud{ inv_ud.begin() };
				//auto ib{ b.begin() };

				/*{
					double temp{ 0 };
					auto riinv_ud_element{ iinv_ud->rbegin() };
					auto rib{ b.rbegin() };

					while (riinv_ud_element != iinv_ud->rend())
					{
						temp += *riinv_ud_element * *rib;

						riinv_ud_element++;
						rib++;
					}

					C.push_back(temp);
				}

				iinv_ud++;

				x.push_back(C.back());*/

				while (iM != M.end())
				{
					{
						double temp{ 0 };
						auto riinv_ud_element{ iinv_ud->rbegin() };
						auto rib{ b.rbegin() };

						while (riinv_ud_element != iinv_ud->rend())
						{
							temp += *riinv_ud_element * *rib;

							riinv_ud_element++;
							rib++;
						}

						C.push_back(temp);
					}

					{
						double temp{ 0 };
						auto riM{ iM->rbegin() };
						auto rix0{ x0.rbegin() };

						while (riM != iM->rend())
						{
							temp += *riM * *rix0;

							riM++;
							rix0++;
						}

						x.push_back(temp + C.back());

					}

					iM++;
					iinv_ud++;
				}

			}

			std::vector<double> x_prev;
			x_prev.reserve(matrix_size);

			unsigned int iter{ 0 };
			const unsigned int maxIter{ 400 };
			double maxE{ -1 };
			const double etol{ 1e-15 };
			double maxRes{ -1 };
			const double restol{ 1e-15 };

			do
			{
				std::swap(x, x_prev);
				x.clear();

				{
					auto iM{ M.begin() };
					auto iC{ C.begin() };

					/*x.push_back(*iC);
					iC++;*/

					while (iM != M.end())
					{
						double temp{ 0 };
						auto riM{ iM->rbegin() };
						auto rix_prev{ x_prev.rbegin() };

						while (riM != iM->rend())
						{
							temp += *riM * *rix_prev;
								
							riM++;
							rix_prev++;
						}

						x.push_back(temp + *iC);

						iM++;
						iC++;
					}
				}

				iter++;

				maxE = -1;
				{
					auto ix_prev{ x_prev.begin() };
					auto ix{ x.begin() };
					while (ix != x.end())
					{
						double e{ std::abs(*ix - *ix_prev) };
						if (e > maxE) maxE = e;

						ix++;
						ix_prev++;
					}
				}

				maxRes = -1;
				{
					std::vector<double> Ax;
					Ax.reserve(matrix_size);

					auto iu{ u.begin() };
					auto id{ d.begin() };
					auto il{ l.begin() };
					auto ix1{ x.begin() };
					auto ix2{ std::next(x.begin()) };
					auto ix3{ std::next(x.begin(), 2) };

					Ax.push_back((*ix1** id) + (*ix2 * *iu));
					iu++;
					id++;
					while (iu != u.end())
					{
						Ax.push_back((*ix1 * *il) + (*ix2 * *id) + (*ix3 * *iu));

						iu++;
						id++;
						il++;
						ix1++;
						ix2++;
						ix3++;
					}
					Ax.push_back((*ix1 * *il) + (*ix2 * *id));

					auto iAx{ Ax.begin() };
					auto ib{ b.begin() };
					while (ib != b.end())
					{
						double res{ std::abs(*ib - *iAx) };
						if (res > maxRes) maxRes = res;

						ib++;
						iAx++;
					}
				}
				if (iter%10 == 0) std::cout << "\tIteracja " << iter << std::endl;
			} while ((iter < maxIter) && (maxE > etol) && (maxRes > restol));

			return x;

		}

	}
}