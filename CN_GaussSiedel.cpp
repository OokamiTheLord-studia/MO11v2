#include "CN_GaussSiedel.h"

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
				auto solutions = solveLinearEquation(u, d, l, b, matrix->front());

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
					i->push_back(*id);
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
			M.resize(matrix_size - 1);
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

				{
					double temp{ 0 };
					auto riinv_ud_element{ iinv_ud->rbegin() };
					auto rib{ b.rbegin() };

					while (riinv_ud_element != iinv_ud->rend())
					{
						temp += *riinv_ud_element * *rib;
					}

					C.push_back(temp);
				}

				iinv_ud++;

				x.push_back(C.back());

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

			do
			{
				std::swap(x, x_prev);
				x.clear();



			} while ()

		}

	}
}