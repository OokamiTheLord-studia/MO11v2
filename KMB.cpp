#include "KMB.h"

namespace MO
{
	void MO::KMB::solveNet(Net* net)
	{
		auto matrix = net->getMatrix();
		const size_t t_count = matrix->size();
		const size_t x_count = matrix->front().size();
		for (size_t t = 1; t < t_count; t++)
		{
			auto& current_time = matrix->at(t);
			auto& past_time = matrix->at(t - 1);
			for (size_t x = 1; x < x_count - 1; x++)
			{
				current_time.at(x) = (lambda * past_time.at(x - 1)) + (lambda * past_time.at(x + 1)) + ((1 - (2 * lambda)) * past_time.at(x));
			}

			net->solveLeftEdgeCondition(t);
			net->solveRightEdgeCondition(t);
		}
	}
}