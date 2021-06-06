#include "KMB.h"

namespace MO
{
	void MO::KMB::solveNet(Net* net)
	{
		auto matrix = net->getMatrix();
		unsigned int t_count = matrix->size();
		unsigned int x_count = matrix->front().size();
		for (unsigned int t = 1; t < t_count; t++)
		{
			auto& current_time = matrix->at(t);
			auto& past_time = matrix->at(t - 1);
			for (unsigned int x = 1; x < x_count - 1; x++)
			{
				current_time.at(x) = (lambda * past_time.at(x - 1)) + (lambda * past_time.at(x + 1)) + ((1 - (2 * lambda)) * past_time.at(x));
			}

			net->solveLeftEdgeCondition(t);
			net->solveRightEdgeCondition(t);
		}
	}
}