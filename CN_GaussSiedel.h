#pragma once
#include "SolvingMethod.h"
#include <vector>

namespace MO
{

    namespace Crank_Nicolson
    {
        using namespace MO;

        class CN_GaussSiedel :
            public SolvingMethod
        {

            const double lambda;

            std::vector<double> solveLinearEquation(std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>);

        public:
            void solveNet(Net*) override;
            CN_GaussSiedel(double lambda) : lambda(lambda) {};

        };

    }
}