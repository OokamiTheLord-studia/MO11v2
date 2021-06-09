#pragma once
#include "SolvingMethod.h"

namespace MO
{

    namespace Crank_Nicolson
    {
        using namespace MO;

        class CN_GaussSiedel :
            public SolvingMethod
        {
        public:
            void solveNet(Net*) override;

        };

    }
}