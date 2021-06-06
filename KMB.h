#pragma once
#include "SolvingMethod.h"
namespace MO
{

    class KMB :
        public SolvingMethod
    {
        void solveNet(MO::Net);
    };

}
