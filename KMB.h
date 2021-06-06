#pragma once
#include "SolvingMethod.h"
namespace MO
{

    class KMB :
        public SolvingMethod
    {
        double lambda;
    public:
        KMB(double lambda) : lambda(lambda) {};

        void solveNet(Net*);
    };

}
