#pragma once
#include "SolvingMethod.h"
namespace MO
{

    class KMB :
        public SolvingMethod
    {
    public:
        void solveNet(MO::Net*);
    };

}
