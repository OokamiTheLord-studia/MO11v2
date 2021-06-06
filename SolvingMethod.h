#pragma once
#include "net.h"

namespace MO
{

	class SolvingMethod
	{
	public:
		virtual void solveNet(MO::Net*) = 0;
	};

}
