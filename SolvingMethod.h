#pragma once
#include "net.h"

namespace MO
{

	class SolvingMethod
	{
		virtual void solveNet(MO::Net) = 0;
	};

}
