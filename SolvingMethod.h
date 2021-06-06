#pragma once
#include "net.h"

namespace MO
{
	class Net;

	class SolvingMethod
	{
	public:
		virtual void solveNet(Net*) = 0;
	};

}
