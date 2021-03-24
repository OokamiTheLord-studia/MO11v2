#pragma once
#ifdef _DEBUG_WITH_LOGS
#define LOG(x) std::cout << x << std::endl;
#else
#define LOG(x)
#endif // _DEBUG_WITH_LOGS