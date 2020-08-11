
#include "Others/Timer.h"

/*Define the namespace*/
namespace ROPTLIB{

	unsigned long getTickCount(void)
	{
		unsigned long currentTime = 0;

#ifdef _WIN64
		LARGE_INTEGER li;
		QueryPerformanceFrequency(&li);
		long long dff = li.QuadPart;
		QueryPerformanceCounter(&li);
		currentTime = static_cast<unsigned long> (li.QuadPart * 1000000 / dff);
#elif _WIN32
		LARGE_INTEGER li;
		QueryPerformanceFrequency(&li);
		long long dff = li.QuadPart;
		QueryPerformanceCounter(&li);
		currentTime = static_cast<unsigned long> (li.QuadPart * 1000000 / dff);
#elif __APPLE__
#include "TargetConditionals.h"
#if TARGET_OS_MAC
		struct timeval current;
		gettimeofday(&current, NULL);
		currentTime = current.tv_sec * 1000000 + current.tv_usec;
#endif
#elif __linux
		struct timeval current;
		gettimeofday(&current, NULL);
		currentTime = current.tv_sec * 1000000 + current.tv_usec;
#endif

		return currentTime;
	}
}; /*end of ROPTLIB namespace*/
