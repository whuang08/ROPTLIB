
#ifndef RANDGEN_CPP
#define RANDGEN_CPP

#include "Others/randgen.h"

void genrandseed(unsigned int s)
{
    srand(s);
}

realdp genrandreal(void)
{
    return static_cast<realdp> (rand()) / RAND_MAX;
}

realdp genrandnormal(void)
{
	static realdp rand1, rand2;
	realdp tmp = genrandreal();
	while (tmp == 1.0)
		tmp = genrandreal();
    rand1 = static_cast<realdp> (-2) * log(static_cast<realdp> (1.0) - tmp);
    rand2 = (static_cast<realdp> (1) - genrandreal()) * static_cast<realdp> (6.2831853071795864769252866);
	return sqrt(rand1) * cos(rand2);
}

#endif
