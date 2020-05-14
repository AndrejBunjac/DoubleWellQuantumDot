#include "BasicFunctions.h"
#include <assert.h>


int BasicFunctions::KroneckerDelta(const int& n, const int& m)
{
	if (n == m) { return 1; }
	else		{ return 0;	}
}

int BasicFunctions::Factorial(int n)
{
	assert(n >= 0);
	return (n == 1 || n == 0) ? 1 : Factorial(n - 1) * n;
}

