#pragma once
#include <cmath>

#define M_PI1         3.141592653589793238462643383279502884L /* pi */

namespace Parameters {

	//Potential parameters
	static double constexpr omega = 1.0;
	static double constexpr omega0 = 0.5;
	static double constexpr distance = 7.0;

	//Quantum number setup
	static int constexpr m1 = 0;
	static int constexpr m2 = 0;
	static int am1 = std::abs(m1);
	static int am2 = std::abs(m2);
	static int constexpr mmax = 4;

	const int nzmax = 3;
	const int nrmax = 3;
};