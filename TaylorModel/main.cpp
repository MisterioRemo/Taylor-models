#include <iostream>
#include "odu.h"
using std::cout;
using std::endl;

int main() {
	// /////////////////////////////////////////////////////////////////
	// u' = v
	// v' = u^2
	// u(0) = 1 + a
	// v(0) = 1 + b
	// a, b = [-0.05; 0.05]
	// t = [0; 6]
	// /////////////////////////////////////////////////////////////////

	vector<interval<double> > initPoint;
	initPoint.push_back(interval<double>(0.95, 1.05));
	initPoint.push_back(interval<double>(-1.05, -0.95));

	equation<double> odu(2, 0, 18);
	odu.initialFlow(&initPoint);

	odu.RungeKutta(0, 6, 0.01);

	return 0;
}
