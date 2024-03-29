#include <iostream>
#include <math.h>
#include <iomanip>
#include "Forsythe.h"


void f(double t, double x[], double dxdt[]);
void rangeKuttaOrder3(double t, double tout, double Zn[], double h);


int main() {
	const double BEGIN_TIME = 0.0, END_TIME = 1.5, X1 = 2, X2 = 0.5;
	double step = 0.075;

	// rkf45
	rkf RKF;
	unsigned char work[6 * (2 * sizeof(Float)) + sizeof(struct rkf_inside)];
	double x[2] = { X1, X2 };
	RKF.f = f; RKF.neqn = 2; RKF.re = 0.0001; RKF.ae = 0.0001;
	RKF.work = work; RKF.flag = 1; RKF.Y = x; RKF.t = 0;

	std::cout.setf(std::ios::fixed);
	std::cout << "RKF45\n";
	std::cout << std::setw(1) << "t"
		<< std::setw(13) << "x1"
		<< std::setw(15) << "x2"
		<< std::setw(13) << "Flag"
		<< "\n";

	for (double h = BEGIN_TIME; h < END_TIME; h += step) {
		RKF.tout = h;
		rkf45(&RKF);
		std::cout << std::setw(5) << std::setprecision(3) << RKF.t
			<< std::setw(15) << std::setprecision(6) << x[0]
			<< std::setw(15) << std::setprecision(6) << x[1]
			<< std::setw(4) << RKF.flag
			<< "\n";
	}

	// Runge–Kutta method h = 0.075
	std::cout << "\n3rd order Runge Kutta method\n";
	std::cout << "h = " << step << "\n";
	std::cout << std::setw(1) << "t"
		<< std::setw(13) << "x1"
		<< std::setw(15) << "x2"
		<< "\n";

	for (double h = step; h <= END_TIME; h += step) {
		x[0] = X1; x[1] = X2;
		rangeKuttaOrder3(BEGIN_TIME, h, x, step);
		std::cout << std::setw(5) << std::setprecision(3) << h
			<< std::setw(15) << std::setprecision(6) << x[0]
			<< std::setw(15) << std::setprecision(6) << x[1]
			<< "\n";
	}

	// Runge–Kutta method h - mine
	step  = 0.00075;
	std::cout << "\nh = " << step << "\n";
	std::cout << std::setw(1) << "t"
		<< std::setw(13) << "x1"
		<< std::setw(15) << "x2"
		<< "\n";

	int k = 0;
	for (double h = step; h < END_TIME + step; h += step) {
		k++;
		x[0] = X1; x[1] = X2;
		rangeKuttaOrder3(BEGIN_TIME, h, x, step);
		if (k == 100) {
			std::cout << std::setw(5) << std::setprecision(3) << h
				<< std::setw(15) << std::setprecision(6) << x[0]
				<< std::setw(15) << std::setprecision(6) << x[1]
				<< std::endl;
			k = 0;
		}
	}
	return 0;
}


void f(double t, double x[], double dxdt[]) {
	dxdt[0] = -14 * x[0] + 13 * x[1] + cos(1 + t);
	dxdt[1] = 20 * x[0] - 30 * x[1] + atan(1 + pow(t, 2));
}


void rangeKuttaOrder3(double BEGIN_TIME, double END_TIME, double zn[], double h) {
	double k1[2], k2[2], k3[2], dxdt[2], x[2], tn;

	for (double t = BEGIN_TIME; t <= END_TIME; t += h) {
		//k1
		tn = t;
		x[0] = zn[0]; x[1] = zn[1];
		f(tn, x, dxdt);
		k1[0] = h * dxdt[0]; k1[1] = h * dxdt[1];

		//k2
		tn = t + h / 2;
		x[0] = zn[0] + k1[0] / 2; x[1] = zn[1] + k1[1] / 2;
		f(tn, x, dxdt);
		k2[0] = h * dxdt[0]; k2[1] = h * dxdt[1];

		//k3
		tn = t + 0.75 * h;
		x[0] = zn[0] + 0.75 * k2[0]; x[1] = zn[1] + 0.75 * k2[1];
		f(tn, x, dxdt);
		k3[0] = h * dxdt[0]; k3[1] = h * dxdt[1];

		//Zn+1
		zn[0] = zn[0] + (2 * k1[0] + 3 * k2[0] + 4 * k3[0]) / 9; 
		zn[1] = zn[1] + (2 * k1[1] + 3 * k2[1] + 4 * k3[1]) / 9;
	}
}