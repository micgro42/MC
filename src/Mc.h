/*
 * Mc.h
 *
 *  Created on: Jun 3, 2014
 *      Author: Michael Große
 */

#ifndef MC_H_
#define MC_H_

#include <vector>
#include <random>
#include <cmath>
#include <iostream>

using namespace std;

#ifndef DEBUG
#define at(x) operator[](x)
#endif

class Mc {
public:
	Mc();
	virtual ~Mc();
	double getRandomUni();
	void startRandomGenerator (double seed);
	int setFields(const vector<double> &fieldReal, const vector<double> &fieldIm);

	/** @brief calculates the magnetization per configuration and in the thermodynamic mean
	 *
	 * @param[out] magnetization the first component is the real part, the second component is the imaginary part
	 * @return 0 if everything went well,
	 */
	int calculateMagnetization(vector<double> &magnetization);


	int createNewConfiguration(double delta, double hitsPerPoint, double &acceptance);

private:

	///real part of the field
	vector<double> _fieldReal;

	///imaginary part of the field
	vector<double> _fieldIm;

	///random generator
	default_random_engine _randGenerator;
};

#endif /* MC_H_ */
