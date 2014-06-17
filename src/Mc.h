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
	Mc(int seed);
	Mc();
	virtual ~Mc();

	void setBReal(double breal){_BReal=breal;}
	void setBIm(double bim){_BIm=bim;}
	void setLambda(double lambdaIn){_lambda=lambdaIn;}
	void setKappa(double kappaIn){_kappa=kappaIn;}

	double getBReal(){return _BReal;}
	double getBIm(){return _BIm;}
	double getLambda(){return _lambda;}

	double getRandomUni();
	void startRandomGenerator (double seed);
	int setFields(const vector<double> &fieldReal, const vector<double> &fieldIm);

	/** @brief calculates the magnetization of the entire field, does not change the field
	 *
	 * @param[out] magnetization the first component is the real part, the second component is the imaginary part
	 * @return 0 if everything went well,
	 */
	int calculateMagnetization(vector<double> &magnetization);

	double calculateP(double realv, double imv);

	/**
	 * @brief creates a new configuration of the field with hitsPerPoint tried steps per point
	 *
	 *
	 *
	 *
	 */
	int createNewConfiguration(const double delta, const double hitsPerPoint, double & acceptance);


	/**
	 * @brief thermalizes the field
	 *
	 *
	 */
	int thermalizeField(double &delta);



	/**
	 * @brief creates steps new configurations
	 *
	 *
	 *
	 *
	 */
	int calculateMeanMagnetization(int steps, const double delta, vector<double> & results);

private:

	///real part of the field
	vector<double> _fieldReal;

	///imaginary part of the field
	vector<double> _fieldIm;

	///real part of a value used for test with 1-point fields
	double _BReal;

	///real part of a value used for test with 1-point fields
	double _BIm;

	/// @brief interpolation between Gaussian model and XY-model
	///
	/// @ can be between 0 and +infinity
	double _lambda;

	/// coupling of neighbours
	double _kappa;

	///random generator
	default_random_engine _randGenerator;
};

#endif /* MC_H_ */
