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
#include "global.h"

using namespace std;

#ifndef DEBUG
#define at(x) operator[](x)
#endif

class Mc {
public:
	Mc(int seed);
	Mc();
	virtual ~Mc();

	void setLambda(double lambdaIn){_lambda=lambdaIn;}
	void setKappa(double kappaIn){_kappa=kappaIn;}
	void setHReal(double hRealIn){_hReal=hRealIn;}
	void setHIm(double hImIn){_hIm=hImIn;}
	void setDelta(double delta){_delta=delta;}

	double getLambda(){return _lambda;}
	double getKappa(){return _kappa;}
	double getHReal(){return _hReal;}
	double getHIm(){return _hIm;}
	double getDelta(){return _delta;}
	double getMeanRealMagnetisation(){return _meanRealMagnetisation;}
	double getMeanRealMagnetisationError(){return _meanRealMagnetisationError;}
	double getMeanImMagnetisation(){return _meanImMagnetisation;}
	double getMeanImMagnetisationError(){return _meanImMagnetisationError;}
	double getMeanSquareAbsMagnetisation(){return _meanSquareAbsMagnetisation;}
	double getMeanSquareAbsMagnetisationError(){return _meanSquareAbsMagnetisationError;}





	double getRandomUni();
	void startRandomGenerator (double seed);
	int setFields(const vector<double> &fieldReal, const vector<double> &fieldIm);

	/** @brief calculates the magnetization of the entire field, does not change the field
	 *
	 * @param[out] magnetization the first component is the real part, the second component is the imaginary part
	 * @return 0 if everything went well,
	 */
	int calculateMagnetization(vector<double> &magnetization);

	double calculateP(int position, double realv, double imv);

	/**
	 * @brief creates a new configuration of the field with hitsPerPoint tried steps per point
	 *
	 *
	 *
	 *
	 */
	int createNewConfiguration(const double hitsPerPoint, double & acceptance);


	/**
	 * @brief thermalizes the field
	 *
	 *
	 */
	int thermalizeField();



	/**
	 * @brief creates steps new configurations
	 *
	 *
	 *
	 *
	 */
	int calculateMeanMagnetization(int steps, vector<double> & results);

private:

	///real part of the field
	vector<double> _fieldReal;

	///imaginary part of the field
	vector<double> _fieldIm;

	///step size for new configuration
	double _delta;

	///real component of the external field
	double _hReal;

	///imaginary component of the external field
	double _hIm;

	/// @brief interpolation between Gaussian model and XY-model
	///
	/// @ can be between 0 and +infinity
	double _lambda;

	/// coupling of neighbours
	double _kappa;

	double _meanRealMagnetisation;
	double _meanRealMagnetisationError;
	double _meanImMagnetisation;
	double _meanImMagnetisationError;

	/**
	 *
	 * magnetization.at(0)*magnetization.at(0)+magnetization.at(1)*magnetization.at(1)
	 *
	 * \f$ \left<\left|\phi_x\right|^2\right> = \left<\left|M\right|^2\right> = Re(M)^2 + Im(M)^2\f$
	 */
	double _meanSquareAbsMagnetisation;
	double _meanSquareAbsMagnetisationError;

	///random generator
	default_random_engine _randGenerator;
};

#endif /* MC_H_ */
