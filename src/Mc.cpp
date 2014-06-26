/*
 * Mc.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: Michael Große
 */

#include "Mc.h"
#include <cmath>
#include "stat5.h"

Mc::Mc(int seed) {
	// TODO Auto-generated constructor stub
	startRandomGenerator(seed);
	_BReal=0.3;
	_BIm=0.7;
	_lambda=4;
	_kappa=0.5;
	_hReal=0.3;
	_hIm=0;

}

Mc::Mc() {
	// TODO Auto-generated constructor stub
	startRandomGenerator(0);
	_BReal=0.3;
	_BIm=0.7;
	_lambda=4;
	_kappa=0.5;
	_hReal=0.3;
	_hIm=0;

}

Mc::~Mc() {
	// TODO Auto-generated destructor stub
}

double Mc::getRandomUni(){
    uniform_real_distribution<double> distributiond(0.0,1.0); // to generate the values
    return distributiond(_randGenerator);
}

void Mc::startRandomGenerator (double seed){
    _randGenerator.seed(seed);
}

///@return 0 if everything went ok, 83 if the fields have different dimensions; 84 if the fields are of 0 length;
///@todo extend to several dim using geom_pbc.c
int Mc::setFields(const vector<double> &fieldReal, const vector<double> &fieldIm){
	int exitcode=1;
	if (fieldReal.size()!=fieldIm.size()){
		exitcode=83;
	}else{
		if(fieldReal.size()>0){
			exitcode=0;
			_fieldReal=fieldReal;
			_fieldIm=fieldIm;
		}else{
			exitcode=84;
		}
	}
	return exitcode;
}



int Mc::calculateMagnetization(vector<double> &magnetization){
	int exitcode=1;
	if (_fieldReal.size()!=0){
		magnetization.assign(2,0);
		for (unsigned int i=0;i<_fieldReal.size();++i){
			magnetization.at(0)+=_fieldReal.at(i);
			magnetization.at(1)+=_fieldIm.at(i);
		}
		magnetization.at(0)=magnetization.at(0)/(double)_fieldReal.size();
		magnetization.at(1)=magnetization.at(1)/(double)_fieldReal.size();
		exitcode=0;
	}else{
		exitcode=84;
	}
	return exitcode;
}

///see equation 4.13 and 4.14
double Mc::calculateP(int position, double realv, double imv){
	double BReal=0;
	double BIm=0;
	for (int i=1;i<=ndim;++i){
		BReal+=_fieldReal.at(nn[i][position])+_fieldReal.at(nn[i+ndim][position]);
		BIm+=_fieldIm.at(nn[i][position])+_fieldIm.at(nn[i+ndim][position]);
	}
	BReal=_hReal+_kappa*BReal;
	BIm=_hIm+_kappa*BIm;
	double phisquare=pow(realv,2)+pow(imv,2);
	return exp(2*(BReal*realv+BIm*imv)-phisquare-_lambda*pow(phisquare-1,2));
}

int Mc::createNewConfiguration(const double delta, const double hitsPerPoint, double &acceptance){
	int exitcode=1;
	int numAccepts=0;
	int numHits=0;
	double rReal;
	double rIm;
	double rAccept;
	double p;
	double pnew;
	for (int i=0; i<nvol;++i){
		for (int j=0;j<hitsPerPoint;++j){
			rReal = getRandomUni()*delta*2-delta;
			rIm = getRandomUni()*delta*2-delta;
			rAccept = getRandomUni();
			p = calculateP(i,_fieldReal.at(i),_fieldIm.at(i));
			pnew=calculateP(i,_fieldReal.at(i)+rReal,_fieldIm.at(i)+rIm);
//			cout << "i " << i << "; j " << j << "; p " << p << endl;

			if ((pnew>p)||(pnew/p>rAccept)){
				_fieldReal.at(i)=_fieldReal.at(i)+rReal;
				_fieldIm.at(i)=_fieldIm.at(i)+rIm;
				++numAccepts;
			}
			++numHits;
		}
	}

	acceptance=(double)numAccepts/numHits;
	return exitcode;
}

///@todo check/make sure that fields are initialised 
int Mc::thermalizeField(double & delta){
	double acceptance;
	double meanAcceptance=0;
	int iterthermal=0;
	int deltaNotChanged=0;
	bool deltaChanged;
	int thermalStep=100;
	while (deltaNotChanged<10){
		deltaChanged=false;
		meanAcceptance=0;
		for (int i=0; i<thermalStep; ++i){
			createNewConfiguration(delta, 10, acceptance);
			meanAcceptance+=acceptance;
		}
		meanAcceptance=meanAcceptance/(double)thermalStep;
		if (0.3>meanAcceptance){
			delta=delta*0.95;
			deltaChanged=true;
		}
		if (meanAcceptance>0.5){
			delta=delta*1.05;
			deltaChanged=true;
		}
		if (deltaChanged){
			deltaNotChanged=0;
		}else{
			++deltaNotChanged;
		}
		iterthermal+=thermalStep;
		if (iterthermal>100*thermalStep){
			deltaNotChanged=20;
		}
	}
	return 0;
}

/**
 * @param[out] results Component 0: Real part of the mean magnetisation
 * Component 1: Imaginary part of the mean magnetisation
 * Component 2: mean squared absolute magnetisation
 * Component 3: real part of foo
 * Component 4: imaginary part of foo
 *
 *
 * @details foo \f$ =\left< \Phi \left( \left|\Phi\right|^2 -1 \right) \right> \f$
 */
int Mc::calculateMeanMagnetization(int steps, const double delta, vector<double> & results){
	results.assign(5,0);
	clear5(10,500);
	double meanMagReal=0;
	double meanMagIm=0;
	double meansqrabsMag=0;
	double fooReal=0;
	double fooIm=0;
	double acceptance;
	vector<double> magnetization;

	for (int i=0; i<steps; ++i){
		createNewConfiguration(delta, 10, acceptance);
		calculateMagnetization(magnetization);
		meanMagReal+=magnetization.at(0);
		meanMagIm+=magnetization.at(1);
		double absmagsquare=magnetization.at(0)*magnetization.at(0)+magnetization.at(1)*magnetization.at(1);
		accum5(1,magnetization.at(0));
		accum5(2,magnetization.at(1));
		accum5(3,absmagsquare);
		meansqrabsMag+=absmagsquare;
		fooReal+=magnetization.at(0)*(absmagsquare-1);
		fooIm+=magnetization.at(1)*(absmagsquare-1);
	}
	results.at(0)=meanMagReal/steps;
	results.at(1)=meanMagIm/steps;
	results.at(2)=meansqrabsMag/steps;
	results.at(3)=fooReal/steps;
	results.at(4)=fooIm/steps;
	cout << "_BReal: " << _BReal << endl;
	cout << "_BIm: " << _BIm << endl;
	cout << "1+|B|^2: " << 1+_BReal*_BReal+_BIm*_BIm << endl;
	cout << "stat5: mean Re(magnetization) " << aver5(1) << " +- " << sigma5(1) << endl;
	cout << "stat5: mean Im(magnetization) " << aver5(2) << " +- " << sigma5(2) << endl;
	cout << "stat5: mean abs(magnetization) " << aver5(3) << " +- " << sigma5(3) << endl;


	cout << "stat5: covar5 Re, Im " << covar5(1,2) << endl;
	cout << "stat5: tau5 Re(magnetization) " << tau5(1) << endl;
	cout << "stat5: tau5 Im(magnetization) " << tau5(2) << endl;
	cout << "stat5: tauint5 Re(magnetization) " << tauint5(1) << endl;
	cout << "stat5: tauint5 Im(magnetization) " << tauint5(2) << endl;


	return 0;
}




