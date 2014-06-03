/*
 * Mc.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: Michael Große
 */

#include "Mc.h"
#include <cmath>

Mc::Mc(int seed) {
	// TODO Auto-generated constructor stub
	startRandomGenerator(seed);
	_BReal=0.3;
	_BIm=0.7;
	_lambda=0;

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
		for (int i=0;i<_fieldReal.size();++i){
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

double Mc::calculateP(double realv, double imv){
	double lambda=0.;
	double phisquare=pow(realv,2)+pow(imv,2);
	return exp(2*(_BReal*realv+_BIm*imv)-(phisquare-lambda*pow(phisquare-1,2)));
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

    for (int i=0;i<hitsPerPoint;++i){
    rReal = getRandomUni()*delta*2-delta;
    rIm = getRandomUni()*delta*2-delta;
    rAccept = getRandomUni();
    p = calculateP(_fieldReal.at(0),_fieldIm.at(0));
    pnew=calculateP(_fieldReal.at(0)+rReal,_fieldIm.at(0)+rIm);
//	cout << "p " << p << endl;

	if ((pnew>p)||(pnew/p>rAccept)){
		_fieldReal.at(0)=_fieldReal.at(0)+rReal;
		_fieldIm.at(0)=_fieldIm.at(0)+rIm;
		++numAccepts;
	}
	++numHits;
    }

	acceptance=(double)numAccepts/numHits;
	return exitcode;
}






