#include<vector>
#include<iostream>
#include"Mc.h"
#define DEFINE_GLOBAL
#include "geom_pbc.c"


int main(int argc, char** argv){
	int seed;
	if (argc > 1){
		seed=atoi(argv[1]);
	}else{
		seed=time(NULL);
	}

	ndim=2;//atoi(argv[4]);
	int steps=6;//atoi(argv[2]);
	lsize = (int *) malloc((ndim+1) * sizeof(int));
	for (int i=1; i<=ndim; ++i){
		lsize[i]=steps;
	}
	geom_pbc();
	cout << "nvol " << nvol << endl;

	Mc test(seed);
	vector<double> fieldReal,fieldIm;
    for (int i=1;i<=nvol;++i){
    	fieldReal.push_back(test.getRandomUni());
    	fieldIm.push_back(test.getRandomUni());
    }

	vector<double> magnetization;

	double delta=0.5;
	cout << "set field exitcode: " << test.setFields(fieldReal,fieldIm) << endl;
	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;
	double acceptance=-1;
	cout << "createNewConfiguration exitcode: " << test.createNewConfiguration(delta, 10, acceptance) << endl;
	cout << "acceptance: " << acceptance << endl;
	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;


	test.setLambda(4);
	test.thermalizeField(delta);
	int numUpdates = 100000;//0;
	vector<double> results;
	test.calculateMeanMagnetization(numUpdates, delta, results);


	cout << "mean real component: " << results.at(0) << endl;
	cout << "mean imaginary component: " << results.at(1) << endl;
	cout << "mean absolute magnetization squared is " << results.at(2) << " and should be " << 1+test.getBReal()*test.getBReal()+test.getBIm()*test.getBIm() << endl;
	cout << "if lambda is nonzero:" << endl;
	cout << results.at(0)+2*test.getLambda()*results.at(3) << " should be " << test.getBReal() << endl;
	cout << results.at(1)+2*test.getLambda()*results.at(4) << " should be " << test.getBIm() << endl;

	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;


	free(lsize);
	return 0;
}
