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


	cout << "set field exitcode: " << test.setFields(fieldReal,fieldIm) << endl;
	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;
	double acceptance=-1;
	cout << "createNewConfiguration exitcode: " << test.createNewConfiguration(10, acceptance) << endl;
	cout << "acceptance: " << acceptance << endl;
	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;


	test.setLambda(4);
	test.thermalizeField();
	int numUpdates = 100000;//0;
	vector<double> results;
	test.calculateMeanMagnetization(numUpdates, results);

	free(lsize);
	return 0;
}
