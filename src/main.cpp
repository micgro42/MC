#include<vector>
#include<iostream>
#include"Mc.h"
int main(int argc, char** argv){
	int seed;
	if (argc > 1){
		seed=atoi(argv[1]);
	}else{
		seed=time(NULL);
	}
	Mc test(seed);
	vector<double> fieldReal,fieldIm,magnetization;
	fieldReal.push_back(0.5);
	fieldIm.push_back(0.5);
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
	int steps = 1000000;
	vector<double> results;
	test.calculateMeanMagnetization(steps, delta, results);


	cout << "mean real component: " << results.at(0) << endl;
	cout << "mean imaginary component: " << results.at(1) << endl;
	cout << "mean absolute magnetization squared is " << results.at(2) << " and should be " << 1+test.getBReal()*test.getBReal()+test.getBIm()*test.getBIm() << endl;

	cout << results.at(0)+2*test.getLambda()*results.at(3) << endl;
	cout << results.at(1)+2*test.getLambda()*results.at(4) << endl;

	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;

	return 0;
}
