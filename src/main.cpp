#include<vector>
#include<iostream>
#include"Mc.h"
int main(int argc, char** argv){

	Mc test;
	vector<double> fieldReal,fieldIm,magnetization;
	fieldReal.push_back(0.5);
	fieldIm.push_back(0.5);
	cout << "set field exitcode: " << test.setFields(fieldReal,fieldIm) << endl;
	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;
	double acceptance=-1;
	cout << "createNewConfiguration exitcode: " << test.createNewConfiguration(1, 10, acceptance) << "; acceptance: " << acceptance << endl;
	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;


	for (int i=0; i<1000000; ++i){
		test.createNewConfiguration(1, 10, acceptance);
	}

	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;

	return 0;
}
