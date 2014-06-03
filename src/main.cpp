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


	double meanAcceptance=acceptance;
	int iterthermal=0;
	int deltaNotChanged=0;
	bool deltaChanged;
	int thermalStep=100;
	while (deltaNotChanged<10){
		deltaChanged=false;
		meanAcceptance=0;
		for (int i=0; i<thermalStep; ++i){
			test.createNewConfiguration(delta, 10, acceptance);
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
		cout << "mean Acceptance: " << meanAcceptance << " delta " << delta << endl;
		if (iterthermal>100*thermalStep){
			deltaNotChanged=20;
		}
	}
	cout << "delta: " << delta << " after " << iterthermal << " steps" << endl;

	double meanMagReal=0;
	double meanMagIm=0;
	double meansqrabsMag=0;
	int steps = 1000000;
	for (int i=0; i<steps; ++i){
		test.createNewConfiguration(delta, 10, acceptance);
		test.calculateMagnetization(magnetization);
		meanMagReal+=magnetization.at(0);
		meanMagIm+=magnetization.at(1);
		meansqrabsMag+=magnetization.at(0)*magnetization.at(0)+magnetization.at(1)*magnetization.at(1);
	}

	cout << "mean real component: " << meanMagReal/steps << endl;
	cout << "mean imaginary component: " << meanMagIm/steps << endl;
	cout << "absolute magnetization squared " << meansqrabsMag/steps << endl;

	cout << "calculateMagnetization exitcode: " << test.calculateMagnetization(magnetization) << endl;
	cout << "Real part: " << magnetization.at(0) << "; Imaginary part: " << magnetization.at(1) << endl;

	return 0;
}
