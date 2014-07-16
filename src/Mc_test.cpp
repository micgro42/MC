/*
 * Mc_test.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: michael
 */

#define BOOST_TEST_MODULE mc_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <vector>
#include <cmath>
#include <algorithm> //std::sort
#include "Mc.h"

#define DEFINE_GLOBAL
#include "geom_pbc.c"

/**
 * @file Mc_test.cc
 *
 * @brief this file contains the unit test suite for the Mc class
 *
 *
 *
 *
 */



namespace logging = boost::log;

struct F1Point {

//     F() : i( 0 ) { std::cout << "setup" << std::endl; }
	F1Point(){
		testMc.setDelta(0.5);
        logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::info);
    }
    ~F1Point()          {  }

    Mc testMc;

};



BOOST_AUTO_TEST_SUITE (Mc_test)

BOOST_FIXTURE_TEST_CASE(set_get_test,F1Point){
	double kappa=5.1;
	double lambda = 1.1;
	double hReal=2.3;
	double hIm=3.3;
	double delta=0.7;


	testMc.setKappa(kappa);
	BOOST_CHECK_EQUAL(testMc.getKappa(),kappa);
	testMc.setKappa(-2*kappa);
	BOOST_CHECK_EQUAL(testMc.getKappa(),-2*kappa);

	testMc.setLambda(lambda);
	BOOST_CHECK_EQUAL(testMc.getLambda(),lambda);
	testMc.setLambda(-2*lambda);
	BOOST_CHECK_EQUAL(testMc.getLambda(),-2*lambda);

	testMc.setHReal(hReal);
	BOOST_CHECK_EQUAL(testMc.getHReal(),hReal);
	testMc.setHReal(-2*hReal);
	BOOST_CHECK_EQUAL(testMc.getHReal(),-2*hReal);

	testMc.setHIm(hIm);
	BOOST_CHECK_EQUAL(testMc.getHIm(),hIm);
	testMc.setHIm(-2*hIm);
	BOOST_CHECK_EQUAL(testMc.getHIm(),-2*hIm);

	testMc.setDelta(delta);
	BOOST_CHECK_EQUAL(testMc.getDelta(),delta);
	testMc.setDelta(-2*delta);
	BOOST_CHECK_EQUAL(testMc.getDelta(),-2*delta);

}

BOOST_FIXTURE_TEST_CASE(MagnetizationTest1,F1Point){
	double lambda = 4;
	double kappa=0.5;
	double h=0.3;
	ndim=2;
	int steps=6;
	lsize = (int *) malloc((ndim+1) * sizeof(int));
	for (int i=1; i<=ndim; ++i){
		lsize[i]=steps;
	}
	geom_pbc();
	BOOST_REQUIRE_EQUAL(nvol,36);
	vector<double> fieldReal,fieldIm;
	for (int i=1;i<=nvol;++i){
		fieldReal.push_back(testMc.getRandomUni());
		fieldIm.push_back(testMc.getRandomUni());
	}


	BOOST_REQUIRE_EQUAL(testMc.setFields(fieldReal,fieldIm),0);
	testMc.setLambda(lambda);
	testMc.setKappa(kappa);
	testMc.setHReal(h);
	testMc.setHIm(0.0);
	testMc.thermalizeField();


	int numUpdates = 100000;//0;
	vector<double> results;
	testMc.calculateMeanMagnetization(numUpdates, results);
	double tol=0.0002;
	BOOST_CHECK_CLOSE(testMc.getMeanRealMagnetisation(),0.8616,tol*100);
	BOOST_CHECK(abs(testMc.getMeanImMagnetisation())<testMc.getMeanImMagnetisationError()); // test that the imaginary part is 0 within the error



	free(lsize);


}

/*
 * BOOST_FIXTURE_TEST_CASE( Magnetization1PointLambda0, F1Point ){
    double tol=0.01;
    double delta=0.5;
    int exitcode=1;
    double BReal=0.3;
    double BIm=0.7;
    testMc.setLambda(0);
    testMc.setBReal(BReal);
    testMc.setBIm(BIm);
    exitcode=testMc.thermalizeField(delta);
    BOOST_REQUIRE(exitcode==0);
    int steps = 1000000;
    vector<double> results;
    exitcode=1;
    exitcode=testMc.calculateMeanMagnetization(steps, delta, results);
    BOOST_REQUIRE(exitcode==0);
    BOOST_CHECK_CLOSE( results.at(0) , BReal ,tol*100);//last variable is tolerance in percent
    BOOST_CHECK_CLOSE( results.at(1) , BIm ,tol*100);
    BOOST_CHECK_CLOSE( results.at(2) , 1+BReal*BReal+BIm*BIm ,tol*100);
}*/

/*
 BOOST_FIXTURE_TEST_CASE( Magnetization1PointLambda4, F1Point ){
    double tol=0.01;
    double delta=0.5;
    int exitcode=1;
    double lambda=4;
    double BReal=0.3;
    double BIm=0.7;
    testMc.setLambda(lambda);
    testMc.setBReal(BReal);
    testMc.setBIm(BIm);
    exitcode=testMc.thermalizeField(delta);
    BOOST_REQUIRE(exitcode==0);
    int steps = 1000000;
    vector<double> results;
    exitcode=1;
    exitcode=testMc.calculateMeanMagnetization(steps, delta, results);
    BOOST_REQUIRE(exitcode==0);
    BOOST_CHECK_CLOSE( results.at(0)+2*testMc.getLambda()*results.at(3) , BReal ,tol*100);
    BOOST_CHECK_CLOSE( results.at(1)+2*testMc.getLambda()*results.at(4) , BIm ,tol*100);
}*/

BOOST_AUTO_TEST_SUITE_END( )
