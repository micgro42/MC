/*
 * Mc_test.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: michael
 */

#define BOOST_TEST_MODULE kongrad_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/log/core.hpp>
#include <boost/log/trivial.hpp>
#include <boost/log/expressions.hpp>
#include <vector>
#include <cmath>
#include <algorithm> //std::sort
#include "Mc.h"

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
            vector<double> fieldReal,fieldIm;
            fieldReal.push_back(0.5);
            fieldIm.push_back(2);
            testMc.setFields(fieldReal,fieldIm);
    	cout << endl;
        logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::info);
    }
    ~F1Point()          {  }

    Mc testMc;

};



BOOST_AUTO_TEST_SUITE (Mc_test)


BOOST_FIXTURE_TEST_CASE( Magnetization1PointLambda0, F1Point ){
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
    BOOST_CHECK_CLOSE( results.at(0) , BReal ,tol*100);
    BOOST_CHECK_CLOSE( results.at(1) , BIm ,tol*100);
    BOOST_CHECK_CLOSE( results.at(2) , 1+BReal*BReal+BIm*BIm ,tol*100);
}

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
}

BOOST_AUTO_TEST_SUITE_END( )
