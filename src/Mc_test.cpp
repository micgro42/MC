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
    	cout << endl;
        logging::core::get()->set_filter (logging::trivial::severity >= logging::trivial::info);
    }
    ~F1Point()          {  }

    Mc testMc;

};


/*
BOOST_AUTO_TEST_SUITE (Mc_test)


BOOST_FIXTURE_TEST_CASE( matrixVectorLaplace_m0, F1Point ){

}
*/
