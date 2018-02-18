#include "catch.hpp"
#include "distance.h"

using std::vector;
using std::string;

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

TEST_CASE( "Factorials are computed", "[factorial]" ) {
//    REQUIRE( Factorial(0) == 1 );
    REQUIRE( Factorial(1) == 1 );
    REQUIRE( Factorial(2) == 2 );
    REQUIRE( Factorial(3) == 6 );
    REQUIRE( Factorial(10) == 3628800 );
}

TEST_CASE( "Euclidean distance", "[euclidean_dist]" ) {
    vector<double> x1 {0.0, 0.0};
    vector<double> x2 {1.0, 0.0};
    REQUIRE( euclidean(x1, x2) == 1.0 );
}

TEST_CASE( "Levenstein distance", "[levenstein_dist]" ) {
    string s0 = "abc";
    string s1 = "bcd";
    REQUIRE( levenstein(s0, s1) == 2 );
}