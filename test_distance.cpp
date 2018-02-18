#include "catch.hpp"
#include "distance.h"

using std::vector;
using std::string;

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