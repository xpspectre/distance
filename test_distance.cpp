#include "catch.hpp"
#include "distance.h"

using std::vector;
using std::string;
using std::sqrt;

TEST_CASE( "Euclidean distance", "[euclidean_dist]" ) {
    vector<double> x1 { 0.0, 0.0 };
    vector<double> x2 { 1.0, 0.0 };
    vector<double> x3 { 1.0, 1.0 };

    REQUIRE( euclidean(x1, x2) == 1.0 );
    REQUIRE( euclidean(x2, x1) == 1.0 );

    REQUIRE( euclidean(x1, x3) == Approx(sqrt(2.0)));
    REQUIRE( euclidean(x3, x1) == Approx(sqrt(2.0)));
}

TEST_CASE( "Euclidean distance matrix", "[pdist_euclidean]" ) {
    vector<double> x1 { 0.0, 0.0 };
    vector<double> x2 { 1.0, 0.0 };
    vector<double> x3 { 1.0, 1.0 };
    vector<vector<double>> pts { x1, x2, x3 };

    SECTION( "vector form distance matrix" ) {
        vector<double> v = pdist(pts, euclidean);

        REQUIRE(v.size() == 3);

        vector<double> v_ = {1.0, sqrt(2.0), 1.0};
        REQUIRE(v == v_);  // Approx doesnt work for vectors, but in this case, it should be exactly the same calc running

        SECTION( "matrix form distance matrix" ) {
            vector<vector<double>> M = squareform(v);

            REQUIRE( M.size() == 3 );
            for(std::size_t i = 0; i < M.size(); ++i) {
                REQUIRE( M[i].size() == 3 );
            }

            vector<vector<double>> M_ { {0,1,sqrt(2.0)}, {1,0,1}, {sqrt(2.0),1,0} };
            REQUIRE( M == M_ );
        }
    }
}

TEST_CASE( "Levenstein distance", "[levenstein_dist]" ) {
    string s0 = "abc";
    string s1 = "bcd";
    string s2 = "";  // empty
    string s3 = "abcd";

    string t1 = "g5a23g5a3av34acra2ct";
    string t2 = "bacbzV#35v#V#nk";
    string t3 = "fdafdsafsafsa";

    SECTION( "simple cases that make up a distance mat" ) {
        REQUIRE( levenstein(s0, s1) == 2 );
        REQUIRE( levenstein(s0, s2) == 3 );
        REQUIRE( levenstein(s0, s3) == 1 );
        REQUIRE( levenstein(s1, s2) == 3 );
        REQUIRE( levenstein(s1, s3) == 1 );
        REQUIRE( levenstein(s2, s3) == 4 );
    }

    SECTION( "slight more complex strs" ) {
        REQUIRE( levenstein(t1, t2) == 17 );
        REQUIRE( levenstein(t2, t3) == 15 );
    }

    SECTION( "custom costs equal to default costs" ) {
        REQUIRE( levenstein(s0, s1, 1, 1, 1) == 2 );
        REQUIRE( levenstein(s0, s2, 1, 1, 1) == 3 );
        REQUIRE( levenstein(s0, s3, 1, 1, 1) == 1 );
    }

    SECTION( "new custom costs" ) {
        // costs are: del, ins, sub
        // note indel costs are equal
        REQUIRE( levenstein(s0, s1, 2, 2, 1) == 3 );  // cheaper to sub everything than indel ends
        REQUIRE( levenstein(s0, s2, 2, 2, 1) == 6 );
        REQUIRE( levenstein(s0, s3, 2, 2, 1) == 2 );
    }
}

// Test pdist and ppdist by making sure pdist works, and then making sure ppdist gives the same result

TEST_CASE( "Edit distance matrix calc", "[pdist_edit]" ) {
    string s0 = "abc";
    string s1 = "bcd";
    string s2 = "";  // empty
    string s3 = "abcd";
    vector<string> strs { s0, s1, s2, s3 };

    SECTION( "vector form distance matrix" ) {
        vector<double> v = pdist(strs, levenstein);

        REQUIRE( v.size() == 6 );

        vector<double> v_ = { 2, 3, 1, 3, 1, 4 };
        REQUIRE( v == v_ );

        SECTION( "matrix form distance matrix" ) {
            vector<vector<double>> M = squareform(v);

            REQUIRE( M.size() == 4 );
            for(std::size_t i = 0; i < M.size(); ++i) {
                REQUIRE( M[i].size() == 4 );
            }

            vector<vector<double>> M_ { {0,2,3,1}, {2,0,3,1}, {3,3,0,4}, {1,1,4,0} };
            REQUIRE( M == M_ );
        }
    }
}

TEST_CASE( "Parallel edit distance matrix calc", "[ppdist_edit]" ) {
    string s0 = "abc";
    string s1 = "bcd";
    string s2 = "";  // empty
    string s3 = "abcd";
    vector<string> strs { s0, s1, s2, s3 };

    SECTION( "vector form distance matrix" ) {
        vector<double> v = pdist(strs, levenstein);
        vector<double> vp = ppdist(strs, levenstein, 4);
        REQUIRE( v == vp );

        SECTION( "matrix form distance matrix" ) {
            vector<vector<double>> M = squareform(v);
            vector<vector<double>> Mp = squareform(vp);
            REQUIRE( M == Mp );
        }
    }
}