#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <string>

using std::cout;
using std::endl;
using std::vector;
using std::pow;
using std::sqrt;
using std::min_element;
using std::string;

double euclidean(const vector<double>& x1, const vector<double>& x2) {
    double d2 = 0.0;
    for (int i = 0; i < x1.size(); ++i) {
        d2 += pow(x1[i] - x2[i], 2);
    }
    return sqrt(d2);
}


double levenstein(const string& s1, const string& s2, double del, double ins, double sub) {
    // Edit distance of strings using user-specified costs
    // Source: http://www.python-course.eu/levenshtein_distance.php

    // For all i and j, dist[i,j] will contain the Levenshtein distance between the first i characters of a and the
    // first j characters of b
    auto n_rows = s1.size() + 1;
    auto n_cols = s2.size() + 1;
    // variable length array allowed in GCC
    double dist[n_rows][n_cols]; // don't need to 0-fill since it's filled from NW-> SE when running below


    // Source prefixes can be transformed into empty strings by deletions
    for (int row = 1; row < n_rows; ++row) {
        dist[row][0] = row * del;
    }

    // Target prefixes can be created from an empty source string by inserting the characters
    for (int col = 1; col < n_cols; ++col) {
        dist[0][col] = col * ins;
    }

    double cost;
    vector<double> subcost {0.0, 0.0, 0.0};
    for (int col = 1; col < n_cols; ++col) {
        for (int row = 1; row < n_rows; ++row) {
            if (s1[row-1] == s2[col-1]) {
                cost = 0.0;
            } else {
                cost = sub;
            }
            subcost[0] = dist[row-1][col] + del;
            subcost[1] = dist[row][col-1] + ins;
            subcost[2] = dist[row-1][col-1] + cost;  // substitution
            dist[row][col] = *min_element(subcost.begin(), subcost.end());
        }
    }

    return dist[n_rows-1][n_cols-1];
}

double levenstein(const string& s1, const string& s2) {
    // Edit distance of strings using default costs
    static const double del = 1.0, ins = 1.0, sub = 1.0;
    return levenstein(s1, s2, del, ins, sub);
}

int main() {
    // Test out Euclidean distance of 2 vectors of doubles
    vector<double> x1 {0.0, 0.0};
    vector<double> x2 {1.1, 2.2};

    double d = euclidean(x1, x2);

    cout << d << endl;

    // Test out edit (Levenstein) distance of 2 strings
    string s1 = "abc";
    string s2 = "bcd";
    string s3 = "";  // empty
    string s4 = "abcd";

//    cout << s1 << " " << s2 << endl;

    cout << levenstein(s1, s2) << endl;
    cout << levenstein(s1, s1) << endl;
    cout << levenstein(s1, s3) << endl;
    cout << levenstein(s1, s4) << endl;

    string t1 = "g5a23g5a3av34acra2ct";
    string t2 = "bacbzV#35v#V#nk";
    string t3 = "fdafdsafsafsa";

    cout << levenstein(t1, t2) << endl;  // 17
    cout << levenstein(t2, t3) << endl;  // 15
    cout << levenstein(t1, t2, 0.9, 1.1, 2.0) << endl;

    return 0;
}