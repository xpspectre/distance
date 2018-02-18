#include <iostream>
#include <vector>
#include <string>
#include "distance.h"
#include "util.h"

using std::cout;
using std::endl;
using std::vector;
using std::string;

int main() {
    // Test out Euclidean distance of 2 vectors of doubles
    cout << "Calculating Euclidean distances..." << endl;
    vector<double> x1 {0.0, 0.0};
    vector<double> x2 {1.1, 2.2};

    double d = euclidean(x1, x2);

    cout << d << endl;

    // Test out edit (Levenstein) distance of 2 strings
    cout << "Calculating edit distances..." << endl;
    string s0 = "abc";
    string s1 = "bcd";
    string s2 = "";  // empty
    string s3 = "abcd";

    // Same order as pdist's output below
    cout << levenstein(s0, s1) << endl;
    cout << levenstein(s0, s2) << endl;
    cout << levenstein(s0, s3) << endl;
    cout << levenstein(s1, s2) << endl;
    cout << levenstein(s1, s3) << endl;
    cout << levenstein(s2, s3) << endl;

    // Some more strings
//    cout << "Calculating more edit distances..." << endl;
//    string t1 = "g5a23g5a3av34acra2ct";
//    string t2 = "bacbzV#35v#V#nk";
//    string t3 = "fdafdsafsafsa";
//
//    cout << levenstein(t1, t2) << endl;  // 17
//    cout << levenstein(t2, t3) << endl;  // 15
//    cout << levenstein(t1, t2, 0.9, 1.1, 2.0) << endl;

    // Test out distance matrix calc
    cout << "Calculating distance matrix in vector form..." << endl;
    vector<string> strs {s0, s1, s2, s3};
    vector<double> dv = pdist(strs, levenstein); // CLion incorrectly says dist arg is an error
    for (double dvi : dv) {
        cout << dvi << endl;
    }

    cout << "Calculating distance matrix in vector form in parallel..." << endl;
    int n_threads = 4;
    vector<double> dvp = ppdist(strs, levenstein, n_threads);
    for (double dvi : dvp) {
        cout << dvi << endl;
    }

    // Test out squareform
//    vector<double> v {1, 2, 3, 4, 5, 6};
    // M should be:
    // 0 1 2 3
    // 1 0 4 5
    // 2 4 0 6
    // 3 5 6 0
//    auto M = squareform(v);
//    print_mat(M);

    auto dm = squareform(dv);
    print_mat(dm);

    // Tests on random strings
    cout << "Doing tests on random strings..." << endl;
    vector<string> xs = random_strs(15, 20);
    for (string x : xs) {
        cout << x << endl;
    }

    vector<double> xsv = ppdist(xs, levenstein, n_threads);
    vector<vector<double>> xsm = squareform(xsv);
    print_mat(xsm);

    return 0;
}