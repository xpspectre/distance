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
    // Note: doesn't use default params since function pointers passed as args into other functions don't recognize them
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

template <class T>
vector<double> pdist(const vector<T> &x, double (*dist)(const T &, const T &)) {
    // Distance matrix between all pairs in x, in order of x, in vector form (lower triangle entries), using specified
    // distance function
    auto nx = x.size();
    vector<double> v (nx*(nx-1)/2, 0.0);
    unsigned long ind;
    for (int i = 1; i < nx; ++i) {
        for (int j = 0; j < i; ++j) {
            ind = nx*j - j*(j+1)/2 + i - 1 - j;
            v[ind] = dist(x[i], x[j]);
        }
    }
    return v;
}

vector<vector<double>> squareform(vector<double>& v) {
    // Convert distance matrix from vector form to square form
    // Note: vector of vectors may not be efficient - maybe just return regular array of arrays?
    unsigned long nx = (1 + sqrt(1 + 8*v.size())) / 2; // must be an integer if v is valid
    vector<vector<double>> M (nx, vector<double>(nx, 0.0));
    unsigned long ind;
    for (int i = 1; i < nx; ++i) {
        for (int j = 0; j < i; ++j) {
            ind = nx*j - j*(j+1)/2 + i - 1 - j;
            M[i][j] = v[ind];
            M[j][i] = v[ind];
        }
    }
    return M;
}

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
    cout << levenstein(s1, s0) << endl;
    cout << levenstein(s2, s0) << endl;
    cout << levenstein(s3, s0) << endl;
    cout << levenstein(s2, s1) << endl;
    cout << levenstein(s3, s1) << endl;
    cout << levenstein(s3, s2) << endl;

    // Some more strings
    cout << "Calculating more edit distances..." << endl;
    string t1 = "g5a23g5a3av34acra2ct";
    string t2 = "bacbzV#35v#V#nk";
    string t3 = "fdafdsafsafsa";

    cout << levenstein(t1, t2) << endl;  // 17
    cout << levenstein(t2, t3) << endl;  // 15
    cout << levenstein(t1, t2, 0.9, 1.1, 2.0) << endl;

    // Test out distance matrix calc
    cout << "Calculating distance matrix in vector form..." << endl;
    vector<string> strs {s0, s1, s2, s3};
    vector<double> dv = pdist(strs, levenstein);  // CLion incorrectly says dist arg is an error
    for (double dvi : dv) {
        cout << dvi << endl;
    }

    // Test out squareform
    vector<double> v {1, 2, 3, 4, 5, 6};
    // M should be:
    // 0 1 2 3
    // 1 0 4 5
    // 2 4 0 6
    // 3 5 6 0
    auto M = squareform(v);

    auto dm = squareform(dv);

    return 0;
}