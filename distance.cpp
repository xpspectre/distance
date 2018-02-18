#include "distance.h"
#include <algorithm>
#include <cmath>
#include <random>

using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::pair;
using std::tuple;
using std::get;
using std::size_t;

double euclidean(const vector<double>& x1, const vector<double>& x2) {
    double d2 = 0.0;
    for (size_t i = 0; i < x1.size(); ++i) {
        d2 += std::pow(x1[i] - x2[i], 2);
    }
    return std::sqrt(d2);
}

double levenstein(const string& s1, const string& s2, double del, double ins, double sub) {
    // Edit distance of strings using user-specified costs
    // Note that the del and ins costs should be the same to ensure the distance is symmetric
    // Note: doesn't use default params since function pointers passed as args into other functions don't recognize them
    // Source: http://www.python-course.eu/levenshtein_distance.php

    // For all i and j, dist[i,j] will contain the Levenshtein distance between the first i characters of a and the
    // first j characters of b
    size_t n_rows = s1.size() + 1;
    size_t n_cols = s2.size() + 1;
    vector<vector<double>> dist (n_rows, vector<double>(n_cols, 0.0)); // do this properly; it can get big

    // Source prefixes can be transformed into empty strings by deletions
    for (size_t row = 1; row < n_rows; ++row) {
        dist[row][0] = row * del;
    }

    // Target prefixes can be created from an empty source string by inserting the characters
    for (size_t col = 1; col < n_cols; ++col) {
        dist[0][col] = col * ins;
    }

    double cost;
    vector<double> subcost {0.0, 0.0, 0.0};
    for (size_t col = 1; col < n_cols; ++col) {
        for (size_t row = 1; row < n_rows; ++row) {
            if (s1[row-1] == s2[col-1]) {
                cost = 0.0;
            } else {
                cost = sub;
            }
            subcost[0] = dist[row-1][col] + del;
            subcost[1] = dist[row][col-1] + ins;
            subcost[2] = dist[row-1][col-1] + cost;  // substitution
            dist[row][col] = *std::min_element(subcost.begin(), subcost.end());
        }
    }

    return dist[n_rows-1][n_cols-1];
}

double levenstein(const string& s1, const string& s2) {
    // Edit distance of strings using default costs
    static const double del = 1.0, ins = 1.0, sub = 1.0;
    return levenstein(s1, s2, del, ins, sub);
}

vector<vector<double>> squareform(vector<double>& v) {
    // Convert distance matrix from vector form to square form
    // Upper triangular indexes ind = (nx choose 2) - (nx-1 choose 2) + j - i - 1
    // See: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.distance.squareform.html#scipy.spatial.distance.squareform
    // Note: vector of vectors may not be efficient - maybe just return regular array of arrays?
    size_t nx = (1 + std::sqrt(1 + 8*v.size())) / 2; // must be an integer if v is valid
    vector<vector<double>> M (nx, vector<double>(nx, 0.0));
    unsigned long ind;
    for (size_t i = 0; i < nx-1; ++i) {
        for (size_t j = i+1; j < nx; ++j) {
            ind = j - 1 - i*(3+i-2*nx)/2;
            M[i][j] = v[ind];
            M[j][i] = v[ind];
        }
    }
    return M;
}

void print_mat(vector<vector<double>>& M) {
    // Print matrix in convenient form
    for (vector<double> row : M) {
        for (double e : row) {
            cout << e << " ";
        }
        cout << endl;
    }
}
