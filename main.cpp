#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <string>
#include <thread>

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

template <class T>
vector<double> pdist(const vector<T>& x, double (*dist)(const T&, const T&)) {
    // Distance matrix between all pairs in x, in order of x, in vector form (upper triangle entries), using specified
    // distance function
    size_t nx = x.size();
    size_t ind;
    vector<double> v (nx*(nx-1)/2, 0.0);
    for (size_t i = 0; i < nx-1; ++i) {
        for (size_t j = i+1; j < nx; ++j) {
            ind = j - 1 - i*(3+i-2*nx)/2;
            v[ind] = dist(x[i], x[j]);
        }
    }
    return v;
}

//template <class T>
//vector<double> pdist1(const vector<T>& x, const T& q) {
//    // Distance vector between points in x and single point q
//
//}

template <typename Iterator>
vector<pair<Iterator, Iterator>> split_chunks(Iterator begin, Iterator end, size_t n) {
    // Divide iterable into multiple iterables of similar size
    // Source: https://codereview.stackexchange.com/questions/106773/dividing-a-range-into-n-sub-ranges
    vector<pair<Iterator, Iterator>> ranges;
    if (n == 0) {
        return ranges;
    }
    ranges.reserve(n);

    auto dist = std::distance(begin, end);
    n = std::min<size_t>(n, dist);
    auto chunk = dist / n;
    auto remainder = dist % n;

    for (size_t i = 0; i < n-1; ++i) {
        auto next_end = std::next(begin, chunk + (remainder ? 1 : 0));
        ranges.emplace_back(begin, next_end);

        begin = next_end;
        if (remainder) remainder -= 1;
    }

    // last chunk
    ranges.emplace_back(begin, end);
    return ranges;
}

template <class T>
void pdist_range(vector<tuple<size_t,size_t,size_t>>::iterator begin, vector<tuple<size_t,size_t,size_t>>::iterator end,
                 vector<double>& v, const vector<T> &x, double (*dist)(const T &, const T &)) {
    // Helper function to calculate pairwise distances specified by range
    for (auto it = begin; it != end; ++it) {
        v[get<2>(*it)] = dist(x[get<0>(*it)], x[get<1>(*it)]);
    }
}

template <class T>
vector<double> ppdist(const vector<T> &x, double (*dist)(const T &, const T &), int n_threads) {
    // Parallel version of pdist
    // Get indices of diagonal elements
    size_t nx = x.size();
    size_t nv = nx*(nx-1)/2;
    size_t ind = 0;
    vector<tuple<size_t,size_t,size_t>> inds;
    inds.reserve(nv);
    for (size_t i = 0; i < nx-1; ++i) {
        for (size_t j = i+1; j < nx; ++j) {
            inds.emplace_back(tuple<size_t,size_t,size_t>(i, j, ind));
            ++ind;
        }
    }

    // Split indices between threads
    auto chunks = split_chunks(inds.begin(), inds.end(), n_threads);

    // Debug: show inds assigned to each chunk
//    for (auto chunk : chunks) {
//        for (auto it = chunk.first; it != chunk.second; ++it) {
//            cout << "(" << get<0>(*it) << "," << get<1>(*it) << ") ";
//        }
//        cout << endl;
//    }

    // Run
    vector<double> v (nv, 0.0);
    for (auto chunk : chunks) {
        pdist_range(chunk.first, chunk.second, v, x, dist);
    }

    return v;
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
    vector<double> dv = pdist(strs, levenstein);  // CLion incorrectly says dist arg is an error
    for (double dvi : dv) {
        cout << dvi << endl;
    }

    cout << "Calculating distance matrix in vector form in parallel..." << endl;
    int n_threads = 4;
    vector<double> dvp = ppdist(strs, levenstein, n_threads);  // CLion incorrectly says dist arg is an error
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

    return 0;
}