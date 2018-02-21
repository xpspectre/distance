#ifndef DISTANCE_DISTANCE_H
#define DISTANCE_DISTANCE_H

#include <vector>
#include <tuple>
#include <string>
#include <thread>
#include <iostream>

double euclidean(const std::vector<double>& x1, const std::vector<double>& x2);

double levenstein(const std::string& s1, const std::string& s2, double del, double ins, double sub);
double levenstein(const std::string& s1, const std::string& s2);

std::vector<std::vector<double>> squareform(std::vector<double>& v);

void print_mat(std::vector<std::vector<double>>& M);

template <typename T>
std::vector<double> pdist(const std::vector<T>& x, double (*dist)(const T&, const T&)) {
    // Distance matrix between all pairs in x, in order of x, in vector form (upper triangle entries), using specified
    // distance function
    size_t nx = x.size();
    size_t ind;
    std::vector<double> v (nx*(nx-1)/2, 0.0);
    for (size_t i = 0; i < nx-1; ++i) {
        for (size_t j = i+1; j < nx; ++j) {
            ind = j - 1 - i*(3+i-2*nx)/2;
            v[ind] = dist(x[i], x[j]);
        }
    }
    return v;
}

template <typename Iterator>
std::vector<std::pair<Iterator, Iterator>> split_chunks(Iterator begin, Iterator end, size_t n) {
    // Divide iterable into multiple iterables of similar size
    // Source: https://codereview.stackexchange.com/questions/106773/dividing-a-range-into-n-sub-ranges
    std::vector<std::pair<Iterator, Iterator>> ranges;
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

template <typename T>
void pdist_range(std::vector<std::tuple<std::size_t,std::size_t,std::size_t>>::const_iterator begin,
                 std::vector<std::tuple<std::size_t,std::size_t,std::size_t>>::const_iterator end,
                 std::vector<double>& v, const std::vector<T>& x, double (*dist)(const T&, const T&)) {
    // Helper function to calculate pairwise distances specified by range
    for (auto it = begin; it != end; ++it) {
        v[std::get<2>(*it)] = dist(x[std::get<0>(*it)], x[std::get<1>(*it)]);
    }
}

template <typename T>
std::vector<double> ppdist(const std::vector<T> &x, double (*dist)(const T &, const T &), int n_threads) {
    // Parallel version of pdist
    // Get indices of diagonal elements
    size_t nx = x.size();
    size_t nv = nx * (nx - 1) / 2;
    size_t ind = 0;
    std::vector <std::tuple<size_t, size_t, size_t>> inds;
    inds.reserve(nv);
    for (size_t i = 0; i < nx - 1; ++i) {
        for (size_t j = i + 1; j < nx; ++j) {
            inds.emplace_back(i, j, ind);
            ++ind;
        }
    }

    // Split indices between threads
    auto chunks = split_chunks(inds.begin(), inds.end(), n_threads);
    // Debug: show inds assigned to each chunk
//    for (auto chunk : chunks) {
//        for (auto it = chunk.first; it != chunk.second; ++it) {
//            std::cout << "(" << std::get<0>(*it) << "," << std::get<1>(*it) << ") ";
//        }
//        std::cout << std::endl;
//    }

    // Run
    std::vector<double> v(nv, 0.0);
    // Single-threaded
//    for (auto chunk : chunks) {
//        pdist_range(chunk.first, chunk.second, v, x, dist);
//    }
    // Multi-threaded
    std::vector<std::thread> ths;
    for (auto chunk : chunks) {
        ths.emplace_back(pdist_range<T>, chunk.first, chunk.second, std::ref(v), std::cref(x), dist);  // CLion incorrectly reports wrong # args
    }
    for (auto &th: ths) {
        th.join();
    }
    return v;
}

#endif //DISTANCE_DISTANCE_H
