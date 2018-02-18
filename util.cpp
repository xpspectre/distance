#include "util.h"
#include <random>
#include <algorithm>

using std::string;
using std::vector;
using std::size_t;

string random_str(size_t n) {
    // Generate random string
    // TODO: also make it have a random length, more efficient initialization of the RNG
    const static string chars = "abcdefghijklmnaoqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890";
    std::mt19937_64 gen { std::random_device()() };
    std::uniform_int_distribution<size_t> dist { 0, chars.length()-1 };
    std::string str;
    std::generate_n(std::back_inserter(str), n, [&] { return chars[dist(gen)]; });
    return str;
}

vector<string> random_strs(size_t N, size_t n) {
    // Generate vector of random strings
    vector<string> strs;
    strs.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        strs.emplace_back(random_str(n));
    }
    return strs;
}
