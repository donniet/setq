#include "partition.hpp"

#include <vector>
#include <tuple>

/* 
number of partitions of N identical items into K distinct bins P(n,k)

Recurrence Relation: 

         0 in bin0  1 in bin0    2 in bin0    3 in bin0          n in bin0
P(n,k) = P(n,k-1) + P(n-1,k-1) + P(n-2,k-1) + P(n-3,k-1) + ... + P(0,k-1)
       = P(n,k-1) + P(n-1,k)  <-- pascal!

        P(0,k) = 1   one way to add zero balls to k bins
for n>0 P(n,0) = 0   zero ways to add a non-zero number of balls to zero bins
        
ex:
P(1,k) = 

*/

unsigned long partition_identical_n_distinct_k(unsigned long n, unsigned long k) 
{
    using namespace std;

    static vector<unsigned long> cache;
    static auto get_cache = [](unsigned long n, unsigned long k, vector<unsigned long> & cache) -> unsigned long &
    {
        unsigned long mi = min(n,k);
        unsigned long ma = max(n,k);

        unsigned int diagonal = n + k + 1;                  // zero indexed
        unsigned int start = diagonal * (diagonal + 1) / 2; // sum of all elements in the triangle before this diagonal
        unsigned int next_start = (diagonal + 1) * (diagonal + 2) / 2;
        unsigned int index = start + k;                     // kth row

        if(cache.size() < next_start) {
            cache.resize(next_start, 0);
        }

        return cache[index];
    };

    unsigned long & cached_value = get_cache(n, k, cache);
    if(cached_value != 0 || (k == 0 && n > 0))
        return cached_value;

    vector<tuple<unsigned long, unsigned long, bool>> stack;
    stack.push_back({n, k, false});

    bool returning = false;
    while(!stack.empty()) {
        tie(n, k, returning) = stack.back();
        stack.pop_back();

        if(returning) {
            // values should all be cached at this point
            unsigned long val = 0;
            for(unsigned long i = 0; i <= n; i++) {
                val += get_cache(i, k-1, cache);
            }
            get_cache(n, k, cache) = val;
            continue;
        }

        if(k == 0)
            get_cache(n, k, cache) = 1;
        else if(n == 0)
            get_cache(n, k, cache) = 0;
        else if(k == 1)
            get_cache(n, k, cache) = 1;
        else {
            // need to recurse
            stack.push_back({n, k, true});          // add a step for us to recacluate from cache
            for(unsigned long i = 0; i <= n; i++) {
                stack.push_back({i, k-1, false});   // recurse down the pyramid
            }
        }
    }

    // value should be cached by this point
    return get_cache(n, k, cache);
}