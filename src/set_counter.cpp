#include <iostream>
#include <iomanip>

#include <random>
#include <string>
#include <algorithm>

#include <cmath>

#include <boost/math/special_functions.hpp>


/*
probability that given n trials of p that x or more trials were successful
TODO: could use the incomplete beta function, but since we are dealing with very large values and
    the binomial distribution converges on the gaussian this should be asymptoically accurate
*/
float probability_follows(float n, float p, float x) {
    float var = n * p * (1. - p);
    float exp = n * p;

    return 0.5 * erfc((x - exp) * sqrt(2 / var));
}



int main(int ac, char ** av) {
    using std::cerr, std::cout;
    using namespace std;

    // make random sequences with 3 symbols (a, b, x), and plot the distribution of C(xa) / C(a) and C(xb) / C(b)
    auto rd = mt19937 {random_device{}()};

    // TODO: replace with string from av
    string sym = "abx";

    float replace_ab_ratio = 0;

    cerr << "RUNNING: " << av[0];
    if(ac > 1) { sym = av[1]; cerr << " " << av[1]; }
    if(ac > 2) { rd.seed(atoi(av[2])); cerr << " " << av[2]; }
    if(ac > 3) { 
        replace_ab_ratio = atof(av[3]);
        if(replace_ab_ratio < 0)
            replace_ab_ratio = 0;
        else if(replace_ab_ratio > 1)
            replace_ab_ratio = 1;

        cerr << " " << replace_ab_ratio; 
    }

    if(sym.size() < 2) {
        cerr << "ERROR: must have more than one symbol" << endl;
        return -1;
    }

    auto ab_coin = uniform_real_distribution(0., 1.);

    vector<int> cnt, rem, exp;
    vector<vector<int>> sequence_counts;

    cnt.resize(sym.size(), 0);
    rem.resize(sym.size(), 0);
    sequence_counts.resize(sym.size(), cnt);
    
    int char_index = 0,
        last_char_index = -1,
        sequence_length = 1e6;

    // pick a random proportion of a:b:x
    for(int i = 0, remaining = sequence_length; i < sym.size(); i++) {
        if(i == sym.size() - 1) {
            cnt[i] = remaining;
        } else {
            cnt[i] = rd() % remaining;
            remaining -= cnt[i];
        }
    }

    // counter of remaining symbols
    rem = cnt;
    exp = cnt;
    // reset our counts
    fill(cnt.begin(), cnt.end(), 0);

    cerr << endl << setw(7) << "a^b" << ",";

    for(int i = 0; i < sym.size(); i++)
        cerr << setw(7) << sym[i] << ",";

    for(int i = 0; i < sym.size(); i++) 
        for(int j = 0; j < sym.size(); j++) {
            cerr << setw(6) << sym[i] << sym[j] << ",";
        }

    cout << setw(7) << "E[a]," << setw(7) << "E[b],";
    cout << setw(7) << "total" << endl;

    cout << setw(7) << replace_ab_ratio << ",";

    for(int i = 0; i < sequence_length; i++) {

        // get one element at random
        char_index = 0;
        
        for(int r = rd() % (sequence_length - i); r > rem[char_index];)
            r -= rem[char_index++];


        // swap the first two symbols on the roll of a die
        // this simulates how two symbols can be replaced with one another because
        // they have a similar meaning in the context
        // we are going to want to back-into this `replace_ab_ratio` number from the data
        // always flip the coin though, so the random number generator stays synced
        auto rep_ab = ab_coin(rd);

        if( (char_index == 0 && rep_ab < replace_ab_ratio)
         || (char_index == 1 && rep_ab < replace_ab_ratio) ) 
        {
            char_index = (char_index + 1) % 2;
        }
        else {
            // cerr << "not replacing" << endl;
        }


        // one fewer symbols remaining
        rem[char_index]--;
        

        // increment our actual symbols after replacing
        cnt[char_index]++;
        if(last_char_index >= 0) {
            sequence_counts[last_char_index][char_index]++;
        }

        last_char_index = char_index;
    }

    for(int i = 0; i < sym.size(); i++) 
        cout << setw(7) << cnt[i] << ",";
    
    for(int i = 0; i < sym.size(); i++)
        for(int j = 0; j < sym.size(); j++)
            cout << setw(7) << sequence_counts[i][j] << ",";

    cout << setw(7) << floor((double)exp[0] * (1.-replace_ab_ratio) + (double)exp[1] * replace_ab_ratio) << ",";
    cout << setw(7) << floor((double)exp[1] * (1.-replace_ab_ratio) + (double)exp[0] * replace_ab_ratio) << ",";
        
    cout << sequence_length << endl;


    if(ac > 4 && string(av[4]) == "dump") {
        cout << setw(7) << "E[_]" << ",";
        
        // dump out expected values
        for(int i = 0; i < sym.size(); i++)
            cout << setw(7) << exp[i] << ",";

        for(int i = 0; i < sym.size(); i++) {
            for(int j = 0; j < sym.size(); j++)
                cout << setw(7) << round((float)cnt[i] * (float)cnt[j] / (float)sequence_length) << ",";
        }
        cout << sequence_length << "\n";

        cout << setw(7) << "P[_]" << ",";
        for(int i = 0; i < sym.size(); i++)
            cout << setw(7) << fixed << setprecision(5) << (float)cnt[i]/(float)sequence_length << ",";

        for(int i = 0; i < sym.size(); i++) {
            for(int j = 0; j < sym.size(); j++) {
                float a = max(cnt[i], cnt[j]);
                float b = min(cnt[i], cnt[j]);
                float x = (float)sequence_counts[i][j];

                cout << setw(7) << fixed << setprecision(5) 
                     << (1.0 - probability_follows(a, b / (float)sequence_length, x)) 
                     << ",";
            }
        }

        cout << endl;
        
    }

    return 0;
}