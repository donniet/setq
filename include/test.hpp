#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <iostream>
#include <string>
#include <sstream>

using std::string, std::cout, std::endl, std::stringstream;

template<typename TESTS>
int run_tests(TESTS & tests, int ac, char * av[]) {
    bool all_succeeded = true;

    auto run_all = [&]() {
        for(auto j = tests.begin(); j != tests.end(); j++)
            if(!j->second())
                all_succeeded = false;
    };

    for(int i = 1; i < ac; i++) {
        if(string("all") == av[i]) {
            run_all();
            continue;
        }

        auto j = tests.find(av[i]);
        if(j == tests.end())
            continue;

        if(!j->second())
            all_succeeded = false;
    }

    if(ac <= 1)
        run_all();

    if(all_succeeded) {
        cout << "ALL TESTS SUCCEEDED" << endl;
        return 0;
    }
    
    cout << "TEST FAILED" << endl;
    return -1;
}

struct test {
    bool successful;
    string name;
    string message;
    bool moved;

    test(string n)
        : successful(true), name(n), moved(false)
    { }

    test(test & r) 
        : successful(r.successful), name(r.name), message(r.message),
          moved(false)
    {
        r.moved = true;
    }

    test(test && r) 
        : successful(r.successful), name(r.name), message(r.message),
          moved(false)
    {
        r.moved = true;
    }

    operator bool() const {
        return successful;
    }

    ~test()
    {
        if(moved) return;

        if(successful) 
            cout << "TEST \"" << name << "\" PASSED" << endl;
        else
            cout << "TEST \"" << name << "\" FAILED\n"
                 << " \t" << message << endl;
    }

    test& success() 
    {
        successful = true;
        return *this;
    }

    test& failure(string msg)
    {
        successful = false;
        message = msg;
        return *this;
    }

};


#endif