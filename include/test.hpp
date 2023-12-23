#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <iostream>
#include <string>
#include <sstream>

using std::string, std::cout, std::endl, std::stringstream;

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