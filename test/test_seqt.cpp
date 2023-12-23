
#include "seqt.hpp"
#include "test.hpp"
#include <map>

#include <iostream> 

using std::ostream, std::cout, std::endl, std::cerr, std::map;

test insert()
{
    test t("insert");

    seqt m;
    

    return t.success();
}

test unordered()
{
    test t("insert");

    seqt m;

    return t.success();
}

int main(int ac, char** av) {
    map<string, test(*)()> tests = {
        {"insert", insert},
        {"unordered", unordered},
    };

    return run_tests(tests, ac, av);
}