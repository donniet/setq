
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
    seqt::node n = m.create_atom('n');
    n = m.splay_insert(n, 'm');
    n = m.splay_insert(n, 'l');
    if(!m.find_in_unordered(n, 'n')) {
        return t.failure("could not find inserted symbol in unorderd node.");
    }
    if(m.find_in_unordered(n, 'q')) {
        return t.failure("found non-inserted symbol in unordered node.");
    }
    if(!m.find_in_unordered(n, 'm')) {
        return t.failure("could not find inserted symbol in unorderd node.");
    }
    if(!m.find_in_unordered(n, 'l')) {
        return t.failure("could not find inserted symbol in unorderd node.");
    }

    seqt::node n1 = m.splay_remove(n, 'm');
    if(m.find_in_unordered(n1, 'm'))
        return t.failure("found a removed symbol in unordered node");

    if(!m.find_in_unordered(n, 'm'))
        return t.failure("could not find symbol in original node after we splay_removed");
    
    return t.success();
}

int main(int ac, char** av) {
    map<string, test(*)()> tests = {
        // {"insert", insert},
        {"unordered", unordered},
    };

    return run_tests(tests, ac, av);
}