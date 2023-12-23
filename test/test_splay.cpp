
#include "splay.hpp"
#include "test.hpp"

#include <map>
#include <functional>
#include <vector>

using std::cout, std::endl, std::string, std::map, std::function, std::vector;

test insert();
test remove();

int main(int ac, char** av) {
    map<string, test(*)()> tests = {
        {"insert", insert},
        {"remove", remove},
    };

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

test insert() 
{
    test t("insert");

    splay_tree<char> tree;
    tree.insert('a');

    if(tree.contains('a'))
        return t.success();

    return t.failure("'a' was not found in splay_tree"); 
}

test remove() {
    test t("remove");

    splay_tree<string> tree;
    vector<string> to_insert = {
        "mercury", "venus", "earth", "mars", "jupiter", "saturn", "uranus", "neptune", "pluto"
    };

    for(auto s : to_insert)
        tree.insert(s);

    if(!tree.contains("pluto"))
        return t.failure("tree did not contain inserted data 'pluto'");

    tree.remove("pluto");
    if(tree.contains("pluto"))
        return t.failure("removed item 'pluto' was still found in tree");

    return t.success();
}
