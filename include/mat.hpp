#ifndef __MAT_HPP__
#define __MAT_HPP__

#include <vector>
#include <set>
#include <map>
#include <stdint.h>

using std::vector;
using std::set;
using std::map;

struct network {
    vector<float> charge; // charge in each neuron
    map<size_t, uint32_t> atoms;

    
    struct iterator;

    iterator atom(uint32_t);


    // creates a new neuron that receives charge from the first and
    // sends to the second
    /*

    IN     q
            \.
             x
              \.    
    OUT        u'

        x = join(q, u);

     */
    iterator join(iterator from, iterator to);

    // creates a new neuron that receives charge equally 
    // from the given neurons
    /*       
                z
             /'   '\
            x       y
          /' '\   /' '\
    IN   a     e i     o

        x = merge(a, e)
        y = merge(i, o)
        z = merge(x, y)
     */
    iterator merge(iterator a, iterator b);

    // creates a new neuron that moves charge from 
    /*

    IN    1     0    0    
          0     1    0
          0     0    1

          a     b    c
           \.    \.   \.
            x ->  y -> z     
             \.         \.    \.   
    OUT       b'         c'    d'   
            
        x = join(a, b');
        y = merge(x, b);


         


     */
    iterator branch(iterator a, iterator b);
};

#endif