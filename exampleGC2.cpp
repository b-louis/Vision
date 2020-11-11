#include <iostream>

#include "maxflow/graph.h"

using namespace std;

// Use the library to compute a minimum cut on the following graph:
//
//        SOURCE
//       /       \
//     1/         \6
//     /      4    \
//   node0 -----> node1
//     |   <-----   |
//     |      3     |
//     \            /
//     5\          /1
//       \        /
//          SINK
void testGCuts() {
    const int INF = 1000000;
    const int ZERO = 0;
    Graph<int,int,int> g(2, 4); // estimated # of nodes/edges
    g.add_node(4); 
    g.add_tweights( 0,   /* capacities */  5, ZERO );
    g.add_tweights( 1,   /* capacities */  6, ZERO );

    g.add_tweights( 2,   /* capacities */  ZERO, 10 );
    g.add_tweights( 3,   /* capacities */  ZERO, 1 );

    g.add_edge( 0, 1,    /* capacities */  4, 3 );
    g.add_edge( 2, 3,    /* capacities */  0, 0 );

    g.add_edge( 0, 2,    /* capacities */  4, INF );
    g.add_edge( 1, 3,    /* capacities */  4, INF );

    int flow = g.maxflow();
    cout << "Flow = " << flow << endl;
    for (int i=0;i<4;i++)
        if (g.what_segment(i) == Graph<int,int,int>::SOURCE)
            cout << i << " is in the SOURCE set" << endl;
        else
            cout << i << " is in the SINK set" << endl;
}

int main() {
    testGCuts();
    return 0;
}
