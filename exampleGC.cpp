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

//////////////////////

//   SOURCE
//     |
//    2|
//     |
//   node0
//     |  
//   1 | 
//     | 
//   node1
//     |  
//    5| 
//     |
//    SINK

//////////////////////

//   SOURCE
//     |     \
//    2|        \ 1
//     |           \
//   node0<--2-->node2
//     |            |  
//   1 |          2 |       
//     |            | 
//   node1<--2-->node3
//     |        /
//    10|    / 1
//     | /
//    SINK

void testGCuts() {
    int inf = 10000;
    Graph<int,int,int> g(2, 1); // estimated # of nodes/edges
    g.add_node(4); 
/*
    Source vers noeuds
*/
    g.add_tweights( 0,   /* capacities */  3, 0 );
    g.add_tweights( 2,   /* capacities */  4, 0 );
//
/*
    |0|--1-->|1| et |2|--2-->|3|
*/
    g.add_edge( 0, 1,    /* capacities */  3, inf );
    g.add_edge( 2, 3,    /* capacities */  4, inf );
//
/*
    |0|<--2-->|2| et |1|<--2-->|3|
*/
    g.add_edge( 0, 2,    /* capacities */  3, 3 );
    g.add_edge( 1, 3,    /* capacities */  4, 4 );
//
/*
    noeuds vers Puis
*/
    g.add_tweights( 1,   /* capacities */  0, 10 );
    g.add_tweights( 3,   /* capacities */  0, 1 );

    // g.add_node(2); 
    // g.add_tweights( 0,   /* capacities */  10, 0 );
    // g.add_edge( 0, 1,    /* capacities */  2, inf );
    // g.add_tweights( 1,   /* capacities */  0, 1 );

    // g.add_node(2); 
    // g.add_tweights( 0,   /* capacities */  1, 5 );
    // g.add_tweights( 1,   /* capacities */  6, 1 );
    // g.add_edge( 0, 1,    /* capacities */  4, 1 );
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
