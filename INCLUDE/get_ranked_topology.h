#ifndef RT_H
#define RT_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>

Node * nodeFrRank (Node * p, int rankValueGT);
void outputRankedTopology(int & arg_counter, char* argv[]);

#endif //  RT_H
