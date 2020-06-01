#ifndef BD_H
#define BD_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cmath>
#include <algorithm>

void outputCoalIntervals(int &  arg_counter, char* argv[]);
void writeBranchLengths(unordered_map <string, double> & clade_and_rank, string str, int & N);
void getMSEBranchDistance(int &  arg_counter, char* argv[]);

#endif // BD_H
