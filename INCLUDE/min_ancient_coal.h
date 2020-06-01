#ifndef AC_H
#define AC_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>

int g_i(int & N, int & i, vector <Node *> v, Node * nodeGT);
int totalNumExtraLineages(int & N, vector <Node *> v, Node * nodeGT);
int calcCostForST(int & arg_counter, char* argv[], int & N, vector <Node *> v);
void searchCandidateSpTreeTopology(int &  arg_counter, char* argv[]);


#endif //  AC_H
