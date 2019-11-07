#ifndef UT_H
#define UT_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>


void getClades (Node * p, std::stack <string> & allTaxa);
void finSort(int n, string* str);
void getUnrDescTaxa(Node * p, int N);
int countNumberTaxa(int & indicator, string & str);
void outputUnrankedTopology(int & arg_counter, char* argv[]);

#endif //  UT_H
