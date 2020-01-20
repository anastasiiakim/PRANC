#ifndef WT_H
#define WT_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip> // output more decimals

using namespace std;

void writeRankTreeIntoStr(Node * p, int & temp, int N, std::string & str);
void writeUnrankTreeIntoStr(Node * p, int N, std::string & str);
void treeProcessing(Node* newnode, int & N, ofstream & file);
void writeRankedUnrTreesIntoVect(int & arg_counter, char * argv[], vector<string> & gts_vect);
void writeRankedTree(int & arg_counter, char * argv[]);
void rankUnrTrees(int & arg_counter, char * argv[]);
int calcNumberOfTaxa(string str);

#endif //WT
