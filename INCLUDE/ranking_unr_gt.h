#ifndef RU_H
#define RU_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip> // output more decimals

using namespace std;

void pushNodesUnrankedGT(int & lbl, std::stack <Node *> & stk, string & str, string * strlbl);

void searchTreeForRanks(Node * p, int & nrank, int * ranks, int & i);
void checkRanks(Node * p, int & i);

void swap (int &x, int &y);
//void writeRankTreeIntoStr(Node * p, int & temp, int N, std::string & str);

//Node * getDescTaxa(Node * p, int irank);
void topsNumber(Node * p, int N, int & prod);
string getGTtopology(Node * p, int N);
int numberOfRankings(Node * newnodeGT, int NumTaxa, int prod);


//int generateAllRanksPermutations(Node * newnode, Node * newnodeGT, int & j, int n, int N, int Numtaxa, int * ranks, string * strlbl, int ** ar_y, double * s, double & probability, int & topscounter, const int tops_count, string & temp_str, string strGT,  ofstream & file, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int *** k);

vector<vector<int> > permuteRanks(Node *p);
Node * getNodeFromName(Node * p, int name_val);
//void getNodeFromName(Node * p, Node * & r, int name);

void assignRanks(Node * p, vector<int> v);
double searchOverRanks(Node * newnode, Node * newnodeGT, int Numtaxa, int ** ar_y, double * s, ofstream & file, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int *** k);


void calcUnrankedProb(int & arg_counter, char * argv[], int & N, Node* newnode, double*s, int** ar_y);
void calcProbsUnrankedGtInput(int & arg_counter, char * argv[]);

#endif //RU

