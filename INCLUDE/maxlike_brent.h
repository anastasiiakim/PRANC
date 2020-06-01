#ifndef MLB_H
#define MLB_H
#include <string.h>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <functional>

double fbrent (double a, double b, double * s, int index, int & N, vector <Node *> vlabels, int ** ar_y, int ***k, vector<string> gts_vect, vector<vector<vector<double>>> & vect);
void brent_rankUnrankedTrees(Node * spnode, int Numtaxa, int rounds, vector<int> index_vector, int ** ar_y, int ***k, double * s, vector <Node *> & vlabels, vector<string> gts_vect, double & threshold, string & candidate_str, double x);

string brent_pickInitialCandidate(vector<string> sp_vect, vector<string> gts_vect, Node *& newnode, int & N, vector <Node *> v);
void brent_calcLikeWithNNI(int &  arg_counter, char* argv[]);

#endif //  MLB_H
