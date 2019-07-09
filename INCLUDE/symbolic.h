#ifndef SY_H
#define SY_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>

double symbolicIntervalProbability(int i, int * m_i, int *** array_k, double * s, ofstream & fout);
double symbolicGeneTreeProbability(int * m, int *** k, double * s, double * coal, int n, ofstream & fout);
double getOneSymbolicGeneTreeProb(int N, double * s,  Node * newnode, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k, ofstream & fout, ofstream & histprobs_file);
double getSymbolicGeneTreeProb(int & arg_counter, char* argv[], int & N, Node * newnode, double* s, int** ar_y);
void symbolicProbsRankedGtInput(int &  arg_counter, char* argv[]);


#endif //  SY_H

