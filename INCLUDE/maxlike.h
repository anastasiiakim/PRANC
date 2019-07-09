#ifndef ML_H
#define ML_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <functional>




void speciesTreeNoBeaded(Node* newnode, int & N, double* s_times, double * s);
void unBeadSpeciesTree(Node * p);

void writeNewickTree(Node * p, int & temp, int N, std::string & str);
void swapTwoNodes(Node * a, Node * b, Node * root, int N, vector<string> & str_vect);
void swapLabels(Node * a, Node * b, Node * root, int N, vector<string> & str_vect);
void insertBetweenNodes(Node * p, Node * root, int N, vector<string> & str_vect);
void cladeTwoNodes(Node * p, Node * ch,  Node * root, int N, vector<string> & str_vect);
void searchInternalNodes(Node * p, Node * root, int N, vector<string> & str_vect);

int pickInitialCandidate(vector<string> sp_vect, vector<string> gts_vect, int & N, Node *& newnode);
void writeTreeWithBrLengths(Node * p, int & temp, int N, std::string & str, double * s, double x);

void distFrRootBasedOnNode(Node * p, int N);
void readUnrankedTrees(string str_ST, int N, int rounds, vector<int> index_vector, int ** ar_y, int ***k, double * s, vector<string> gts_vect, double & threshold, string & candidate_str, double x);
void rankUnrankedTrees(Node * spnode, int Numtaxa, int rounds, vector<int> index_vector, int ** ar_y, int *** k, double * s, vector<string> gts_vect, double & threshold, string & candidate_str, double x);
double calcNegativeLike(int & N, Node * newnode, double* s, int** ar_y, vector<string> gts_vect);//, istream & finGT);
double fbrent (double a, double b, double t, double * s, int index, int & N, Node * newnode, int ** ar_y, vector<string> gts_vect, vector<vector<vector<double>>> & vect);
void fisher_shuffle(vector<int> & v);
void calcLikeNoNNI(int &  arg_counter, char* argv[]);
void calcLikeWithNNI(int &  arg_counter, char* argv[]);

double getGT(int & N, double * s,  Node * newnode, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k, vector<vector<double>> & vect);
void initProbsVect(int * m, int *** k, double * s, double * coal, int & n, vector<vector<double>> & vect);
double updateProbability(int & N, double * s,  Node * newnode, Node * newnodeGT, int ** ar_y, Node ** arMrca, int ** ar_rankns, int ***k, vector<vector<vector<double>>> &  vect, int & interval, int & gt_count);
double updateHistoryProbability(vector<vector<vector<double>>> & vect, int & hist_counter, int & n, int & gt_count);
double updateNegativeLike(int & N, Node * newnode, double* s, int** ar_y, int ***k, vector<string> & gts_vect, vector<vector<vector<double>>> & vect, int & interval);
double newCalcNegativeLike(int & N, Node * newnode, double* s, int** ar_y, vector<string> & gts_vect, vector<vector<vector<double>>> & vect);
#endif //  ML_H
