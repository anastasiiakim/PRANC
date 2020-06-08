#ifndef ML_H
#define ML_H
#include <string.h>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <functional>




void speciesTreeNoBeaded(Node* newnode, int & N, double* s_times, double * s, vector <Node *> v);
void unBeadSpeciesTree(Node * p);

void writeNewickTree(Node * p, int & temp, int N, std::string & str);
void swapTwoNodes(Node * a, Node * b, Node * root, int N, vector<string> & str_vect);
void swapLabels(Node * a, Node * b, Node * root, int N, vector<string> & str_vect);
void insertBetweenNodes(Node * p, Node * root, int N, vector<string> & str_vect);
void cladeTwoNodes(Node * p, Node * ch,  Node * root, int N, vector<string> & str_vect);
void searchInternalNodes(Node * p, Node * root, int N, vector<string> & str_vect);

int calcNumberOfTaxaBrL(std::string str);
string pickInitialCandidate(vector<string> sp_vect, vector<string> gts_vect, Node *& newnode, int & N, vector <Node *> v, double init_h, double init_a, double init_b, int N_subset);
void writeTreeWithBrLengths(Node * p, int & temp, int N, std::string & str, double * s, double x);

void distFrRootBasedOnNode(Node * p, int N, vector <Node *> & v);
void distFrRootBasedOnRanks(Node * p, int N, vector <Node *> & v);
void readUnrankedTrees(string str_ST, int N, int rounds, vector<int> index_vector, int ** ar_y, int ***k, double * s, vector <Node *> v, vector<string> gts_vect, double & threshold, string & candidate_str, double x, double init_h, double init_a, double init_b, int N_subset, int N_maxsubset);
void rankUnrankedTrees(Node * spnode, int Numtaxa, int rounds, vector<int> index_vector, int ** ar_y, int ***k, double * s, vector <Node *> & vlabels, vector<string> gts_vect, double & threshold, string & candidate_str, double x, double init_a, double init_b, double init_h, int N_subset, int N_maxsubset);

double calcNegativeLike(int & N, Node * newnode, double* s, vector <Node *> v, int** ar_y, vector<string> gts_vect);//, istream & finGT);
void fisher_shuffle(vector<int> & v);
void calcLikeNoNNI(int &  arg_counter, char* argv[]);
void calcLikeWithNNI(int &  arg_counter, char* argv[]);

double getGT(int & N, double * s,  vector <Node *> v, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k, vector<vector<double>> & vect);
void initProbsVect(int * m, int *** k, double * s, double * coal, int & n, vector<vector<double>> & vect);
double updateProbability(int & N, double * s,  vector <Node *> v, Node * newnodeGT, int ** ar_y, Node ** arMrca, int ** ar_rankns, int ***k, vector<vector<vector<double>>> &  vect, int & interval, int & gt_count);
double updateHistoryProbability(vector<vector<vector<double>>> & vect, int & hist_counter, int & n, int & gt_count);
double updateNegativeLike(int & N, double* s, vector <Node *> v, int** ar_y, int ***k, vector<string> & gts_vect, vector<vector<vector<double>>> & vect, int & interval);
double newCalcNegativeLike(int & N, double* s, vector <Node *> v, int** ar_y, vector<string> & gts_vect, vector<vector<vector<double>>> & gts_probs_vect);
double lbfgs(double h, double * s, int & N, vector <Node *> vlabels, int ** ar_y, vector<string> gts_vect, double init_a, double init_b);

#endif //  ML_H
