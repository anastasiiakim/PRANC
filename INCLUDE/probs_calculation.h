#ifndef PM_H
#define PM_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <node.h>

void pushToArray (Node * p, int & tail, Node ** ar);
void pushToArrayNoSort (Node * p, int & tail, Node ** ar);
Node * popFromArray (int & tail, Node ** ar);


void distFromRoot (Node * p);
void isUltrametric (Node * p, double & tempdist);
void saveDistFrRoot(Node * p, double * ar, int & i);
void printIsUltrametric (double * ar, int N);

void deleteTree (Node * p);
void deleteStack (std::stack <Node *> stk);

void pushNodes (int & lbl, std::stack <Node *> & stk, std::string str);
void getRanks(Node* newnode, int & tail, Node ** ar);
Node * getNodeFromRank (Node * p, int rankValueGT);
void getTaxa (Node * p, std::stack <Node *> & allTaxa);

Node * mrcaST (Node * rch, Node * lch);
Node * mrcaST2 (Node * p, Node * tp, std::string v);
Node * getMrca(vector <Node *> vlabels, std::stack <Node * > & stk);

void writeRankedTreeTopology(Node * p, int N, ofstream & file);
void storeLabels (Node * p, int N, string * ar, int & i);

void makeBeadedTree (Node * root, size_t MaxQ);
void makeBeadedArray (Node * root, int ** y, size_t MaxQ);

void getPath (Node * p, std::stack <Node *> & stk);
void getPathMatchRank(Node * &  p, int rank_val);

int * getNextHistory (int * ar, int * ar_orig, int n, int idx, int & numb);
void arrayFreq(int * arr, int * res, int n);
void get_hseq(int * arseq, Node * nodeGT, vector <Node *> v, int n);
void reverseArraySort(int * arseq, int n);

void getDescTaxa(Node * node, int N);
void returnMrca(vector <Node *> v, Node * nodeGT, int n, Node ** ar);
void get_node(int * rank_hist, int ** temp, int n, Node ** ar);

void calc_k_ijz (int ***k, int n, int * m, int ** y, int ** ar_r);

void writeTree(Node* root);

double factorial (int num);
double calcBinomCoef (int n, int k);

void storeLabels(int n, int & i, string* str, Node * p);
void getS(Node * p, double * s, vector <Node *> & v);

double lambda_ij (int i, int j, int *** k);
double invPartialCoal(int j);
double conditionalProbability(int i, int * m_i, int *** array_k, double * s); // min i = 2 since l[i-2] exists
void getInvPartialCoal(double * ar, int N);
double geneTreeProbability(int * m, int *** k, double * s, vector <Node *> vlabels, double * coal, int n);

double getGeneTreeProb(int N, double * s, vector <Node *> vlabels, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k);



int getNumberOfTaxa(int & arg_counter, char* argv[], Node* &newnode);
void speciesTreeProcessing(Node* newnode, int & N, double* s_times, double* s, vector <Node *> & vlabels, int** ar_y);
double calcRankedProb(int & arg_counter, char* argv[], int & N, double* s, vector <Node *> vlabels, int** ar_y);
void calcProbsRankedGtInput(int &  arg_counter, char* argv[]);
void outputCoalIntervals(int &  arg_counter, char* argv[]);


#endif //  PM_H
