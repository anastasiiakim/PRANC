#ifndef GC_H
#define GC_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cmath>
#include <algorithm>
#include <sstream>

typedef pair<unsigned long, vector<int>> ret_t;

string idToStr(int id);

int strToId(string str);

struct BagOfNodes
{
    vector<Node*> m;
    // Constructor
    BagOfNodes(){};
    BagOfNodes(vector<string> str);
    BagOfNodes(vector<string> str, const BagOfNodes & b);
    void clip(Node *);
};

Node* makeTree(vector<vector<string>> strs);

inline bool cmpt(const int & a, const int & b);

bool sortMapByVal(const pair<string, int> &a, const pair<string, int> &b);

bool sortVectBySize(const pair<string, int> & a, const pair<string, int> &b); //ToDo 

ret_t visit(const int & row, const int & col, const unsigned long & freq, unsigned nsize, unsigned long a[][2], vector<int> path);

void getConsensusTree(int &  arg_counter, char* argv[]);

#endif // GC_H
