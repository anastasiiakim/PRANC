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
    BagOfNodes(vector<string> str)
    {
        for(unsigned i = 0; i < str.size(); ++i)
        {
            m.push_back(new Node(str[i]));
        }
    };
    BagOfNodes(vector<string> str, const BagOfNodes & b)
    {
        for(unsigned i = 0; i < str.size(); ++i)
        {
            for(unsigned j = 0; j < b.m.size(); ++j)
            {
                if (b.m[j]->label == str[i])
                {
                    m.push_back(b.m[j]);
                    break;
                } 
            }
        }
    };
    void clip(Node * rule)
    {
        auto it = m.end();
        it = find(m.begin(), m.end(), rule->left);
        if (it != m.end()) m.erase(it);
        it = find(m.begin(), m.end(), rule->right);
        if (it != m.end())
        {
            m.erase(it);
            m.push_back(rule);
        }
    };
};

Node* makeTree(vector<vector<string>> strs);

inline bool cmpt(const int & a, const int & b);

bool sortMapByVal(const pair<string, int> &a, const pair<string, int> &b);

bool sortVectBySize(const pair<string, int> & a, const pair<string, int> &b); //ToDo 

ret_t visit(const int & row, const int & col, const unsigned long & freq, unsigned nsize, unsigned long a[][2], vector<int> path);

void getConsensusTree(int &  arg_counter, char* argv[]);

#endif // GC_H
