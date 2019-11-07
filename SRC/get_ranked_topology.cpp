#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include "string_manipulations.h"
#include "probs_calculation.h"
#include "get_ranked_topology.h"

using namespace std;


Node * nodeFrRank (Node * p, int rankValueGT)
{
    if(p == NULL) cout << "Tree is empty" << endl;
    Node * k;
    Node * t;
    if(p -> label.size() == 0)
    {   

        if(p -> rank != rankValueGT) 
        {
            t =  nodeFrRank(p -> right, rankValueGT);
            k =  nodeFrRank(p -> left, rankValueGT);
            if(t != 0) return t;
            return k;
        }
        else 
        {
            return p;
        }   
    }
    else return 0;
}   


void outputRankedTopology(int & arg_counter, char* argv[])
{ 
    string str;
    string hash = "|";
    ifstream gt_file(argv[arg_counter]);
    ofstream topo_file("rtopos.txt");

    while(getline(gt_file, str, ';'))
    {
        if(str.size() < 3) break;
        removeSpaces(str);
        std::stack <Node *> stkST;

        int lbl = 0;
        double tempdist = 0.;	
        int labelcount = 0;
        countParentheses(str); 

        pushNodes(lbl, stkST , str);
        Node * newnode;
        newnode = stkST.top();
        newnode -> distFrRoot = 0.;
        int tail = 0;

        distFromRoot(newnode);

        int N = lbl;
        Node ** ar = new Node * [N];
        int ** ar_y = new int * [N];
        for (int i = 0; i < N; i++) ar_y[i] = new int [N];

        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
            {
                ar_y[i][j] = 0;
            }

        pushToArray(newnode, tail, ar);	
        getRanks(newnode, tail, ar);

        int lblscounter = 0;
        string strLabels[N];

        for(int i = 2; i < N; ++i)
        {   
            stack <Node *> allTaxa;
            string tempstr[N];
            int tmpcount = 0;
            getTaxa(nodeFrRank(newnode, i), allTaxa);
            while(!allTaxa.empty())
            {
                tempstr[tmpcount] = allTaxa.top() -> label;
                allTaxa.pop();
                tmpcount++;  
            } 
            sortString(N, tempstr);
            for(int j = 0; j < N; ++j)
            {
                if(tempstr[j].size() != 0)
                {
                    topo_file << tempstr[j] << hash;
                }
            }
            topo_file << "-" << i << "-";
        }
        topo_file << "\n";

        for(int i = 0; i < N; ++i)
        {
            delete[] ar_y[i];	
        }
        delete[] ar_y;	
        delete[] ar;

        deleteStack(stkST);
        deleteTree(newnode);
    }
    ofstream freq_file("rfreqs.txt");
    topo_file.close();

    ifstream infile;
    infile.open("rtopos.txt");
    if (infile.is_open()) countUniqueStrings(infile, freq_file);
    freq_file.close();
    gt_file.close();
    infile.close();

} 

