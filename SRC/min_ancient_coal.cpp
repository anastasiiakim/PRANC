#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <iomanip> // output more decimals
#define EPS 1e-8
#define MAX_EXPONENT_VALUE 1024
#include "string_manipulations.h"
#include "node.h"
#include "queue.h"
#include "probs_calculation.h"


using namespace std;




int g_i(int & N, int & i, vector <Node *> v, Node * nodeGT)
{
   std::stack <Node *> allTaxa;
   int res = 0;
   int prod;
   Node * mrcanode;
   for(int j = i + 1; j < N; ++j)
   {
    prod = 1;
    for(int k = j; k < N; ++k)
    {
        getTaxa(getNodeFromRank(nodeGT, k), allTaxa);
        mrcanode  = getMrca(v, allTaxa);

        if(mrcanode -> rank <= i) prod = 0;
    }
        res += prod;
   }

   return N - res;
}

int totalNumExtraLineages(int & N, vector <Node *> v, Node * nodeGT)
{

    int sum = 0;
    for(int t = 1; t < N; ++t)
    {
      //  cout << t << " " << g_i(N, t, nodeST, nodeGT) << endl;
        sum += g_i(N, t, v, nodeGT) - (t + 1);
    }
    return sum;
}



int calcCostForST(int & arg_counter, char* argv[], int & N, vector <Node *> v)
{
    string strGT="";
    ifstream finGT(argv[arg_counter]); //gtuniqtrees.txt
    ++arg_counter;


    Node * newnodeGT;
    int cost = 0;

    while(getline(finGT, strGT, ';'))
    {
        if(strGT.size() < 3) break;
        removeSpaces(strGT);


        std::stack <Node *> stkGT;
        int lblGT = 0;
        Node **arGT = new Node * [N];
        pushNodes(lblGT, stkGT , strGT);
        newnodeGT = stkGT.top();
        newnodeGT -> distFrRoot = 0.;
        int tailGT = 0;
        distFromRoot(newnodeGT);
        pushToArray(newnodeGT, tailGT, arGT);
        getRanks(newnodeGT, tailGT, arGT);
        getDescTaxa(newnodeGT, N);
        
//        cout << totalNumExtraLineages(N, newnode, newnodeGT) << endl;
        cost += totalNumExtraLineages(N, v, newnodeGT);

        delete [] arGT;

        deleteTree(newnodeGT);
        deleteStack(stkGT);
    }

    finGT.close();
    return cost;

}


void searchCandidateSpTreeTopology(int &  arg_counter, char* argv[])
{
    ofstream fout("outMacScore.txt");
    Node * newnode;
    int N = getNumberOfTaxa(arg_counter, argv, newnode);
    double * s_times = new double [N-1];
    double * s = new double [N-2];
    vector <Node *> v;
    int ** ar_y = new int * [N];
    for (int i = 0; i < N; i++) ar_y[i] = new int [N];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            ar_y[i][j] = 0;
        }

    speciesTreeProcessing(newnode, N, s_times, s, v, ar_y);

    fout << calcCostForST(arg_counter, argv, N, v) << endl;

    for(int i = 0; i < N; ++i)
    {
        delete[] ar_y[i];
    }
    delete[] ar_y;
    delete[] s_times;
    delete[] s;
    deleteTree(newnode);
    fout.close();
}

