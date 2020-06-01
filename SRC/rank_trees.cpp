#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip> // output more decimals
#define eps 0.0e-10
#include "string_manipulations.h"
#include "node.h"
#include "queue.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"
#include "write_ranked_tree.h"
#include "rank_trees.h"

using namespace std;


void searchWriteRanks(Node * newnodeGT, int numTaxa, ofstream & file)
{
    int temp = 0;
    string out_str = "";
    vector<vector<int>> seq = permuteRanks(newnodeGT);
    for(int i = 0; i < seq.size(); ++i)
    {
	out_str = "";
	assignRanks(newnodeGT, seq[i]);
	writeRankTreeIntoStr(newnodeGT, temp, numTaxa, out_str);
	file << out_str << ";" << endl;
     }        
}


void rankUnrankedTrees(int & arg_counter, char * argv[])
{ 
    string strGT="";
    ifstream finGT(argv[arg_counter]); 
    ++arg_counter;
    ofstream outUnr(argv[arg_counter], ios::out | ios::app);
    Node * newnodeGT;
    
    int N = 0;
    int lblUnrGT = 0;

    while(getline(finGT, strGT, '\n'))
    {
        if(strGT.size() < 3) break;

        std::stack <Node *> stkGTunr;
        lblUnrGT = 0;
        pushNodesUnrankedGT(lblUnrGT, stkGTunr , strGT);
        N = lblUnrGT;
        newnodeGT = stkGTunr.top();
        getDescTaxa(newnodeGT, N);
        newnodeGT -> rank = 1;
        searchWriteRanks(newnodeGT, N, outUnr);
       
        deleteStack(stkGTunr);
        deleteTree(newnodeGT);
    }
 
    finGT.close();
    outUnr.close();
} 

