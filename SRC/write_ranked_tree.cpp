#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iomanip> // output more decimals
#include "string_manipulations.h"
#include "node.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"


using namespace std;

void writeUnrankTreeIntoStr(Node * p, int N, std::string & str)
{
    if(p -> label.size() != 0) str += p -> label;
    else
    {
        str += '(';
        writeUnrankTreeIntoStr(p -> left, N, str);
        str += ',';
        writeUnrankTreeIntoStr(p -> right, N, str);
        str += ')';
    }
}


void writeRankTreeIntoStr(Node * p, int & temp, int N, std::string & str)
{
    if(p -> label.size() != 0) str += p -> label;
    else
    {
        str += '(';
        writeRankTreeIntoStr(p -> left, temp, N, str);
        str += ',';
        writeRankTreeIntoStr(p -> right, temp, N, str);
        str += ')';
    }
    if(p -> rank != 1)
    {
        if(p -> label.size() == 0)
        {
            temp = p -> rank - p -> parent -> rank;
            str += ':';
            str +=  to_string(temp);
        }
        else
        {
            temp = p -> rank - p -> parent -> rank + N;
            str += ':';
            str +=  to_string(temp);
        }
    }
}


void treeProcessing(Node* newnode, int & N, ofstream & file)
{
    Node ** ar = new Node * [N];
    int tail = 0;
    pushToArray(newnode, tail, ar);
    getRanks(newnode, tail, ar);

    int temp = 0;
    string str="";
    writeRankTreeIntoStr(newnode, temp, N, str);
    file << str << ";" << endl;
    cout << str << ";" << endl;


    delete[] ar;
}


void writeRankedTree(int & arg_counter, char * argv[])
{
    Node * newnode;
    int N = getNumberOfTaxa(arg_counter, argv, newnode);
    ofstream tree_file("outRankTree.txt");
   // ofstream tree_file(argv[arg_counter]);
   // ++arg_counter;
    treeProcessing(newnode, N, tree_file);
    tree_file.close();
    deleteTree(newnode);
}


int calcNumberOfTaxa (std::string str)
{
    int lbl = 0;
    int i = 0;
    while(i < (int) str.size())
    {
	if(((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= 'A' && str[i] <= 'Z')) && ((i+1) < (int) str.size()))
        {
            
	    if ((static_cast<int> (str[i+1] - '0') >= 0 && static_cast<int> (str[i+1] - '0') <= 9)) ++i; //allow t1,t2,... taxa names
            if(str[i+1] == ',' || str[i+1] == ')')
            {
                ++lbl;
            }
        }
        ++i;
    }
    return lbl;
}



void writeRankedUnrTreesIntoVect(int & arg_counter, char * argv[], vector<string> & gts_vect)
{
    Node * newnode;
    ifstream tree_file(argv[arg_counter]);
    //++arg_counter;
    string str = "";
    getline(tree_file, str, '\n');
    tree_file.close();

    int N = calcNumberOfTaxa(str);
    string input_str = "";
   // ofstream out_file("res.txt");

    tree_file.open(argv[arg_counter]);
    while(getline(tree_file, input_str, '\n'))
    {
        if(input_str.size() < 3) break;
        std::stack <Node *> stkGTunr;
        int lblUnrGT = 0;
        string * ar_strlbl = new string[N-1];
        pushNodesUnrankedGT(lblUnrGT, stkGTunr , input_str, ar_strlbl);
        newnode = stkGTunr.top();
        getDescTaxa(newnode, N);
        newnode -> rank = 1;

        int prod = 1;
        str = "";
        int temp = 0;
        int tops_count = numberOfRankings(newnode, N, prod);
        vector<vector<int>> seq = permuteRanks(newnode);
        if(seq.size() != tops_count) cout << "ERROR: not all ranks!" << endl;
        for(int i = 0; i < seq.size(); ++i)
        {
            str = "";
            assignRanks(newnode, seq[i]);
            writeRankTreeIntoStr(newnode, temp, N, str);
      //    out_file << str << ";" << endl;
            gts_vect.push_back(str);
        }

        delete [] ar_strlbl;
        deleteStack(stkGTunr);
        deleteTree(newnode);

    }
 //   out_file.close();
    tree_file.close();
}





void rankUnrTrees(int & arg_counter, char * argv[])
{
    vector<string> gts_vect;
    writeRankedUnrTreesIntoVect(arg_counter, argv, gts_vect);
    ofstream out_file("outRankTopos.txt");
    for(int i = 0; i < gts_vect.size(); ++i)
    {
        out_file << gts_vect[i] << ";" << endl;
 
    }
    out_file.close();
}






