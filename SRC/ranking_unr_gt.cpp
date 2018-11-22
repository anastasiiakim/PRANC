#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector>
#include <iomanip> // to print more decimals
#include "string_manipulations.h"
#include "node.h"
#include "queue.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"


using namespace std;

void pushNodesUnrankedGT(int & lbl, std::stack <Node *> & stk, string & str, int N, string * strlbl)
{
    int i = 0;
    int j = 0;
    while(i < (int) str.size())
    {
        if((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= 'A' && str[i] <= 'Z'))
        {
            std::string strtemp;
            strtemp += str[i];    
            while(((i+1) < (int) str.size()) && str[i+1] != '(' && str[i+1] != ',' && str[i+1] != ')')
            {
                strtemp += str[i+1];		   
                ++i;
            }
            Node * newnode = new Node();
            newnode -> label = strtemp;
            ++lbl;	

            stk.push(newnode);
        }

        else if (str[i] == ')')
        {	
            Node * newnode = new Node();
            newnode -> right = stk.top();
            newnode -> right -> parent = newnode;
            newnode -> right -> isRight = 1;
            stk.pop();
            newnode -> left = stk.top();
            newnode -> left -> parent = newnode;
            newnode -> left -> isRight = 0;
            stk.pop();

            stack <Node *> allTaxa;
            getTaxa(newnode, allTaxa);
            string temp;
            while(!allTaxa.empty())
            {
                temp += allTaxa.top() -> label;
                allTaxa.pop();
            }
            strlbl[j] = temp;
            ++j;
            stk.push(newnode);
        } 
        ++i;
    }

}	

void searchTreeForRanks(Node * p, int & nrank, int * ranks, int & i)
{
    if(p -> label.size() == 0)
    {
        if(p -> rank != 1)
        {
            if((p -> parent -> desctaxa).find(p -> desctaxa) != string::npos)
            {
                p -> rank = ranks[nrank];
                ++nrank;
                if(p -> parent -> rank < p -> rank) ++i;
                else return;  
            }
        }
        searchTreeForRanks(p -> left, nrank, ranks, i);
        searchTreeForRanks(p -> right, nrank, ranks, i);
    }
}

void checkRanks(Node * p, int & i)
{
    if(p -> label.size() == 0)
    {
        if(p -> rank != 1)
        {
            if(p -> parent -> rank < p -> rank)  
            {
                ++i;
            } 
        }
        checkRanks(p -> right, i);
        checkRanks(p -> left, i);
    }
}

void swap (int &x, int &y)
{
    int temp = x;
    x = y;
    y = temp;
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


void topsNumber(Node * p, int N, int & prod)
{
    if(p -> label.size() == 0)
    {
        prod *= (p -> leavesnum - 1);
        topsNumber(p -> right, N, prod);   
        topsNumber(p -> left, N, prod);   
    }
    else return;
}

string getGTtopology(Node * p, int N)
{
    string res_str="";
    for(int i = 2; i < N; ++i)
    {   
        stack <Node *> allTaxa;
        string tempstr[N];
        int tmpcount = 0;
        getTaxa(getNodeFromRank(p, i), allTaxa);
        while(!allTaxa.empty())
        {
            tempstr[tmpcount] = allTaxa.top() -> label;
            allTaxa.pop();
            tmpcount++;   
        } 
        sortString(N, tempstr);
        for(int j = 0; j < N; ++j)
        {
            res_str += tempstr[j];
        }
        res_str += "-";
        res_str += to_string(i);
        res_str += "-";
    }
    return res_str;
}

int numberOfRankings(Node * newnodeGT, int NumTaxa, int prod)
{
    topsNumber(newnodeGT, NumTaxa, prod);
    return factorial(NumTaxa-1)/prod; 
}
/*
int generateAllRanksPermutations(Node * newnode, Node * newnodeGT, int & j, int n, int N, int Numtaxa, int * ranks, string * strlbl, int ** ar_y, double * s, double & probability, int & topscounter, const int tops_count, string & temp_str, string strGT,  ofstream & file, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int *** k)
{
    int ret = 2;
    if(n == 1) 
    {    
        int icounter = 0;
        int dummy = 0;

        searchTreeForRanks(newnodeGT, dummy, ranks, icounter);
        if (icounter == N)
        {
            temp_str="";
            temp_str = getGTtopology(newnodeGT, Numtaxa);
            ++topscounter;
            double probabilitytmp = getGeneTreeProb(Numtaxa, s, newnode, newnodeGT, ar_y, array_invcoal, arMrca, ar_rankns, k);
            probability += probabilitytmp;

            file << temp_str << " " << setprecision(15) << probabilitytmp << endl;
            ++j;
            for(int irank = 2; irank < Numtaxa; ++irank)
            {
                Node * tempnode = getNodeFromRank(newnodeGT, irank);
            }

            if (topscounter > tops_count) return 3;

        }   
    }
    else
    {
        for(int i = 0; i < n - 1; ++i)
        {
            ret = generateAllRanksPermutations(newnode, newnodeGT, j, n-1, N, Numtaxa, ranks, strlbl, ar_y, s, probability,  topscounter, tops_count, temp_str, strGT, file, array_invcoal, arMrca, ar_rankns, k);
            if (ret == 3) return 3;
            if(n % 2 == 0)
                swap(ranks[i], ranks[n-1]);
            else
                swap(ranks[0], ranks[n-1]);
        }
        ret = generateAllRanksPermutations(newnode, newnodeGT, j, n-1, N, Numtaxa, ranks, strlbl, ar_y, s, probability,  topscounter, tops_count, temp_str, strGT, file, array_invcoal, arMrca, ar_rankns, k);
        if (ret == 3) return 3;
    }

    return 2;
}
*/
vector<vector<int> > permuteRanks(Node * p)
{
    if (p == NULL) {
        vector<int> rank_sequence;
        vector<vector<int> > v;
        v.push_back(rank_sequence);
        return v;
    }

    vector<vector<int> > results, left_part, right_part;
    left_part  = permuteRanks(p -> left);
    right_part = permuteRanks(p -> right);
    int size = right_part[0].size() + left_part[0].size() + 1;

    vector<bool> flags(left_part[0].size(), 0);
    for (int i = 0; i < right_part[0].size(); ++i)
        flags.push_back(1);

    for (int i = 0; i < left_part.size(); ++i) {
        for (int j = 0; j < right_part.size(); ++j) {
            do {
                vector<int> vtemp(size);
                vtemp[0] = p -> debug_name;
                if(p -> label.size() != 0) vtemp.pop_back();
                int l = 0, r = 0;

                for (int t = 0; t < flags.size(); ++t) {
                    if(flags[t])
                    {
                        vtemp[t+1] = right_part[j][r++]; 
                    }
                    else
                    {
                        vtemp[t+1] = left_part[i][l++];
                    }
                }
                results.push_back(vtemp);
          
            } while (next_permutation(flags.begin(), flags.end()));
        }
    }
    return results;
}


Node * getNodeFromName (Node * p, int rankValueGT)
{
    if(p == NULL) cout << "Tree is empty" << endl;
    Node * k;
    Node * t;
    if(p -> label.size() == 0)
    {       
        if(p -> debug_name != rankValueGT) 
        {
            t =  getNodeFromName(p -> right, rankValueGT);
            k =  getNodeFromName(p -> left, rankValueGT);
            if(t != NULL) return t;
            return k;
        }
        else 
        {
            return p;
        }   
    }
    else return 0;
}   

void assignRanks(Node * newnodeGT, vector<int> v)
{
        Node * p;
        for(int j = 1; j < v.size(); ++j) 
        {   
            int ret = 0;
            p = getNodeFromName(newnodeGT, v[j]);
            p -> rank = j+1;
        }
}


double searchOverRanks(Node * newnode, Node * newnodeGT, int Numtaxa, int ** ar_y, double * s,  ofstream & file, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int *** k)
{
    double probability = 0.;
    int prod = 1;
    int tops_count = numberOfRankings(newnodeGT, Numtaxa, prod);
    //cout << "Number of topologies: " << tops_count << endl;
    vector<vector<int>> rank_seq = permuteRanks(newnodeGT);
    if(rank_seq.size() != tops_count) cout << "ERROR: not all ranks!" << endl;

    for(int i = 0; i < rank_seq.size(); ++i)
    {
        assignRanks(newnodeGT, rank_seq[i]);
        double probabilitytmp = getGeneTreeProb(Numtaxa, s, newnode, newnodeGT, ar_y, array_invcoal, arMrca, ar_rankns, k);
        probability += probabilitytmp;

        file << getGTtopology(newnodeGT, Numtaxa) << " " << setprecision(15) << probabilitytmp << endl;
    }        
    return probability;    
}


void calcUnrankedProb(int & arg_counter, char * argv[], int & N, Node * newnode, double * s, int ** ar_y)
{ 
    string strGT="";
    ifstream finGT(argv[arg_counter]); //"gtuniqtrees.txt"
    ++arg_counter;
    ofstream outGT("outUnrGT.txt", ios::out | ios::app);
    double TOTAL_PROBABILITY = 0.;

    int prod;
    int count;
    int topcounter;
    string temp_str="";
   
    Node ** arMrca = new Node * [N-1];
    int ** ar_rankns = new int * [N-1];
    for (int i = 0; i < N-1; ++i) ar_rankns[i] = new int [N-1];
    for(int i = 0; i < N-1; ++i)
    {
        for(int j = 0; j < N - 1; ++j)
        {
            ar_rankns[i][j] = 0; 
        }
    }	

    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; ++i)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }

    double *  array_invcoal = new double[N-1];
    getInvPartialCoal(array_invcoal, N);
    Node * newnodeGT;
    
    int * ranks = new int[N-2];
    while(getline(finGT, strGT, '\n'))
    {
        if(strGT.size() < 3) break;
        count = 0;
        topcounter = 1;
        prod = 1;
        double oneGTPr = 0.;

        std::stack <Node *> stkGTunr;
        int lblUnrGT = 0;
        int M = N;//lbl;
        string * ar_strlbl = new string[M-1];
        pushNodesUnrankedGT(lblUnrGT, stkGTunr , strGT, M, ar_strlbl);

        newnodeGT = stkGTunr.top();
        getDescTaxa(newnodeGT, N);
        newnodeGT -> rank = 1;
    
        oneGTPr = searchOverRanks(newnode, newnodeGT, N, ar_y, s, outGT, array_invcoal, arMrca, ar_rankns, k);

        TOTAL_PROBABILITY += oneGTPr;
        delete [] ar_strlbl;
        deleteStack(stkGTunr);
        deleteTree(newnodeGT);

    }
    
    delete [] ranks;
    delete [] array_invcoal;
    
    for(int i = 0; i < N-1; ++i)
        delete[] ar_rankns[i];
    delete [] ar_rankns;

    delete[] arMrca;
    
    for(int i = 0; i < N+1; ++i)
    {
        for(int j = 0; j < N; ++j)
            delete [] k[i][j];
        delete [] k[i];
    }
    delete [] k;

    
    finGT.close();
    outGT.close();
   
} 

void calcProbsUnrankedGtInput(int & arg_counter, char * argv[])
{
    Node * newnode;
    int N = getNumberOfTaxa(arg_counter, argv, newnode);
    double * s = new double [N-1];
    int ** ar_y = new int * [N];
    for (int i = 0; i < N; ++i) ar_y[i] = new int [N];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            ar_y[i][j] = 0;
        }

    speciesTreeProcessing(newnode, N, s, ar_y);
    calcUnrankedProb(arg_counter, argv, N, newnode, s, ar_y);

    for(int i = 0; i < N; ++i)
    {
        delete[] ar_y[i];	
    }
    delete[] ar_y;	
    delete[] s;
    deleteTree(newnode);
}
