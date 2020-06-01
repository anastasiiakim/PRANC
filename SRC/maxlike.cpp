#include <string.h>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <iomanip> // output more decimals
#include <functional>
#include <vector>
#include <ctime>

#include "string_manipulations.h"
#include "node.h"
#include "queue.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"
#include "write_ranked_tree.h"
#include "maxlike.h"
#include "lbfgsb.h"

#define MAX_EXPONENT_VALUE 1024
#define init_a 0.001
#define init_b 6
#define init_h 1e-10
#define N_rounds 0
#define N_nni 5
using namespace std;


double calcNegativeLike(int & N, double* s, vector <Node *> v, int** ar_y, vector<string> gts_vect)
{

    string strGT="";
    Node ** arMrca = new Node * [N-1];
    int ** ar_rankns = new int * [N-1];
    for (int i = 0; i < N-1; i++) ar_rankns[i] = new int [N-1];
    for(int i = 0; i < N-1; ++i)
    {
        for(int j = 0; j < N - 1; ++j)
        {
            ar_rankns[i][j] = 0;
        }
    }

    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; i++)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }

    double *  array_invcoal = new double[N-1];
    getInvPartialCoal(array_invcoal, N);
    Node * newnodeGT;
    double prec_total_prob = 0.;
    double one_gt_prob = 0.;
    for(int count = 0; count < gts_vect.size(); ++count)
    {
        strGT = gts_vect[count];
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
        one_gt_prob = getGeneTreeProb(N, s, v, newnodeGT, ar_y, array_invcoal, arMrca, ar_rankns, k);
        prec_total_prob += log(one_gt_prob);

        delete [] arGT;

        deleteTree(newnodeGT);
        newnodeGT = NULL;
        deleteStack(stkGT);

    }

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
    return -prec_total_prob;
}

void fisher_shuffle(vector<int> & v)
{
    int N = v.size();
    srand(time(NULL));
    for(int i = N - 1; i > 0; --i)
    {
        int pick = rand() % (i+1); // rand number from 0 to i
        swap(v[i], v[pick]);
    }

}


double newCalcNegativeLike(int & N, double* s, vector <Node *> v, int** ar_y, vector<string> & gts_vect, vector<vector<vector<double>>> & gts_probs_vect)
{

    string strGT="";
    Node ** arMrca = new Node * [N-1];
    int ** ar_rankns = new int * [N-1];
    for (int i = 0; i < N-1; i++) ar_rankns[i] = new int [N-1];
    for(int i = 0; i < N-1; ++i)
    {
        for(int j = 0; j < N - 1; ++j)
        {
            ar_rankns[i][j] = 0;
        }
    }

    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; i++)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }

    double *  array_invcoal = new double[N-1];
    getInvPartialCoal(array_invcoal, N);
    Node * newnodeGT;
    double prec_total_prob = 0.;
    double total_prob = 0.;
    double one_gt_prob = 0.;
    for(int count = 0; count < gts_vect.size(); ++count)
    {
        //cout << "---------------------------------------------------countgt = " << count << endl;
        strGT = gts_vect[count];
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
        vector<vector<double>> vect;
        one_gt_prob = getGT(N, s, v, newnodeGT, ar_y, array_invcoal, arMrca, ar_rankns, k, vect);
        gts_probs_vect.push_back(vect);
        prec_total_prob += log(one_gt_prob);
        total_prob += one_gt_prob;

        delete [] arGT;

        deleteTree(newnodeGT);
        deleteStack(stkGT);
    }


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
    return -prec_total_prob;
}



void writeNewickTree(Node * p, int & temp, int N, std::string & str)
{
    if(p -> label.size() != 0) str += p -> label;
    else
    {
        str += '(';
        writeNewickTree(p -> left, temp, N, str);
        str += ',';
        writeNewickTree(p -> right, temp, N, str);
        str += ')';
    }
}

void swapTwoNodes(Node * a, Node * b, Node * root, int N, vector<string> & str_vect)
{
    Node * tmp;
    if(a -> parent -> left -> debug_name == a -> debug_name)
    {
        a -> parent -> left = b;
        if(b -> debug_name == b -> parent -> right -> debug_name)
        {
            b -> parent -> right = a;
        }
        else if(b -> debug_name == b -> parent -> left -> debug_name)
        {
            b -> parent -> left = a;
        }
    }
    else if(a -> parent -> right -> debug_name == a -> debug_name)
    {
        a -> parent -> right = b;
        if(b -> debug_name == b -> parent -> right -> debug_name)
            b -> parent -> right = a;
        else if(b -> debug_name == b -> parent -> left -> debug_name)
            b -> parent -> left = a;

    }
    tmp = b -> parent;
    b -> parent = a -> parent;
    a -> parent = tmp;


    std::string temp_str = "";
    int temp = 0;
    writeNewickTree(root, temp, N, temp_str);
    temp_str += ";";
    str_vect.push_back(temp_str);

    if(a -> parent -> left -> debug_name == a -> debug_name)
    {
        a -> parent -> left = b;
        if(b -> debug_name == b -> parent -> right -> debug_name)
        {
            b -> parent -> right = a;
        }
        else if(b -> debug_name == b -> parent -> left -> debug_name)
        {
            b -> parent -> left = a;
        }
    }
    else if(a -> parent -> right -> debug_name == a -> debug_name)
    {
        a -> parent -> right = b;
        if(b -> debug_name == b -> parent -> right -> debug_name)
            b -> parent -> right = a;
        else if(b -> debug_name == b -> parent -> left -> debug_name)
            b -> parent -> left = a;

    }
    tmp = b -> parent;
    b -> parent = a -> parent;
    a -> parent = tmp;


}

void swapLabels(Node * a, Node * b, Node * root, int N, vector<string> & str_vect)
{
    std::string temp_str = "";
    int temp = 0;
    std::string tmp = "";


    tmp = b -> label;
    b -> label = a -> label;
    a -> label = tmp;

    writeNewickTree(root, temp, N, temp_str);
    temp_str += ";";
    str_vect.push_back(temp_str);
    tmp = b -> label;
    b -> label = a -> label;
    a -> label = tmp;

}


void insertBetweenNodes(Node * p, Node * root, int N, vector<string> & str_vect)
{
    Node * tmp;
    std::string temp_str = "";
    int temp = 0;

    if (p -> parent -> left -> debug_name == p -> debug_name)
    {
        tmp = p -> parent -> right;
        p -> left -> parent = p -> parent;
        p -> parent -> right -> parent = p;
        p -> parent -> left = p -> left;
        p -> parent -> right = p;
        p -> left = tmp;
        writeNewickTree(root, temp, N, temp_str);
        temp_str += ";";
        str_vect.push_back(temp_str);

        //to original tree
        p -> parent -> right = p -> left;
        p -> left -> parent = p -> parent;
        p -> left = p -> parent -> left;
        p -> parent -> left -> parent = p;
        p -> parent -> left = p;

        tmp = p -> parent -> right;
        p -> right -> parent = p -> parent;
        p -> parent -> right -> parent = p;
        p -> parent -> left = p -> right;
        p -> parent -> right = p;
        p -> right = tmp;

        temp_str="";
        writeNewickTree(root, temp, N, temp_str);
        temp_str += ";";
        str_vect.push_back(temp_str);
        //to original tree
        p -> parent -> right = p -> right;
        p -> right -> parent = p -> parent;
        p -> right = p -> parent -> left;
        p -> parent -> left -> parent = p;
        p -> parent -> left = p;

    }
    else
    {
        tmp = p -> parent -> left;
        p -> right -> parent = p -> parent;
        p -> parent -> left -> parent = p;
        p -> parent -> right = p -> right;
        p -> parent -> left = p;
        p -> right = tmp;

        writeNewickTree(root, temp, N, temp_str);
        temp_str += ";";
        str_vect.push_back(temp_str);
        //to original tree
        p -> parent -> left = p -> right;
        p -> right -> parent = p -> parent;
        p -> right = p -> parent -> right;
        p -> parent -> right -> parent = p;
        p -> parent -> right = p;



        tmp = p -> parent -> left;
        p -> left -> parent = p -> parent;
        p -> parent -> left -> parent = p;
        p -> parent -> right = p;
        p -> parent -> left = p -> left;
        p -> left = tmp;


        temp_str="";
        writeNewickTree(root, temp, N, temp_str);
        temp_str += ";";
        str_vect.push_back(temp_str);
        //to original tree
        tmp = p -> parent -> left;
        p -> parent -> left = p -> left;
        p -> left -> parent = p -> parent;
        p -> left = tmp;
        tmp -> parent = p;
        p -> parent -> right = p;


    }

}


void cladeTwoNodes(Node * p, Node * ch,  Node * root, int N, vector<string> & str_vect)
{
    Node * tmp;
    int idx = 0;
    if (p -> left -> label.size() == 0)
    {
        tmp = p -> left;
        tmp -> parent = p -> parent;
        p -> left = ch;
        ch -> parent = p;
        p -> parent -> right = p;
        p -> parent -> left = tmp;

    }
    else
    {
        tmp = p -> right;
        tmp -> parent = p -> parent;
        p -> right = ch;
        ch -> parent = p;
        p -> parent -> right = p;
        p -> parent -> left = tmp;
        idx = 1;

    }
    std::string temp_str = "";
    int temp = 0;
    writeNewickTree(root, temp, N, temp_str);
    temp_str += ";";
    str_vect.push_back(temp_str);

    if (idx == 0)
    {
        tmp = p -> parent -> left;
        if(p -> parent -> left -> debug_name == p -> debug_name)
        {
            p -> parent -> right = ch;
        }
        else
        {
            p -> parent -> left = ch;
        }
        ch -> parent = p -> parent;
        p -> left = tmp;
        tmp -> parent = p;
    }
    else
    {
        tmp = p -> parent -> left;
        if(p -> parent -> left -> debug_name == p -> debug_name)
        {
            p -> parent -> right = ch;
        }
        else
        {
            p -> parent -> left = ch;
        }
        ch -> parent = p -> parent;
        p -> right = tmp;
        tmp -> parent = p;
    }



}

void searchInternalNodes(Node * p, Node * root, int N, vector<string> & str_vect)
{
    if (p -> label.size() == 0)
    {
        if(p -> rank != 1)
        {

            if(p -> parent -> left -> debug_name == p -> debug_name)
            {

                if(p -> parent -> right -> label.size() !=0)
                {
                    if (p -> left -> label.size() != 0 && p -> right -> label.size() != 0)
                    {
                        swapLabels(p -> left, p -> parent -> right, root, N, str_vect);
                        swapLabels(p -> right, p -> parent -> right, root, N, str_vect);
                    }
                    else if (p -> left -> label.size() != 0)
                    {
                        swapLabels(p -> left, p -> parent -> right, root, N, str_vect);
                        cladeTwoNodes(p, p -> parent -> right, root, N, str_vect);

                    }
                    else if (p -> right -> label.size() != 0)
                    {
                        swapLabels(p -> right, p -> parent -> right, root, N, str_vect);
                        cladeTwoNodes(p, p -> parent -> right, root, N, str_vect);
                    }
                    else
                    {
                        insertBetweenNodes(p, root, N, str_vect);

                    }

                }
                else
                {
                    swapTwoNodes(p, p -> parent -> right -> left, root, N, str_vect);
                    swapTwoNodes(p, p -> parent -> right -> right, root, N, str_vect);


                }
            }
            else if(p -> parent -> right -> debug_name == p -> debug_name)
            {
                if(p -> parent -> left -> label.size() != 0)
                {
                    if (p -> left -> label.size() != 0 && p -> right -> label.size() != 0)
                    {
                        swapLabels(p -> left, p -> parent -> left, root, N, str_vect);
                        swapLabels(p -> right, p -> parent -> left, root, N, str_vect);

                    }
                    else if (p -> left -> label.size() != 0)
                    {
                        swapLabels(p -> left, p -> parent -> left, root, N, str_vect);
                        cladeTwoNodes(p, p -> parent -> left, root, N, str_vect);
                    }
                    else if (p -> right -> label.size() != 0)
                    {
                        swapLabels(p -> right, p -> parent -> left, root, N, str_vect);
                        cladeTwoNodes(p, p -> parent -> left, root, N, str_vect);



                    }
                    else
                    {
                        insertBetweenNodes(p, root, N, str_vect);

                    }

                }
                else
                {

                    swapTwoNodes(p, p -> parent -> left -> left, root, N, str_vect);
                    swapTwoNodes(p, p -> parent -> left -> right, root, N, str_vect);


                }
            }
        }
        searchInternalNodes(p -> right, root, N, str_vect);
        searchInternalNodes(p -> left, root, N, str_vect);
    }
}

void distFrRootBasedOnRanks(Node * p, int N, vector <Node *> & v)
{
    if(p -> label.size() == 0)
    {
        p -> distFrRoot = p -> rank - 1;
        p -> right -> original_parent = p;
        p -> left -> original_parent = p;
        distFrRootBasedOnRanks(p -> left, N, v);
        distFrRootBasedOnRanks(p -> right, N, v);



    }
    else
    {
        p -> distFrRoot = N - 1;
        v.push_back(p);
    }
}


void rankUnrankedTrees(Node * spnode, int Numtaxa, int rounds, vector<int> index_vector, int ** ar_y, int ***k, double * s, vector <Node *> & vlabels, vector<string> gts_vect, double & threshold, string & candidate_str, double x)
{
    int prod = 1;

    size_t maxQue = (size_t) Numtaxa*(1+Numtaxa)/2;

    int tops_count = numberOfRankings(spnode, Numtaxa, prod);
    vector<vector<int>> seq = permuteRanks(spnode);
    double s_array[seq.size()][Numtaxa-2];
    double res[seq.size()];

    if(seq.size() != tops_count) cout << "ERROR: not all ranks!" << endl;
    //    int NR = 2*Numtaxa; // how many rankings to select

 //   cout << "overal  ranks: " << seq.size() << endl;
    int sample_NR = 2*Numtaxa; // how many rankings to select
   // int sample_NR = 2*Numtaxa; // how many rankings to select
    //int check_NR = (int) (Numtaxa/2.); // how many rankings to select
    int check_NR = (int) (Numtaxa); // how many rankings to select
   // int check_NR = 2*Numtaxa; // how many rankings to select
    int diff_NR = sample_NR - check_NR; // how many rankings to select
    double neg_loglike = 0.;
    int th_counter = 0;
//    cout << "threshold: " << threshold << endl;


    vector<int> shuffle_rankings;  
    for(int i = 0; i < seq.size(); ++i)
        shuffle_rankings.push_back(i+1);
    fisher_shuffle(shuffle_rankings);
    //select subset of ranks
    int NR;
    if(seq.size() > sample_NR)
    {
        NR = diff_NR;
    }
    else
    {
        NR = seq.size() - check_NR;
        if(check_NR >= seq.size()) 
        {
            check_NR = seq.size();
            NR = check_NR;
        }
    }  

  //  cout << "**************************** NR: " << NR << endl;
    for(int i = 0; i < check_NR; ++i)
    {
        if(vlabels.size() != 0)
        {
            while (!vlabels.empty())
            {
                vlabels.pop_back();
            }
        }

        assignRanks(spnode, seq[shuffle_rankings[i]-1]);
        distFrRootBasedOnRanks(spnode, Numtaxa, vlabels);

        makeBeadedTree(spnode, maxQue);
        makeBeadedArray(spnode, ar_y, maxQue);

        //   vector<vector<vector<double>>> vect;
        //   neg_loglike = newCalcNegativeLike(Numtaxa, s, vlabels,  ar_y, gts_vect, vect);

        double min_temp = 0.;
        neg_loglike = lbfgs(init_h, s, Numtaxa, vlabels, ar_y, gts_vect);
    //    cout << i << " " << neg_loglike << endl;

        unBeadSpeciesTree(spnode);

        res[i] = neg_loglike;
        for(int t = 0; t < Numtaxa-2; ++t)
            s_array[i][t] = s[t];
        if(neg_loglike < threshold)
        {
            th_counter++;
        }
    }

    double min_res = res[0];
//    cout << "----- min_res = " << min_res << " " << th_counter  << " " << check_NR<< endl; 
    if(th_counter != 0)
    {

        if(NR == check_NR)
        {
            int j = 0;
            for(int i = 1; i < NR; ++i)
            {
                if(res[i] < min_res)
                {
                    min_res =  res[i];
                    j = i;
                }
            }
   //         cout << "+++++best ranking j: " << j << endl;
  //          cout << "+++++min_res = " << min_res << " threshold = " << threshold << endl;

            if(min_res < threshold)
            {

                if(vlabels.size() != 0)
                {
                    while (!vlabels.empty())
                    {
                        vlabels.pop_back();
                    }
                }

                threshold = min_res;
 //               cout << "NR fin_res = " << threshold << endl;
                assignRanks(spnode, seq[shuffle_rankings[j]-1]);
                distFrRootBasedOnRanks(spnode, Numtaxa, vlabels);
                int temp = 0;
                candidate_str = "";
                for(int t = 0; t < Numtaxa - 2; ++t)
                    s[t] = s_array[j][t];
                writeTreeWithBrLengths(spnode, temp, Numtaxa, candidate_str, s, x);
            }
        }

        else
        {    
            for(int i = check_NR; i < NR; ++i)
            {
                // cout << "ranking ----------------------------- " << i << endl;
                if(vlabels.size() != 0)
                {
                    while (!vlabels.empty())
                    {
                        vlabels.pop_back();
                    }
                }

                assignRanks(spnode, seq[shuffle_rankings[i]-1]);
                distFrRootBasedOnRanks(spnode, Numtaxa, vlabels);

                makeBeadedTree(spnode, maxQue);
                makeBeadedArray(spnode, ar_y, maxQue);

                //   vector<vector<vector<double>>> vect;
                //   neg_loglike = newCalcNegativeLike(Numtaxa, s, vlabels,  ar_y, gts_vect, vect);

                double min_temp = 0.;
                neg_loglike = lbfgs(init_h, s, Numtaxa, vlabels, ar_y, gts_vect);
     //           cout << i << " " << neg_loglike << endl;

                unBeadSpeciesTree(spnode);
                res[i] = neg_loglike;
                for(int t = 0; t < Numtaxa-2; ++t)
                    s_array[i][t] = s[t];

            }




            int j = 0;
            for(int i = 1; i < NR; ++i)
            {
                if(res[i] < min_res)
                {
                    min_res =  res[i];
                    j = i;
                }
            }
    //        cout << "-best ranking j: " << j << endl;
   //         cout << "-min_res = " << min_res << " threshold = " << threshold << endl;

            if(min_res < threshold)
            {
                if(vlabels.size() != 0)
                {
                    while (!vlabels.empty())
                    {
                        vlabels.pop_back();
                    }
                }
                threshold = min_res;
                assignRanks(spnode, seq[shuffle_rankings[j]-1]);
                distFrRootBasedOnRanks(spnode, Numtaxa, vlabels);
                int temp = 0;
                candidate_str = "";
                for(int t = 0; t < Numtaxa - 2; ++t)
                    s[t] = s_array[j][t];
                writeTreeWithBrLengths(spnode, temp, Numtaxa, candidate_str, s, x);
            }

        }
    }
}


void calcLikeNoNNI(int & arg_counter, char* argv[])
{
    //	double a = 0.0001;
    //	double b = 6.;
    //	double t = 1e-10;
    int temp = 0;
    string candidate_str="";
    double x = 0.1;

    Node * newnode;
    int N = getNumberOfTaxa(arg_counter, argv, newnode);

    double * s_times = new double [N-1];
    double * s = new double [N-2];
    int ** ar_y = new int * [N];
    vector <Node *> v;
    for (int i = 0; i < N; i++) ar_y[i] = new int [N];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            ar_y[i][j] = 0;
        }

    speciesTreeProcessing(newnode, N, s_times, s, v, ar_y);
    vector<string> gts_vect;

    if(strcmp(argv[arg_counter],"-ugt") == 0)
    {
        ++arg_counter;
        writeRankedUnrTreesIntoVect(arg_counter, argv, gts_vect);
    }
    else if(strcmp(argv[arg_counter],"-rgt") == 0)
    {
        ++arg_counter;
        ifstream finGT(argv[arg_counter]);
        string strR = "";

        while(getline(finGT >> std::ws, strR, ';'))
        {
            if(strR.size() < 3) break;
            removeSpaces(strR);
            gts_vect.push_back(strR);
        }
        finGT.close();
    }


    vector<int> index_vector;
    for(int i = 0; i < N - 2; ++i)
        index_vector.push_back(i);


    double initial_br_length[N];
    for(int i = 0; i < N - 2; ++i)
        initial_br_length[i] = s[i];

    double res = 0.;

    vector<vector<vector<double>>> gts_probs_vect;
    double neg_loglike = newCalcNegativeLike(N, s, v, ar_y, gts_vect, gts_probs_vect);

    cout << "Negative log-likelihood  = " << neg_loglike << endl;
    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; i++)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }






    /*

       for(int index = N-3; index >= 0; --index)
       {
       for(int j = 0; j < 1; ++ j)
       {
       res = fbrent(a, b, t, s, index, N, newnode, ar_y, k, gts_vect, gts_probs_vect);
       for(int i = 0; i < N - 2; ++i)
       cout << s[i] << " ";
       cout << endl;
       }

       cout << "---------- " << index << " res " << res << endl;
       }
       */
    //        for(int i = 0; i < N - 2; ++i) cout << s[i] << " ";
    //	    cout << endl;



    int rounds = N_rounds;
    //int rounds = (int) (N/2.);
    //   for(int j = 0; j < rounds; ++j)
    //   {
    //     fisher_shuffle(index_vector);


    //   for(int index = 0; index < N - 2; ++index)
    // {
    res = lbfgs(init_h, s, N, v, ar_y, gts_vect);
    // cout << res << " " << res << endl;


    //		cout << j << " " << index_vector[index] << " res " << res << endl;
    //  }

    //     for(int i = 0; i < N - 2; ++i) cout << s[i] << " ";
    //   cout << endl;
    //  }


    double mse_res = 0.;
    for(int i = 0; i < N - 2; ++i)
        mse_res += fabs((s[i] - initial_br_length[i]) * (s[i] - initial_br_length[i]));

    cout << "mse: " << sqrt(mse_res) << endl;

    cout << "initial interval lengths" << endl;
    for(int i = 0; i < N - 2; ++i) cout << initial_br_length[i] << " ";
    cout << endl;

    cout << "estimated interval lengths" << endl;
    for(int i = 0; i < N - 2; ++i) cout << s[i] << " ";
    cout << endl;

    cout << "abs difference in interval lengths" << endl;
    for(int i = 0; i < N - 2; ++i) cout << fabs(s[i] - initial_br_length[i]) << " ";
    cout << endl;



    unBeadSpeciesTree(newnode);

    writeTreeWithBrLengths(newnode, temp, N, candidate_str, s, x);
    ofstream ftop("outNoNniMLTopo.txt", ios::out | ios::app);
    ftop << candidate_str << ";" << endl;
    ftop.close();

    ofstream ftopi("outNoNniMse.txt", ios::out | ios::app);
    ftopi << sqrt(mse_res) << endl;
    ftopi.close();

    ofstream ftopj("outNoNniDiffIntLengths.txt", ios::out | ios::app);
    for(int i = 0; i < N - 2; ++i) 
        ftopj << fabs(s[i] - initial_br_length[i]) << " ";
    ftopj << endl;
    ftopj.close();



    for(int i = 0; i < N+1; ++i)
    {
        for(int j = 0; j < N; ++j)
            delete [] k[i][j];
        delete [] k[i];
    }
    delete [] k;
    for(int i = 0; i < N; ++i)
    {
        delete[] ar_y[i];
    }
    delete[] ar_y;
    delete[] s_times;
    delete[] s;
    deleteTree(newnode);

}

void unBeadSpeciesTree(Node * p)
{
    if(p -> label.size() == 0)
    {

        if(p -> left != NULL)
        {

            unBeadSpeciesTree(p -> left);
            Node * temp = p -> left;
            if(p -> left -> left == NULL && p->left->right !=NULL)
            {
                p -> left -> right -> parent = p -> left -> parent;
                p -> left -> parent -> right = p -> left -> right;

                delete temp;

            }
            if(p -> left -> right == NULL && p->left->left !=NULL)
            {
                p -> left -> left -> parent = p -> left -> parent;
                p -> left -> parent -> left = p -> left -> left;

                delete temp;

            }



        }
        if(p -> right != NULL)
        {

            unBeadSpeciesTree(p -> right);

            Node * temp = p -> right;
            if(p -> right -> left == NULL  && p->right->right !=NULL)
            {

                p -> right -> right -> parent = p -> right -> parent;
                p -> right -> parent -> right = p -> right -> right;

                delete temp;
            }
            if(p -> right -> right == NULL  && p->right->left !=NULL)
            {

                p -> right -> left -> parent = p -> right -> parent;
                p -> right -> parent -> left = p -> right -> left;

                delete temp;
            }

        }

    }

}
/*
   void speciesTreeNoBeaded(Node* newnode, int & N, double* s_times, double * s, vector <Node *> v)
   {
   Node ** ar = new Node * [N];
   int tail = 0;


   pushToArray(newnode, tail, ar);
   getRanks(newnode, tail, ar);

   getS(newnode, s_times, v);
   for(int j = 2; j < N; ++j)
   {
   s[j-2] = s_times[j-2] - s_times[j-1];
   }
   int itemp = 0;
   double * arDistFrRoot = new double [N];
   saveDistFrRoot(newnode, arDistFrRoot, itemp);

   delete[] ar;
   delete[] arDistFrRoot;
   }
   */

void readUnrankedTrees(string str_ST, int N, int rounds, vector<int> index_vector, int ** ar_y, int ***k,  double * s, vector <Node *> v, vector<string> gts_vect, double & threshold, string & candidate_str, double x)
{
    std::stack <Node *> stkGTunr;
    int lblUnrGT = 0;
    string ar_strlbl[N-1];
    pushNodesUnrankedGT(lblUnrGT, stkGTunr , str_ST, ar_strlbl);
    Node * pnode = stkGTunr.top();
    getDescTaxa(pnode, N);
    pnode -> rank = 1;
    rankUnrankedTrees(pnode, N, rounds, index_vector, ar_y, k, s, v, gts_vect, threshold, candidate_str, x);

    deleteTree(pnode);
}


int calcNumberOfTaxaBrL(std::string str)
{
    int i = 0;
    int lbl = 0;
    while(i < (int) str.size())
    {

        // last term added to avoid treating e-05 as a label
        if(((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= 'A' && str[i] <= 'Z')) && ((i+1) < (int) str.size()) && (str[i+1] != '-' && str[i+1] != '+'))
        {
            while(((i+1) < (int) str.size()) && str[i+1] != ':')
            {
                ++i;
            }
            ++i;
            while(((i+1) < (int) str.size()) && (str[i+1] != ',' || str[i+1] == ')'))
            {			
                ++i;
            }

            ++lbl;

        }
        else if ( str[i] == ')' )
        {
            ++i;
        }
        ++i;
    }

    return lbl;
}





string pickInitialCandidate(vector<string> sp_vect, vector<string> gts_vect, Node *& newnode, int & N, vector <Node *> v)
{
    int lbl = 0;
    int itemp = 0;
    int tail = 0;
    size_t maxQue;
    string str = "";
    double tempdist = 0.;
    double neg_loglike;
    vector <string> species_vect;
    str = sp_vect[0];
    N = calcNumberOfTaxa(str);
    //N = calcNumberOfTaxaBrL(str);
    //cout << "N " << N << endl;
    for(int i = 0; i < sp_vect.size(); ++i)
    {
        Node * pt;
        vector <string> ranked_vect;
        str = sp_vect[i];
        std::stack <Node *> stkGTunr;
        int lblUnrGT = 0;
        string * ar_strlbl = new string[N-1];
        pushNodesUnrankedGT(lblUnrGT, stkGTunr, str, ar_strlbl);
        pt = stkGTunr.top();
        getDescTaxa(pt, N);
        pt -> rank = 1;
        int prod = 1;
        str = "";
        int temp = 0;
        int tops_count = numberOfRankings(pt, N, prod);
        vector<vector<int>> seq = permuteRanks(pt);
        if(seq.size() != tops_count) cout << "ERROR: not all ranks!" << endl;
        for(int i = 0; i < seq.size(); ++i)
        {
            str = "";
            assignRanks(pt, seq[i]);
            writeRankTreeIntoStr(pt, temp, N, str);
            ranked_vect.push_back(str);
        }

        delete [] ar_strlbl;
        deleteStack(stkGTunr);
        deleteTree(pt);



        vector<int> index_vector;
        vector<int> fin_num_cands;
        int v_size = ranked_vect.size();
        for(int i = 0; i < ranked_vect.size(); ++i)
            index_vector.push_back(i);

        //select subset of rankings	
        //int N_subset = 2*N;
        int N_subset = index_vector.size();
        if(index_vector.size() > N_subset)
        {
            fisher_shuffle(index_vector);
            for(int i = 0; i < N_subset; ++i)
                fin_num_cands.push_back(index_vector[i]);
            v_size = N_subset;
            for(int j = 0;j < N_subset; ++j)
                species_vect.push_back(ranked_vect[fin_num_cands[j]]);
        }
        else
        {
            species_vect = ranked_vect;
        }

    }

    double tmp_array[species_vect.size()];
    for(int i = 0; i < species_vect.size(); ++i)
    {
        str = species_vect[i];
        stack<Node*> stkST;
        lbl = 0;
        pushNodes(lbl, stkST , str);
        Node * newnode = stkST.top();
        newnode -> distFrRoot = 0.;
        tempdist = 0.;
        distFromRoot(newnode);
        isUltrametric(newnode, tempdist);

        tail = 0;
        int ** ar_y = new int * [N];
        for (int i = 0; i < N; i++) ar_y[i] = new int [N];
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                ar_y[i][j] = 0;
            }
        }
        double * arDistFrRoot = new double [N];
        Node ** ar = new Node * [N];

        double s_times[N-1];
        double s[N-2];

        pushToArray(newnode, tail, ar);
        getRanks(newnode, tail, ar);
        if(v.size() != 0)
        {
            while (!v.empty())
            {
                v.pop_back();
            }
        }

        getS(newnode, s_times, v);
        for(int j = 2; j < N; ++j)
        {
            s[j-2] = s_times[j-2] - s_times[j-1];
        }
        itemp = 0;
        saveDistFrRoot(newnode, arDistFrRoot, itemp);

        maxQue = (size_t) N*(1+N)/2;
        makeBeadedTree(newnode, maxQue);
        makeBeadedArray(newnode, ar_y, maxQue);

        vector<int> index_vector;
        for(int i = 0; i < N - 2; ++i)
            index_vector.push_back(i);


        double initial_br_length[N];
        for(int i = 0; i < N - 2; ++i)
            initial_br_length[i] = s[i];

        double res = 0.;
        //int rounds = (int) (N/2.);
        //     int rounds = N_rounds;



        //   vector<vector<vector<double>>> gts_probs_vect;
        //   neg_loglike = newCalcNegativeLike(N, s, v, ar_y, gts_vect, gts_probs_vect);

        // cout << "pick init nll: " << neg_loglike << endl;
        int ***k;
        k = new int **[N+1];
        for(int i = 0; i < N+1; i++)
        {
            k[i] = new int *[N];
            for(int j = 0; j < N; ++j)
                k[i][j] = new int[N+1];
        }

        // for(int j = 0; j < rounds; ++j)
        //  {
        //     fisher_shuffle(index_vector);
        //   for(int index = 0; index < N - 2; ++index)
        // {
        neg_loglike = lbfgs(init_h, s, N, v, ar_y, gts_vect);
        //cout << "pick init neg_loglike = " << neg_loglike << endl;
        // }

        // if(res > neg_loglike) j = rounds;
        //  else neg_loglike = res;
        // }
        tmp_array[i] = neg_loglike;
        for(int i = 0; i < N+1; ++i)
        {
            for(int j = 0; j < N; ++j)
                delete [] k[i][j];
            delete [] k[i];
        }
        delete [] k;
        for(int i = 0; i < N; ++i)
        {
            delete[] ar_y[i];
        }
        delete[] ar_y;
        deleteTree(newnode);

        delete[] ar;
        delete[] arDistFrRoot;
    }

    double tmp = tmp_array[0];
    int min_i = 0;
    for(int i = 1; i < species_vect.size(); ++i)
    {
        if(tmp_array[i] < tmp)
        {
            tmp = tmp_array[i];
            min_i = i;
        }
    }
    //cout << min_i << endl;
    return species_vect[min_i];
}



void writeTreeWithBrLengths(Node * p, int & temp, int N, std::string & str, double * s, double x)
{
    if(p -> label.size() != 0) str += p -> label;
    else
    {
        str += '(';
        writeTreeWithBrLengths(p -> left, temp, N, str, s, x);
        str += ',';
        writeTreeWithBrLengths(p -> right, temp, N, str, s, x);
        str += ')';
    }
    if(p -> rank != 1)
    {
        if(p -> label.size() == 0)
        {
            temp = p -> rank - p -> parent -> rank;
            double sum = 0.;
            for(int i = 0; i < temp; ++i)
                sum += s[p -> rank - 2 - i];

            str += ':';
            str +=  to_string(sum);
        }
        else
        {

            temp = N - p -> parent -> rank - 1;
            for(int i = 0; i < temp; ++i)
                x += s[N-3-i];

            str += ':';
            str +=  to_string(x);

        }
    }

}




void calcLikeWithNNI(int & arg_counter, char* argv[])
{
    //	clock_t sp_cand_start, sp_cand_stop;
    double threshold = 0.;
    Node * newnode;



    ifstream st_file(argv[arg_counter]);
    ++arg_counter;
    vector<string> sp_vect;
    string strR = "";
    while(getline(st_file >> std::ws, strR, ';'))
    {
        if(strR.size() < 3) break;
        removeSpaces(strR);
        sp_vect.push_back(strR);
    }
    st_file.close();
    /*
       ifstream finGT(argv[arg_counter]);
       ++arg_counter;
       vector<string> gts_vect;


       strR = "";
       while(getline(finGT >> std::ws, strR, ';'))
       {
       if(strR.size() < 3) break;
       removeSpaces(strR);
       gts_vect.push_back(strR);
       }
       finGT.close();
       */

    vector<string> gts_vect;

    if(strcmp(argv[arg_counter],"-ugt") == 0)
    {
        ++arg_counter;
        writeRankedUnrTreesIntoVect(arg_counter, argv, gts_vect);
    }
    else if(strcmp(argv[arg_counter],"-rgt") == 0)
    {
        ++arg_counter;
        ifstream finGT(argv[arg_counter]);
        string strR = "";

        while(getline(finGT >> std::ws, strR, ';'))
        {
            if(strR.size() < 3) break;
            removeSpaces(strR);
            gts_vect.push_back(strR);
        }
        finGT.close();
    }

    vector <Node *> v;
    int N = 0;
    string initial_tree_string = pickInitialCandidate(sp_vect, gts_vect, newnode, N, v);
    ofstream ftopi("outInitialTopo.txt", ios::out | ios::app);
    ftopi << initial_tree_string << ";" << endl;
    ftopi.close();


    int ** ar_y = new int * [N];
    for (int i = 0; i < N; i++) ar_y[i] = new int [N];
    for (int i = 0; i < N; ++i)
    {
        for (int j = 0; j < N; ++j)
        {
            ar_y[i][j] = 0;
        }
    }
    double * arDistFrRoot = new double [N];
    Node ** ar = new Node * [N];

    double s_times[N-1];
    double s[N-2];

    double tempdist = 0.;
    int lbl = 0;
    string str = initial_tree_string;
    stack<Node*> stkST;
    pushNodes(lbl, stkST , str);
    newnode = stkST.top();
    newnode -> distFrRoot = 0.;
    distFromRoot(newnode);
    isUltrametric(newnode, tempdist);
    int tail = 0;
    pushToArray(newnode, tail, ar);
    getRanks(newnode, tail, ar);
    if(v.size() != 0)
    {
        while (!v.empty())
        {
            v.pop_back();
        }
    }

    getS(newnode, s_times, v);
    for(int j = 2; j < N; ++j)
    {
        s[j-2] = s_times[j-2] - s_times[j-1];
    }
    int itemp = 0;
    saveDistFrRoot(newnode, arDistFrRoot, itemp);

    int tempo = 0;
    string out_str="";
    double x = 0.1;
    writeTreeWithBrLengths(newnode, tempo, N, out_str, s, x);
    size_t maxQue = (size_t) N*(1+N)/2;
    makeBeadedTree(newnode, maxQue);
    makeBeadedArray(newnode, ar_y, maxQue);

    delete[] arDistFrRoot;

    delete[] ar;

    vector<int> index_vector;
    for(int i = 0; i < N - 2; ++i)
        index_vector.push_back(i);



    //int rounds = (int) (N/2.);
    //
    int rounds = N_rounds;


    //   vector<vector<vector<double>>> gts_probs_vect;
    //    double neg_loglike = newCalcNegativeLike(N, s, v, ar_y, gts_vect, gts_probs_vect);

    // cout << "init neg_loglike = " << neg_loglike << endl;
    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; i++)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }

    //  for(int j = 0; j < rounds; ++j)
    //   {
    //       fisher_shuffle(index_vector);
    //       cout << "round " << j << endl;
    //      for(int index = 0; index < N - 2; ++index)
    //     {

    threshold = lbfgs(init_h, s, N, v, ar_y, gts_vect);
    //cout  << " " << threshold << endl;
    // }

    //   if(threshold > neg_loglike) j = rounds;
    // else neg_loglike = threshold;
    //  }

    double neg_loglike = threshold;
    int temp = 0;
    string candidate_str = "";
    unBeadSpeciesTree(newnode);

    writeTreeWithBrLengths(newnode, temp, N, candidate_str, s, x);
    //cout << "after br lengths optimization tree: " << candidate_str << endl;

    double old_threshold = 0.;
    int nni_counts = 0;
    while(fabs(old_threshold - threshold) > 0.1 && nni_counts < N_nni)
    {
        //cout << "*************************** NNI **************************** " << nni_counts << endl;
        old_threshold = threshold;
        vector<string> str_vect;
        searchInternalNodes(newnode, newnode, N, str_vect);
        deleteTree(newnode);


        for(int i = 0; i < str_vect.size(); ++i)
        {

            //cout << "*************************** nni **************************** " << i << endl;
            readUnrankedTrees(str_vect[i], N, rounds, index_vector, ar_y, k, s, v, gts_vect, threshold, candidate_str, x);
        }
        tempdist = 0.;
        lbl = 0;

        Node ** ar_c = new Node * [N];
        stack<Node*> stkST;
        pushNodes(lbl, stkST , candidate_str);
        newnode = stkST.top();
        newnode -> distFrRoot = 0.;
        distFromRoot(newnode);
        isUltrametric(newnode, tempdist);
        tail = 0;
        pushToArray(newnode, tail, ar_c);
        getRanks(newnode, tail, ar_c);
        if(v.size() != 0)
        {
            while (!v.empty())
            {
                v.pop_back();
            }
        }

        getS(newnode, s_times, v);
        nni_counts++;

        delete[] ar_c;
        //cout << "LIKE: " << threshold << endl;
        //cout << "DIFFERENCE: " << fabs(old_threshold-threshold) << endl;
    }
    // after the first round of nni, select subset of nni's
    /*
       while(fabs(old_threshold - threshold) > 0.1 && nni_counts < 5)
       {

       old_threshold = threshold;
       vector<string> str_vect;
       searchInternalNodes(newnode, newnode, N, str_vect);
       deleteTree(newnode);

       vector<int> shuffle_nni;
       for(int i = 0; i < 2*N - 4; ++i)
       shuffle_nni.push_back(i);

       fisher_shuffle(shuffle_nni);


       for(int i = 0; i < N - 2; ++i)//keep half nni's
       {
       cout << "------------------species tree " << i << endl;
       readUnrankedTrees(str_vect[shuffle_nni[i]], N, rounds, index_vector, ar_y, k, s, gts_vect, threshold, candidate_str, x);
       }
       tempdist = 0.;
       lbl = 0;

       Node ** ar_c = new Node * [N];
       stack<Node*> stkST;
       pushNodes(lbl, stkST , candidate_str);
       newnode = stkST.top();
       newnode -> distFrRoot = 0.;
       distFromRoot(newnode);
       isUltrametric(newnode, tempdist);
       tail = 0;
       pushToArray(newnode, tail, ar_c);
       getRanks(newnode, tail, ar_c);
       nni_counts++;

       delete[] ar_c;
       cout << "----------NNI_count: " << nni_counts << " diff: " << fabs(old_threshold-threshold) << endl;

       }
       */
    //  cout << "total nni: " << nni_counts << endl;
    cout << "Negative log-likelihood  = " << threshold << endl;
    ofstream ftop("outWithNniMLTopo.txt", ios::out | ios::app);
    ftop << candidate_str << ";" << endl;
    ftop.close();
    //cout << "inferred tree: " << candidate_str << ";" << endl;
    //	for(int t = 0; t < N - 2; ++t)
    //  {
    //    cout << "s[" << t << "] = " << s[t] << endl;
    //  }

    /*
       ofstream fin("s_i.txt");
       for(int i = 0; i < N - 2; ++i)
       {
       fin <<  s[i] << "\t";
       }
       fin << endl;
       fin.close();
       */
    for(int i = 0; i < N; ++i)
    {
        delete[] ar_y[i];
    }
    delete[] ar_y;

    for(int i = 0; i < N+1; ++i)
    {
        for(int j = 0; j < N; ++j)
            delete [] k[i][j];
        delete [] k[i];
    }
    delete [] k;

    deleteTree(newnode);

}




void initProbsVect(int * m, int *** k, double * s, double * coal, int & n, vector<vector<double>> & vect)
{
    vector<double> res;
    for(int i = 2; i < n; ++i)
    {
        res.push_back(conditionalProbability(i, m, k, s));
    }
    res.push_back(coal[m[0] - 1]);
    vect.push_back(res);
}

double getGT(int & N, double * s, vector <Node *> v, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k, vector<vector<double>> & vect)
{

    double arhists[MAX_EXPONENT_VALUE];
    for(int i = 0; i < MAX_EXPONENT_VALUE; ++i)
        arhists[i] = 0.;
    returnMrca(v, newnodeGT, N, arMrca);

    int * arRankHistory = new int [N-1];
    get_hseq(arRankHistory, newnodeGT, v, N);
    reverseArraySort(arRankHistory, N);
    /*  cout << "Max Ranked History" << endl;
        for(int i = 0; i < N-1; ++i) cout << arRankHistory[i] << " ";
        cout << endl;
        */

    int *  arRankHistoryCopy = new int[N-1];
    for(int i = 0; i < N-1; ++i)
    {
        arRankHistoryCopy[i] = arRankHistory[i];
    }
    int * m_i = new int[N-1];
    arrayFreq(arRankHistoryCopy, m_i, N-1);//1st arg changes after running this func
    int max_m_i = *max_element(m_i, m_i+N-1);
    get_node(arRankHistory, ar_rankns, N, arMrca);
    for(int i = 0; i < N+1; i++)
        for(int j = 0; j < N; ++j)
            for(int z = 0; z < N+1 ; ++z)
            {
                k[i][j][z] = 0;
            }
    calc_k_ijz(k, N, m_i, ar_y, ar_rankns);
    int hist_counter = 0;
    initProbsVect(m_i, k, s, array_invcoal, N, vect);
    double onehistprob = 1.;
    for(int i = 0; i < N - 1; ++i)
    {
        // cout << hist_counter << " " << i << " " << vect[hist_counter][i] << endl;
        onehistprob *= vect[hist_counter][i];
    }
    //cout << "99999-*_**_*_*_*_*_*_**_*_" <<  (int) floor(fabs(log2(fabs(onehistprob)))) << endl;
    if(onehistprob != 0 && ((int) floor(fabs(log2(fabs(onehistprob))))) < MAX_EXPONENT_VALUE)
        arhists[ (int) floor(fabs(log2(fabs(onehistprob))))] +=  fabs(onehistprob);



    int * arRHcopy = new int[N-1];
    for(int i = 0; i < N - 1; ++i)
        arRHcopy[i] = arRankHistory[i];

    int numb = 0;
    int * next_history;
    int history_counter = 1;
    while(*(arRHcopy + N - 2) != 1)
    {
        hist_counter++;
        next_history = getNextHistory (arRHcopy, arRankHistory, N, history_counter, numb);
        /* cout << "---------- Next Label history ---------- " << endl;
           for(int i = 0; i < N - 1; ++i)
           cout << next_history[i] << " ";
           cout << endl;
           */
        for(int i = 0; i < N-1; ++i)
            arRankHistoryCopy[i] = next_history[i];
        arrayFreq(arRankHistoryCopy, m_i, N-1);
        max_m_i = *max_element(m_i, m_i+N-1);
        get_node(next_history, ar_rankns, N, arMrca);

        for(int i = 0; i < N+1; i++)
        {
            for(int j = 0; j <= max_m_i ; ++j)
            {
                for(int z = 0; z < N+1 ; ++z)
                {
                    k[i][j][z] = 0;
                }
            }
        }
        calc_k_ijz(k, N, m_i, ar_y, ar_rankns);
        initProbsVect(m_i, k, s, array_invcoal, N, vect);
        onehistprob = 1.;
        for(int i = 0; i < N - 1; ++i)
        {
            onehistprob *= vect[hist_counter][i];
        }
        if(onehistprob != 0 && ((int) floor(fabs(log2(fabs(onehistprob))))) < MAX_EXPONENT_VALUE)
            arhists[ (int) floor(fabs(log2(fabs(onehistprob))))] +=  fabs(onehistprob);
        history_counter++;
    }
    double fin_prob = 0.;
    for(int i = MAX_EXPONENT_VALUE - 1; i >= 0 ; --i)
    {
        if(arhists[i] != 0)
        {
            fin_prob += arhists[i];
        }
    }
    delete[] arRankHistory;
    delete[] arRHcopy;
    delete[] arRankHistoryCopy;
    delete[] m_i;
    return fin_prob;
}


double updateNegativeLike(int & N, double* s, vector <Node *> v, int** ar_y, int ***k, vector<string> & gts_vect, vector<vector<vector<double>>> & vect, int & interval)
{
    string strGT="";
    Node ** arMrca = new Node * [N-1];
    int ** ar_rankns = new int * [N-1];
    for (int i = 0; i < N-1; i++) ar_rankns[i] = new int [N-1];
    for(int i = 0; i < N-1; ++i)
    {
        for(int j = 0; j < N - 1; ++j)
        {
            ar_rankns[i][j] = 0;
        }
    }
    Node * newnodeGT;
    double prec_total_prob = 0.;
    double one_gt_prob = 0.;
    //double total = 0.;
    for(int count = 0; count < gts_vect.size(); ++count)
    {
        strGT = gts_vect[count];
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
        one_gt_prob = updateProbability(N, s, v, newnodeGT, ar_y, arMrca, ar_rankns, k, vect, interval, count);
        // total += one_gt_prob;
        prec_total_prob += log(one_gt_prob);

        delete [] arGT;

        deleteTree(newnodeGT);
        newnodeGT = NULL;
        deleteStack(stkGT);

    }
    // cout << "total" << " " << -prec_total_prob << endl;

    for(int i = 0; i < N-1; ++i)
        delete[] ar_rankns[i];
    delete [] ar_rankns;

    delete[] arMrca;


    return -prec_total_prob;
}


double updateHistoryProbability(vector<vector<vector<double>>> & vect, int & hist_counter, int & n, int & gt_count)
{
    double res = 1.;
    for(int i = 2; i < n+1; ++i)
    {
        res *= vect[gt_count][hist_counter][i-2];
    }
    return res;
}


double updateProbability(int & N, double * s, vector <Node *> v, Node * newnodeGT, int ** ar_y, Node ** arMrca, int ** ar_rankns, int ***k, vector<vector<vector<double>>> & vect, int & interval, int & gt_count)
{

    double arhists[MAX_EXPONENT_VALUE];
    for(int i = 0; i < MAX_EXPONENT_VALUE; ++i)
        arhists[i] = 0.;

    returnMrca(v, newnodeGT, N, arMrca);

    int * arRankHistory = new int [N-1];
    get_hseq(arRankHistory, newnodeGT, v, N);
    reverseArraySort(arRankHistory, N);
    /* cout << "Max Ranked History" << endl;
       for(int i = 0; i < N-1; ++i) cout << arRankHistory[i] << " ";
       cout << endl;
       */

    int *  arRankHistoryCopy = new int[N-1];

    for(int i = 0; i < N-1; ++i)
    {
        arRankHistoryCopy[i] = arRankHistory[i];
    }

    int * m_i = new int[N-1];
    arrayFreq(arRankHistoryCopy, m_i, N-1);//1st arg changes after running this func
    int max_m_i = *max_element(m_i, m_i+N-1);
    get_node(arRankHistory, ar_rankns, N, arMrca);
    for(int i = 0; i < N+1; i++)
        for(int j = 0; j < N; ++j)
            for(int z = 0; z < N+1 ; ++z)
            {
                k[i][j][z] = 0;
            }
    calc_k_ijz(k, N, m_i, ar_y, ar_rankns);
    int hist_counter = 0;

    vect[gt_count][hist_counter][interval-2] = conditionalProbability(interval, m_i, k, s);
    double res = 1.;
    for(int i = 2; i < N+1; ++i)
    {
        res *= vect[gt_count][hist_counter][i-2];
    }
    double onehistprob = res;

    if(onehistprob != 0 && ((int) floor(fabs(log2(fabs(onehistprob))))) < MAX_EXPONENT_VALUE)
        arhists[ (int) floor(fabs(log2(fabs(onehistprob))))] +=  fabs(onehistprob);
    int * arRHcopy = new int[N-1];
    for(int i = 0; i < N - 1; ++i)
        arRHcopy[i] = arRankHistory[i];

    int numb = 0;
    int * next_history;
    int history_counter = 1;
    while(*(arRHcopy + N - 2) != 1)
    {
        hist_counter++;
        next_history = getNextHistory (arRHcopy, arRankHistory, N, history_counter, numb);
        /*		  cout << "---------- Next Label history ---------- " << endl;
                  for(int i = 0; i < N - 1; ++i)
                  cout << next_history[i] << " ";
                  cout << endl;
                  */

        for(int i = 0; i < N-1; ++i)
            arRankHistoryCopy[i] = next_history[i];
        arrayFreq(arRankHistoryCopy, m_i, N-1);
        max_m_i = *max_element(m_i, m_i+N-1);
        get_node(next_history, ar_rankns, N, arMrca);

        for(int i = 0; i < N+1; i++)
        {
            for(int j = 0; j <= max_m_i ; ++j)
            {
                for(int z = 0; z < N+1 ; ++z)
                {
                    k[i][j][z] = 0;
                }
            }
        }
        calc_k_ijz(k, N, m_i, ar_y, ar_rankns);


        vect[gt_count][hist_counter][interval-2] = conditionalProbability(interval, m_i, k, s);
        //onehistprob = updateHistoryProbability(vect, hist_counter, N, gt_count);
        res = 1.;
        for(int i = 2; i < N+1; ++i)
        {
            res *= vect[gt_count][hist_counter][i-2];
        }
        onehistprob = res;
        if(onehistprob != 0 && ((int) floor(fabs(log2(fabs(onehistprob))))) < MAX_EXPONENT_VALUE)
            arhists[ (int) floor(fabs(log2(fabs(onehistprob))))] +=  fabs(onehistprob);
        history_counter++;
    }
    double fin_prob = 0.;
    for(int i = MAX_EXPONENT_VALUE - 1; i >= 0 ; --i)
    {
        if(arhists[i] != 0)
        {
            fin_prob += arhists[i];
        }
    }

    delete[] arRankHistory;
    delete[] arRHcopy;
    delete[] arRankHistoryCopy;
    delete[] m_i;

    return fin_prob;
}

double lbfgs(double h, double * s, int & N, vector <Node *> vlabels, int ** ar_y, vector<string> gts_vect)
{
    /* Local variables */
    static double f, g[1024];
    static integer i__;
    static double l[1024];
    static integer m, n;
    static double u[1024], x[1024], t1, t2, wa[43251];
    static integer nbd[1024], iwa[3072];
    static integer taskValue;
    static integer *task=&taskValue; /* must initialize !! */
    static double factr;
    static integer csaveValue;
    static integer *csave=&csaveValue;
    static double dsave[29];
    static integer isave[44];
    static logical lsave[4];
    static double pgtol;
    static integer iprint;


    /*     We wish to have output at every iteration. */
    iprint = -1; 
    /*     We specify the tolerances in the stopping criteria. */
    factr = 1e7;
    pgtol = 1e-5;
    /*        m of limited memory corrections stored.  (n and m should not */
    /*        exceed the limits nmax and mmax respectively.) */
    n = N-2;
    m = 5;
    /*                    l   specifies the lower bounds, */
    /*                    u   specifies the upper bounds. */
    for (int i = 0; i < n; i++) 
    {
        x[i] = (init_a + init_b)/2.; // init val, can put any
        nbd[i] = 2;
        l[i] = init_a;
        u[i] = init_b;
    }

    *task = (integer)START;
L111:
    setulb(&n, &m, x, l, u, nbd, &f, g, &factr, &pgtol, wa, iwa, task, &
            iprint, csave, lsave, isave, dsave);
    if ( IS_FG(*task) ) {

        //s[index] = x[0];

        f = calcNegativeLike(N, x, vlabels, ar_y, gts_vect);

        // TODO REMOVE int interval = index + 2;

        for (int i = 0; i < n; i++) 
        {
            for (int j = 0; j < n; j++) 
            {
                s[j] = x[j];
            }
            s[i] = x[i] + h;
            g[i] = calcNegativeLike(N, s, vlabels, ar_y, gts_vect);
            g[i] = (g[i] - f) / h;
        }
        goto L111;
    }

    if ( *task==NEW_X ) {
        goto L111;
    }
    return f;
} 


/*
   double cpu_brent_improved(double a, double b, double * s, int index, int & N, vector <Node *> vlabels, int ** ar_y, int ***k, vector<string> gts_vect, vector<vector<vector<double>>> & vect)
   {
   double q, c, f_q, f_a, f_b, f_c;

   s[index] = a;
   int interval = index + 2;
   f_a = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
   s[index] = b;
   f_b = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
   cout << f_a << " " << f_b << endl; 
//   if(a>b || f_a * f_b >= 0)   // The algorithm assumes a<b
//   {
//       return -1; // Root not bound by the given guesses
//   }

do{
c = (a+b)/2;
s[index] = c;
f_c = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
cout << interval << " " << c << " " << f_c << endl;
if(fabs(f_a-f_c) > brent_tolerance && fabs(f_b-f_c) > brent_tolerance)  // f(a)!=f(c) and f(b)!=f(c)
{
// Inverse quadratic interpolation
q = a*f_b*f_c/((f_a-f_b)*(f_a-f_c)) + b*f_a*f_c/((f_b-f_a)*(f_b-f_c)) + c*f_a*f_b/((f_c-f_a)*(f_c-f_b));
}else{
// Secant rule
q = b - f_b*(b-a)/(f_b-f_a);
}


s[index] = q;
cout << "q: " << q << endl;
f_q = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
if(c>q)
{
swap(q, c);
}
if(f_c*f_q < 0)
{
a = q;
b = c;
}else if(f_q*f_b < 0){
a = c;
}else{
b = q;
}
}while(f_q < brent_tolerance || fabs(b-a)<brent_tolerance); // Convergence conditions
cout << "f_q: " << f_q << endl; 
return f_q;
}


// fp = (f(x + h) - f(x)) / h
// fpp = (f(x+h) - 2*f(x) + f(x-h)) / h /h;
double newtonRapson(double a, double h, double * s, int index, int & N, vector <Node *> vlabels, int ** ar_y, int ***k, vector<string> gts_vect, vector<vector<vector<double>>> & vect)
{ 
int interval = index + 2;
// s[index] = a;
double fx;
//  s[index] = a + h;
double fxph;
// s[index] = a - h;
double fxmh;
double fp;
double fpp;
double val = 0;
double new_val = a; 
while(fabs(val-new_val) >= brent_tolerance)
{

val = new_val;
s[index] = val;
fx = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
s[index] = val + h;
fxph = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
s[index] = val - h;
fxmh = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);

fp = (fxph - fx)/h;
fpp =  (fxph - 2*fx + fxmh)/h/h;
new_val = val - fp/fpp;
cout << fabs(val-new_val) << " " << fp/fpp << " " << new_val << endl;
}

s[index] = new_val;
return updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
}
*/


