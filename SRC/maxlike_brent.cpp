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
#include "maxlike_brent.h"

#define MAX_EXPONENT_VALUE 1024
#define brent_epsilon 1e-06
#define brent_a 0.001
#define brent_b 6
#define brent_tolerance 1e-06

using namespace std;


//brent() original FORTRAN77 version by Richard Brent. Richard Brent, Algorithms for Minimization Without Derivatives, Dover, 2002

double fbrent (double a, double b, double * s, int index, int & N, vector <Node *> vlabels, int ** ar_y, int ***k, vector<string> gts_vect, vector<vector<vector<double>>> & vect)
{
    double c;
    double d;
    double e;
    double fu;
    double fv;
    double fw;
    double fx;
    double m;
    double p;
    double q;
    double r;
    double sa;
    double sb;
    double t2;
    double tol;
    double u;
    double v;
    double w;
    double x;

    c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );


    sa = a;
    sb = b;
    x = sa + c * ( b - a );
    w = x;
    v = w;
    e = 0.0;

    s[index] = x;
    int interval = index + 2;
    fx = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
    fw = fx;
    fv = fw;

    for ( ; ; )
    {
        m = 0.5 * ( sa + sb ) ;
        tol = brent_epsilon * fabs ( x ) + brent_tolerance;
        t2 = 2.0 * tol;

        if ( fabs ( x - m ) <= t2 - 0.5 * ( sb - sa ) )
        {
            break;
        }

        r = 0.0;
        q = r;
        p = q;

        if ( tol < fabs ( e ) )
        {
            r = ( x - w ) * ( fx - fv );
            q = ( x - v ) * ( fx - fw );
            p = ( x - v ) * q - ( x - w ) * r;
            q = 2.0 * ( q - r );
            if ( 0.0 < q )
            {
                p = - p;
            }
            q = fabs ( q );
            r = e;
            e = d;
        }

        if ( fabs ( p ) < fabs ( 0.5 * q * r ) &&
                q * ( sa - x ) < p &&
                p < q * ( sb - x ) )
        {

            d = p / q;
            u = x + d;

            if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
            {
                if ( x < m )
                {
                    d = tol;
                }
                else
                {
                    d = - tol;
                }
            }
        }

        else
        {
            if ( x < m )
            {
                e = sb - x;
            }
            else
            {
                e = sa - x;
            }
            d = c * e;
        }

        if ( tol <= fabs ( d ) )
        {
            u = x + d;
        }
        else if ( 0.0 < d )
        {
            u = x + tol;
        }
        else
        {
            u = x - tol;
        }

        s[index] = u;
        fu = updateNegativeLike(N, s, vlabels, ar_y, k, gts_vect, vect, interval);
        //        cout << u << " " << fu << endl;
        if ( fu <= fx )
        {
            if ( u < x )
            {
                sb = x;
            }
            else
            {
                sa = x;
            }
            v = w;
            fv = fw;
            w = x;
            fw = fx;
            x = u;
            fx = fu;
        }
        else
        {
            if ( u < x )
            {
                sa = u;
            }
            else
            {
                sb = u;
            }

            if ( fu <= fw || w == x )
            {
                v = w;
                fv = fw;
                w = u;
                fw = fu;
            }
            else if ( fu <= fv || v == x || v == w )
            {
                v = u;
                fv = fu;
            }
        }
    }
    return fx;
}



void brent_rankUnrankedTrees(Node * spnode, int Numtaxa, int rounds, vector<int> index_vector, int ** ar_y, int ***k, double * s, vector <Node *> & vlabels, vector<string> gts_vect, double & threshold, string & candidate_str, double x)
{
    int prod = 1;
    //	double a = 0.0001;
    //	double b = 6.;
    //	double t = 1e-10;

    size_t maxQue = (size_t) Numtaxa*(1+Numtaxa)/2;

    int tops_count = numberOfRankings(spnode, Numtaxa, prod);
    vector<vector<int>> seq = permuteRanks(spnode);
    double s_array[seq.size()][Numtaxa-2];
    double res[seq.size()];

    if(seq.size() != tops_count) cout << "ERROR: not all ranks!" << endl;
    // int NR = (int) (Numtaxa/2.); // how many rankings to select
    int NR = Numtaxa; // how many rankings to select
    double neg_loglike;

    if(seq.size() > NR)
    {
        //select subset of ranks
        vector<int> shuffle_rankings;
        for(int i = 0; i < seq.size(); ++i)
            shuffle_rankings.push_back(i+1);

        fisher_shuffle(shuffle_rankings);

        for(int i = 0; i < NR; ++i)
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

            double temp;

            vector<vector<vector<double>>> vect;
            neg_loglike = newCalcNegativeLike(Numtaxa, s, vlabels,  ar_y, gts_vect, vect);

            double min_temp = 0.;
            for(int j = 0; j < rounds; ++j)
            {
                fisher_shuffle(index_vector);
                for(int index = 0; index < Numtaxa - 2; ++index)
                {
                    temp = fbrent(brent_a, brent_b, s, index_vector[index], Numtaxa, vlabels, ar_y, k, gts_vect, vect);

                }

                //		if(temp > neg_loglike)
                //			j = rounds;
                //		else
                //			neg_loglike = temp;
                //  else
                // {
                if(j == 0)
                {
                    min_temp = neg_loglike;
                }
                if(temp <= min_temp)
                {
                    min_temp = temp;
                    neg_loglike = min_temp;
                }
                // }
            }




            unBeadSpeciesTree(spnode);
            res[i] = neg_loglike;
            for(int t = 0; t < Numtaxa-2; ++t)
                s_array[i][t] = s[t];
        }




        double min_res = res[0];
        int j = 0;
        for(int i = 0; i < NR; ++i)
        {
            if(res[i] < min_res)
            {
                min_res =  res[i];
                j = i;
            }
        }
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
    else
    {
        for(int i = 0; i < seq.size(); ++i)
        {

            if(vlabels.size() != 0)
            {
                while (!vlabels.empty())
                {
                    vlabels.pop_back();
                }
            }


            assignRanks(spnode, seq[i]);
            distFrRootBasedOnRanks(spnode, Numtaxa, vlabels);



            makeBeadedTree(spnode, maxQue);
            makeBeadedArray(spnode, ar_y, maxQue);

            double temp;


            vector<vector<vector<double>>> vect;
            neg_loglike = newCalcNegativeLike(Numtaxa, s, vlabels, ar_y, gts_vect, vect);

            double min_temp = 0.;
            for(int j = 0; j < rounds; ++j)
            {
                fisher_shuffle(index_vector);
                for(int index = 0; index < Numtaxa - 2; ++index)
                {
                    temp = fbrent(brent_a, brent_b, s, index_vector[index], Numtaxa, vlabels, ar_y, k, gts_vect, vect);
                }


                /*
                   if(temp > neg_loglike)
                   j = rounds;
                   else
                   neg_loglike = temp;
                   */
                // cout << "round: " << j << " brent: " << temp << " negLL: " << neg_loglike << endl;
                //   if(temp > neg_loglike)
                //   {
                //     j = rounds;
                //   }
                //		else
                //  {
                if(j == 0)
                {
                    min_temp = neg_loglike;
                }
                if(temp <= min_temp)
                {
                    min_temp = temp;
                    neg_loglike = min_temp;
                }
                //  }

            }

            unBeadSpeciesTree(spnode);
            res[i] = neg_loglike;
            for(int t = 0; t < Numtaxa-2; ++t)
            {
                s_array[i][t] = s[t];
                //  cout << "s[" << t << "] = " << s[t] << endl;
            }

            //  cout << "fin negll for ith tree " << res[i] << endl;
        }


        double min_res = res[0];
        int j = 0;
        for(int i = 0; i < seq.size(); ++i)
        {
            if(res[i] < min_res)
            {
                min_res =  res[i];
                j = i;
            }
        }

        // cout << "best ranking j: " << j << endl;
        // cout << "min_res = " << min_res << " threshold = " << threshold << endl;
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
            assignRanks(spnode, seq[j]);
            distFrRootBasedOnRanks(spnode, Numtaxa, vlabels);
            int temp = 0;
            candidate_str = "";
            for(int t = 0; t < Numtaxa - 2; ++t)
            {
                s[t] = s_array[j][t];
                //   cout << "s[" << t << "] = " << s[t] << endl;
            }
            writeTreeWithBrLengths(spnode, temp, Numtaxa, candidate_str, s, x);
            //  cout << "fin_tree: " << candidate_str << endl;
        }
    }
}




string brent_pickInitialCandidate(vector<string> sp_vect, vector<string> gts_vect, Node *& newnode, int & N, vector <Node *> v)
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
        if(index_vector.size() > N)
        {
            fisher_shuffle(index_vector);
            for(int i = 0; i < N; ++i)
                fin_num_cands.push_back(index_vector[i]);
            v_size = N;
            for(int j = 0;j < N; ++j)
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
        int rounds = N;
        vector<vector<vector<double>>> gts_probs_vect;
        neg_loglike = newCalcNegativeLike(N, s, v, ar_y, gts_vect, gts_probs_vect);
        // cout << "nll: " << neg_loglike << endl;
        int ***k;
        k = new int **[N+1];
        for(int i = 0; i < N+1; i++)
        {
            k[i] = new int *[N];
            for(int j = 0; j < N; ++j)
                k[i][j] = new int[N+1];
        }

        for(int j = 0; j < rounds; ++j)
        {
            fisher_shuffle(index_vector);
            for(int index = 0; index < N - 2; ++index)
            {
                res = fbrent(brent_a, brent_b,  s, index_vector[index], N, v, ar_y, k, gts_vect, gts_probs_vect);

            }

            if(res > neg_loglike) j = rounds;
            else neg_loglike = res;
        }
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
    return species_vect[min_i];
}





void brent_calcLikeWithNNI(int & arg_counter, char* argv[])
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
    int rounds = N;
    vector<vector<vector<double>>> gts_probs_vect;
    double neg_loglike = newCalcNegativeLike(N, s, v, ar_y, gts_vect, gts_probs_vect);
    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; i++)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }

    for(int j = 0; j < rounds; ++j)
    {
        fisher_shuffle(index_vector);

        for(int index = 0; index < N - 2; ++index)
        {

            threshold = fbrent(brent_a, brent_b, s, index_vector[index], N, v, ar_y, k, gts_vect, gts_probs_vect);

        }

        if(threshold > neg_loglike) j = rounds;
        else neg_loglike = threshold;
    }

    int temp = 0;
    string candidate_str = "";
    unBeadSpeciesTree(newnode);

    writeTreeWithBrLengths(newnode, temp, N, candidate_str, s, x);
    //cout << "after br lengths optimization tree: " << candidate_str << endl;

    double old_threshold = 0.;
    int nni_counts = 0;
    while(fabs(old_threshold - threshold) > 0.1 && nni_counts < 5)
    {
        old_threshold = threshold;
        vector<string> str_vect;
        searchInternalNodes(newnode, newnode, N, str_vect);
        deleteTree(newnode);


        for(int i = 0; i < str_vect.size(); ++i)
        {
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
        cout << "LIKE: " << threshold << endl;
        cout << "DIFFERENCE: " << fabs(old_threshold-threshold) << endl;
    }
    ofstream ftop("outWithNniMLTopo.txt", ios::out | ios::app);
    ftop << candidate_str << ";" << endl;
    ftop.close();

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



