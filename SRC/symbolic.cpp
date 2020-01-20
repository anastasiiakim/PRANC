#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <iomanip> // output more decimals
#define EPS 1e-8
#include "string_manipulations.h"
#include "node.h"
#include "queue.h"
#include "probs_calculation.h"
#include "symbolic.h"


using namespace std;



double symbolicIntervalProbability(int i, int * m_i, int *** array_k, double * s, ofstream & fout) 
{
    double res = 0.;
    double prod;
    double m = m_i[i-1];
    for(int j = 0; j < m + 1; ++j)
    {
        prod = 1.;
        for(int k = 0; k < m + 1; ++k)
        {
            if(k != j)
            {
                prod *= lambda_ij(i, k, array_k) - lambda_ij(i, j, array_k);
            }
        }

    if(j == 0 && i == 2) fout << " + ";
    if(j == 0) fout << "(";
      fout << "exp("  << -lambda_ij(i, j, array_k) << "*(s" << i-1 << "-s" << i << "))*1/(" << prod << ")";
    if (j < m) fout << " + ";
    if (j == m) fout << ") ";
        
        res += (double) exp(-lambda_ij(i, j, array_k)*s[i-2])/prod;
    }

    fout << " * " << endl;
    return res;
}

double symbolicGeneTreeProbability(int * m, int *** k, double * s, double * coal, int n, ofstream & fout)
{
  double res = 1.;
  for(int i = 2; i < n; ++i)
  {
    res *= symbolicIntervalProbability(i, m, k, s, fout);
  }
  fout << pow(2, m[0]) << "/" << factorial(m[0])*factorial(m[0]+1) << endl;

  return res*coal[m[0] - 1];
}



double getOneSymbolicGeneTreeProb(int N, double * s,  Node * newnode, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k, ofstream & fout, ofstream & histprobs_file)
{
  
  returnMrca(newnode, newnodeGT, N, arMrca);

  int * arRankHistory = new int [N-1];
  get_hseq(arRankHistory, newnodeGT, newnode, N);
  reverseArraySort(arRankHistory, N);


  int *  arRankHistoryCopy = new int[N-1];
  for(int i = 0; i < N-1; ++i)
  {
    arRankHistoryCopy[i] = arRankHistory[i];
  }
  int * m_i = new int[N-1];
  arrayFreq(arRankHistoryCopy, m_i, N-1);
  int max_m_i = *max_element(m_i, m_i+N-1);
  get_node(arRankHistory, ar_rankns, N, arMrca);
  for(int i = 0; i < N+1; i++)
    for(int j = 0; j < N; ++j)
      for(int z = 0; z < N+1 ; ++z)
      {
        k[i][j][z] = 0;
      }
  calc_k_ijz(k, N, m_i, ar_y, ar_rankns);
  double onehistprob = symbolicGeneTreeProbability(m_i, k, s, array_invcoal, N, fout);
  double prob_val = onehistprob;

  for(int i = 0; i < N-1; ++i) histprobs_file << arRankHistory[i];
    histprobs_file << '\t';
  histprobs_file << onehistprob << endl;

  int * arRHcopy = new int[N-1];
  for(int i = 0; i < N - 1; ++i)
    arRHcopy[i] = arRankHistory[i];

  int history_counter = 1;
  int numb = 0;
  int * next_history;

  while(*(arRHcopy + N - 2) != 1)
  {
    next_history = getNextHistory (arRHcopy, arRankHistory, N, history_counter, numb);
    fout << endl;


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
    onehistprob = symbolicGeneTreeProbability(m_i, k, s, array_invcoal, N, fout);
    prob_val += onehistprob;

    for(int i = 0; i < N-1; ++i) histprobs_file << next_history[i];
      histprobs_file << '\t';        
    histprobs_file << onehistprob << endl;
    history_counter++;
  }
 

  delete[] arRankHistory;
  delete[] arRHcopy;
  delete[] arRankHistoryCopy;
  delete[] m_i;

  return prob_val;
}




double getSymbolicGeneTreeProb(int & arg_counter, char* argv[], int & N, Node * newnode, double* s, int** ar_y)
{
  string strGT="";
  ifstream finGT(argv[arg_counter]); 
  ++arg_counter;
  ofstream fout("outSymbolic.txt");
  ofstream histprobs_file("outHistProbs.txt");

  double total = 0.;
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
    total += getOneSymbolicGeneTreeProb(N, s,  newnode, newnodeGT, ar_y, array_invcoal, arMrca, ar_rankns, k, fout, histprobs_file);

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

  finGT.close();
  histprobs_file.close();
  fout.close();

  return total;

}



void symbolicProbsRankedGtInput(int &  arg_counter, char* argv[])
{
  Node * newnode;
  int N = getNumberOfTaxa(arg_counter, argv, newnode);

  double * s_times = new double [N-1];
  double * s = new double [N-2];
  int ** ar_y = new int * [N];
  for (int i = 0; i < N; i++) ar_y[i] = new int [N];

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
    {
      ar_y[i][j] = 0;
    }
  
  speciesTreeProcessing(newnode, N, s_times, s, ar_y);
  cout << "Total: " << getSymbolicGeneTreeProb(arg_counter, argv, N, newnode, s, ar_y) << endl;

  for(int i = 0; i < N; ++i)
  {
    delete[] ar_y[i];
  }
  delete[] ar_y;
  delete[] s_times;
  delete[] s;
  deleteTree(newnode);
}



