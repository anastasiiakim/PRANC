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


void pushToArray (Node * p, int & tail, Node ** ar)
{
  int i = tail;
  ar[i] = p;
  while((i > 0) && (ar[i] -> time < ar[i - 1] -> time))
  {
    Node * tempe;
    tempe = ar[i];
    ar[i] = ar[i-1];
    ar[i-1] = tempe;
    i--;
  }

  tail ++;
}

void pushToArrayNoSort (Node * p, int & tail, Node ** ar)
{
  int i = tail;
  ar[i] = p;
  tail ++;
}



Node * popFromArray (int & tail, Node ** ar)
{
  tail --;
  return ar[tail];
}



void distFromRoot (Node * p)
{
  if (p -> label.size() == 0)
  {
    p -> right -> distFrRoot = p -> distFrRoot + p -> right -> length;
    p -> left -> distFrRoot = p -> distFrRoot + p -> left -> length;
    p -> right -> parent = p;
    p -> left -> parent = p;

    distFromRoot(p -> right);
    distFromRoot(p -> left);
  }
  else // for leaves
  {
    p -> outdegree = 0;
    p -> population = 1;
  }

}


void isUltrametric (Node * p, double & tempdist)
{
  if (p == NULL) cout << "Tree is empty" << endl;
  if (p -> label.size() == 0)
  {
    isUltrametric(p -> right, tempdist);
    isUltrametric(p -> left, tempdist);
  }

  else
  {
    tempdist = p -> distFrRoot;
  }
}

void saveDistFrRoot(Node * p, double * ar, int & i)
{
  if (p -> label.size() == 0)
  {
    saveDistFrRoot(p -> right, ar, i);
    saveDistFrRoot(p -> left, ar, i);
  }
  else
  {
    ar[i] = p -> distFrRoot;
    ++i;
  }

}

void printIsUltrametric (double * ar, int N)
{
  for(int i = 0; i < N - 1; ++i)
    if (fabs(ar[i] - ar[i+1]) > EPS)
    {
      cout << "WARNING: The species tree is not ultrametric. At least two branches differ by eps of 1e-08." << endl;
      return;
    }
}




void deleteTree (Node *  p)
{
  if(p == NULL )
  {
    return;
  }

  if(p -> left != NULL)
  {
    deleteTree(p -> left);
  }

  if(p -> right != NULL)
  {
    deleteTree(p -> right);
  }
  delete p;
}

void deleteStack (std::stack <Node *> stk)
{
  while(!stk.empty())
  {
    stk.pop();
  }
}




void pushNodes (int & lbl, std::stack <Node *> & stk, std::string str)
{
  int i = 0;
  while(i < (int) str.size())
  {

    // last term added to avoid treating e-05 as a label
    if(((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= 'A' && str[i] <= 'Z')) && ((i+1) < (int) str.size()) && (str[i+1] != '-' && str[i+1] != '+'))
    {
      std::string strtemp;
      strtemp += str[i];
      while(((i+1) < (int) str.size()) &&  str[i+1] != ':')
      {
        strtemp += str[i+1];
        ++i;
      }
      Node * newnode = new Node();
      newnode -> label = strtemp;
      newnode -> time = 0.;
      ++lbl;
      ++i;
      if(((i+1) < (int) str.size()) && (static_cast<int> (str[i+1] - '0') >= 0 && static_cast<int> (str[i+1] - '0') <= 9))
      {

        newnode -> length = getBranchLengths(i+1, str);
      }
      stk.push(newnode);

    }

    else if ( str[i] == ')' )
    {
      Node * newnode = new Node();
      newnode -> right = stk.top();
      newnode -> right -> isRight = 1;
      stk.pop();
      newnode -> left = stk.top();
      newnode -> left -> isRight = 0;
      stk.pop();
      ++i;

      if(((i+1) < (int) str.size()) && (static_cast<int> (str[i+1] - '0') >= 0 && static_cast<int> (str[i+1] - '0') <= 9))
      {

        newnode -> length = getBranchLengths(i+1, str);
      }

      newnode -> time = newnode -> right -> time + newnode -> right -> length;
      stk.push(newnode);

    }

    ++i;

  }

}

void getRanks(Node* newnode, int & tail, Node ** ar)
{
  int k = 1, k1 = 1;
  while(newnode -> time != 0)
  {
    Node * tempnode;
    tempnode = popFromArray(tail, ar);//ar[tail]
    if(tempnode -> rank == 0) tempnode -> rank = k; //added if cond to avoid segfault
    k++;
    k1++;
    if((tail > 0) && (tempnode -> time == ar[tail - 1] -> time) && (tempnode -> time > 0))
    {
      int x;

      //    cout << "Warning: some ranks are tied. Rand assignment was applied" << endl;
      k--;
      srand(5);
      x = rand()%(k1-k+1)+k;
      tempnode -> rank = x;
      if(x == k) ar[tail-1] -> rank = k1;
      else ar[tail-1] -> rank = k;

    }

    else if((tail > 0) && (tempnode -> time != ar[tail - 1] -> time))
    {
      k = k1;
    }

    pushToArray(tempnode -> right, tail, ar);
    pushToArray(tempnode -> left, tail, ar);

    newnode = ar[tail - 1];
  }

}


Node * getNodeFromRank (Node * p, int rankValueGT)
{
  if(p == NULL) cout << "Tree is empty" << endl;
  Node * k;
  Node * t;
  if(p -> label.size() == 0)
  {
    if(p -> rank != rankValueGT)
    {
      t =  getNodeFromRank(p -> right, rankValueGT);
      k =  getNodeFromRank(p -> left, rankValueGT);
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




void getTaxa (Node * p, std::stack <Node *> & allTaxa)
{
  if(p -> label.size() == 0)
  {
    getTaxa(p -> right, allTaxa);
    getTaxa(p -> left, allTaxa);
  }
  else
  {
    allTaxa.push(p);
  }
}


Node * mrcaST (Node * p, Node * rch, Node * lch)
{
  if(p == NULL) return NULL;
  if(p -> label == rch -> label || p -> label == lch -> label)
  {
    return p;
  }
  else
  {
    Node * r = mrcaST(p -> right, rch, lch);
    Node * l = mrcaST(p -> left, rch, lch);

    if(l!=NULL && r!=NULL)
    {
      return p;
    }
    else
    {
      if(l!=NULL)
      {
        return l;
      }
      return r;
    }
  }
}

Node * mrcaST2 (Node * p, Node * tp, std::string v)
{
  if(p == NULL) return NULL;
  if( p == tp || p -> label == v) return p;
  else
  {
    Node * l = mrcaST2(p -> left, tp, v);
    Node * r = mrcaST2(p -> right, tp, v);

    if(l&&r) return p;
    else
    {
      if(l) return l;
      return r;
    }
  }

}



Node * getMrca(Node * p, std::stack <Node * > & stk)
{
  Node * a = stk.top();
  stk.pop();
  Node * b = stk.top();
  stk.pop();
  Node * ptemp = mrcaST(p, a, b);
  while(!stk.empty())
  {
    if(ptemp -> rank == 1) return ptemp;
    Node * t = stk.top();
    stk.pop();
    ptemp = mrcaST2(p, ptemp, t -> label);
  }
  return ptemp;
}

void storeLabels (Node * p, int N, string * ar, int & i)
{
  if(p -> label.size() != 0)
  {
    storeLabels(p -> right, N, ar, i);
    storeLabels(p -> left, N, ar, i);
  }
  else
  {
    ar[i] = p -> label;
    ++i;
  }
}


void addClone(Node * popnode, Queue & q)
{
  for (int i = 0; i < q.size(); ++i)
  {

    Node * inode = q[i];
    Node * parnd = inode -> parent;
    Node * snode = new Node();
    snode -> length = popnode -> length;
    inode -> length -= popnode -> length;
    inode -> parent = snode;
    snode -> parent = parnd;
    snode -> rank = popnode -> rank;
    snode -> distFrRoot = popnode -> distFrRoot;
    snode -> time = popnode -> time;
    if (inode -> isRight == 0)
    {
      snode -> left = inode;
      parnd -> left = snode;
    }
    else
    {
      snode -> right = inode;
      parnd -> right = snode;
    }

  }

}

void makeBeadedTree (Node * root, size_t MaxQ)
{
  Queue q(MaxQ);
  q.push(root);
  while(!q.isEmpty())
  {
    Node * nd = q.pop();
    if(nd -> rank != 1 && nd -> rank != 0) addClone(nd, q);
    if(nd -> left != NULL && nd -> right != NULL)
    {
      q.push(nd -> left);
      q.push(nd -> right);
    }
  }

}


void makeBeadedArray (Node * root, int ** y, size_t MaxQ)
{
  Queue q(MaxQ);
  q.push (root);
  int j = 1;
  int t = 1;
  root -> population = 1;
  root -> outdegree = 2;
  // stk.push(root);
  y[root -> rank][1] = 2;
  while (!q.isEmpty())
  {
    Node * nd = q.pop();
    if (j >= nd -> rank && t != nd -> rank) j = 1;
    if (nd -> rank != 0 && nd -> rank != 1)
    {
      if (nd -> left != NULL && nd -> right != NULL)
      {
        y[nd -> rank][j] = 2;
        nd -> outdegree = 2;
        nd -> population = j;
      }
      else if (nd -> left != NULL || nd -> right != NULL)
      {
        y[nd -> rank][j] = 1;
        nd -> outdegree = 1;
        nd -> population = j;
      }
      j++;
      t = nd -> rank;
    }

    if (nd -> left != NULL) q.push(nd -> left);
    if (nd -> right != NULL) q.push(nd -> right);
  }
}


void getPath (Node * p, std::stack <Node *> & stk)
{
  stk.push(p);
  if(p -> rank != 1)
  {
    p = p -> parent;
    getPath(p, stk);
  }
}

void getPathMatchRank(Node * &  p, int rank_val)
{
  if (p -> rank != rank_val && p -> rank != 1)
  {
    p = p -> parent;
    getPathMatchRank(p, rank_val);

  }

}


int * getNextHistory (int * ar, int * ar_orig, int n, int idx, int & numb)
{
  if (*(ar + n - 2) == 1) return ar;
  for(int t = n - 2; t >=0; t--)
  {
    if(*(ar + t - 1) < *(ar + t))
    {
      *(ar + t) -= 1;
      for(int j = t + 1; j <= n - 2; ++j)
      {
        *(ar + j) = *(ar_orig + j);
      }
      numb++;
      if (numb == idx) return ar;
      else
      {
        ar = getNextHistory(ar, ar_orig, n, idx, numb);
      }
    }
  }
  return ar;
}



void arrayFreq(int * arr, int * res, int n)
{
  for (int j =0; j<n; j++)
    arr[j] = arr[j]-1;
  for (int i=0; i<n; i++)
    arr[arr[i]%n] = arr[arr[i]%n] + n;
  for (int i =0; i<n; i++)
    res[i] = arr[i]/n;
}


void get_hseq(int * arseq, Node * nodeGT, Node * nodeST, int n)
{
  arseq[0] = 1;
  for(int i = 2; i < n; ++i)
  {
    std::stack <Node *> allTaxa;
    getTaxa(getNodeFromRank(nodeGT, i), allTaxa);
    Node * mrcanode  = getMrca(nodeST, allTaxa);
    if (mrcanode -> rank < i)
    {
      arseq[i-1] = mrcanode -> rank;
    }
    else arseq[i-1] = i;
  }
}

void reverseArraySort(int * arseq, int n)
{
  for(int i = n - 2; i > 0; --i)
  {
    if(arseq[i-1] > arseq[i])
      arseq[i-1] = arseq[i];
  }
}



void getDescTaxa(Node * p, int N)
{
  if(p -> label.size() == 0)
  {
    stack <Node *> allTaxa;
    getTaxa(p, allTaxa);
    string temp;
    p -> leavesnum = allTaxa.size();
    while(!allTaxa.empty())
    {
      temp += allTaxa.top() -> label;
      allTaxa.pop();
    }
    p -> desctaxa = temp;

    getDescTaxa(p -> right, N);
    getDescTaxa(p -> left, N);
  }
  else
  {
    p -> desctaxa = p -> label;
  }
}


void returnMrca(Node * nodeST, Node * nodeGT, int n, Node ** ar)
{
  for(int i = 1; i < n; ++i)
  {
    std::stack <Node *> allTaxa;
    getTaxa(getNodeFromRank(nodeGT, i), allTaxa);
    Node * mrcanode  = getMrca(nodeST, allTaxa);
    ar[i-1] = mrcanode;
  }
}

void get_node(int * rank_hist, int ** temp, int n, Node ** ar)
{
  int tmp = 0;
  int j = 0;
  Node * mrcanode;

  for(int i = 1; i < n; ++i)
  {
    mrcanode = ar[i-1];
    getPathMatchRank(mrcanode, rank_hist[i-1]);
    if (tmp != mrcanode -> rank) j = 1;

    temp[mrcanode->rank-1][j-1] = mrcanode -> population;
    tmp = mrcanode -> rank;
    ++j;
  }
}

void calc_k_ijz (int ***k, int n, int * m, int ** y, int ** ar_r/*, int*  ar_j*/)
{
  for(int i = 1; i <= n; ++i)
  {
    k[n][0][i] = 1;

  }
  for(int i = n - 1; i >= 1; --i)
  {
    int tempz = 1;
    for(int z = 1; z <= i; ++ z)
    {

      for(int t = 1; t <= y[i][z]; ++t) //outdegree
      {
        k[i][m[i-1]][z] += k[i+1][0][tempz];
        ++tempz;

      }

      for(int j = m[i-1] - 1; j >= 0; j--)
      {
        for(int h = j; h >=0; h--)
        {
          if (ar_r[i-1][h] == z)
            k[i][h][z] = k[i][h+1][z] - 1;
          else
            k[i][h][z] = k[i][h+1][z];

        }

      }
    }
  }
}


void writeTree(Node* root)
{
  cout << "Name " << root->debug_name;
  if (root->left) cout << " left child: " << root->left->debug_name;
  if (root->right) cout << " right child: " << root->right->debug_name;
  if (root->parent) cout << " parent " << root->parent->debug_name;
  cout << " rank " << root->rank;
  cout << " distFrRoot " << root->distFrRoot;
  cout << " label " << root -> label;
  cout << " desctaxa " << root -> desctaxa;
  cout << endl;
  //	cout << " time " << root -> time;
  //	cout << " population " << root -> population << endl;
  //	cout << " outdegree " << root -> outdegree << endl;


  if (root->left) writeTree(root->left);
  if (root->right) writeTree(root->right);
}

double factorial (int num)
{
  double accum = 1.;
  while (num > 0)
  {
    accum *= num;
    num--;
  }
  return accum;
}


//change bin_coef to int
double calcBinomCoef (int n, int k)
{
  double prod = 1.;
  for(int i = 1; i <= k; ++i)
    prod *= (double) (n + 1 - i)/i;
  return prod;
}


void writeRankedTreeTopology(Node * p, int N, ofstream & file)
{
  string * tempstr = new string[N];

  for(int i = 2; i < N; ++i)
  {
    stack <Node *> allTaxa;
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
      file << tempstr[j];
    file << "-" << i << "-";

    for(int t = 0; t < N; ++t)
      tempstr[t] = "";
  }
  file << "\n";
  delete [] tempstr;
}


void storeLabels(int n, int & i, string* str, Node * p)
{
  if(p -> label.size() == 0)
  {
    storeLabels(n, i, str, p -> right);
    storeLabels(n, i, str, p -> left);
  }
  else
  {
    str[i] = p -> label;
    ++i;
  }
}

void getS(Node * p, double * s)
{
  if(p -> label.size() == 0)
  {
    s[p -> rank - 1] = p -> time;
    getS(p -> right, s);
    getS(p -> left, s);
  }
}
/*
   int g_i(int N, int i, Node * nodeST, Node * nodeGT)
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
   mrcanode  = getMrca(nodeST, allTaxa);
   if(mrcanode -> rank <= i) prod = 0;
   }
   res += prod;
   }
   return N - res;
   }
   */


double lambda_ij (int i, int j, int *** k)
{
  double res = 0.;
  for(int z = 1; z < i + 1; ++z)
  {
    res += calcBinomCoef(k[i][j][z], 2);
  }
  return res;
}

double invPartialCoal(int j)
{
  return (double) (pow(2, j-1) / factorial(j)) / factorial(j-1);
}

double conditionalProbability(int i, int * m_i, int *** array_k, double * s) // min i = 2 since l[i-2] exists
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
    res += (double) exp(-lambda_ij(i, j, array_k)*s[i-2])/prod;
  }

  return res;
}




void getInvPartialCoal(double * ar, int N)
{
  for(int i = 2; i < N + 1; ++ i)
  {
    ar[i-2] = (double) (pow(2, i - 1) / factorial(i)) / factorial(i - 1);
  }
}


double geneTreeProbability(int * m, int *** k, double * s, double * coal, int n)
{
  double res = 1.;
  for(int i = 2; i < n; ++i)
  {
    res *= conditionalProbability(i, m, k, s);
  }
  return res*coal[m[0] - 1];
}



double getGeneTreeProb(int N, double * s,  Node * newnode, Node * newnodeGT, int ** ar_y, double * array_invcoal, Node ** arMrca, int ** ar_rankns, int ***k)
{

  double arhists[MAX_EXPONENT_VALUE];
  for(int i = 0; i < MAX_EXPONENT_VALUE; ++i)
    arhists[i] = 0.;

  returnMrca(newnode, newnodeGT, N, arMrca);

  int * arRankHistory = new int [N-1];
  get_hseq(arRankHistory, newnodeGT, newnode, N);
  reverseArraySort(arRankHistory, N);
  //   cout << "Max Ranked History" << endl;
  //  for(int i = 0; i < N-1; ++i) cout << arRankHistory[i] << " ";
  //  cout << endl;

  int *  arRankHistoryCopy = new int[N-1];
  for(int i = 0; i < N-1; ++i)
  {
    arRankHistoryCopy[i] = arRankHistory[i];
  }
  int * m_i = new int[N-1];
  arrayFreq(arRankHistoryCopy, m_i, N-1);//1st arg changes after running this func
  /*    cout << "m_i[i]: " << endl;
        for(int i=0; i< N-1; ++i)
        cout << m_i[i] << " ";
        cout << endl;
        cout <<"rank_hist " << endl;
        for(int i = 0; i < N-1; ++i)
        {
        cout  << arRankHistory[i] << " ";
        }
        cout << endl; */
  int max_m_i = *max_element(m_i, m_i+N-1);
  get_node(arRankHistory, ar_rankns, N, arMrca);
  for(int i = 0; i < N+1; i++)
    for(int j = 0; j < N; ++j)
      for(int z = 0; z < N+1 ; ++z)
      {
        k[i][j][z] = 0;
      }
  calc_k_ijz(k, N, m_i, ar_y, ar_rankns);//*, m_i);
  double onehistprob = geneTreeProbability(m_i, k, s, array_invcoal, N);
  double prob_val = onehistprob;

  if(onehistprob != 0 && ((int) floor(fabs(log2(fabs(onehistprob))))) < MAX_EXPONENT_VALUE)
    arhists[ (int) floor(fabs(log2(fabs(onehistprob))))] +=  fabs(onehistprob);
  //  cout<< floor(fabs(log2(fabs(onehistprob)))) << endl;
  int * arRHcopy = new int[N-1];
  for(int i = 0; i < N - 1; ++i)
    arRHcopy[i] = arRankHistory[i];

  int history_counter = 1;
  int numb = 0;
  int * next_history;

  while(*(arRHcopy + N - 2) != 1)
  {
    next_history = getNextHistory (arRHcopy, arRankHistory, N, history_counter, numb);

    /*             cout << "---------- Next Label history ---------- " << endl;
                   for(int i = 0; i < N - 1; ++i)
                   cout << next_history[i] << " ";
                   cout << endl;
                   */
    //repeat
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
    calc_k_ijz(k, N, m_i, ar_y, ar_rankns);//*, m_i);
    onehistprob = geneTreeProbability(m_i, k, s, array_invcoal, N);
    prob_val += onehistprob;

    if(onehistprob != 0 && ((int) floor(fabs(log2(fabs(onehistprob))))) < MAX_EXPONENT_VALUE)
      arhists[ (int) floor(fabs(log2(fabs(onehistprob))))] +=  fabs(onehistprob);
    //     cout << "get_prob = " << geneTreeProbability(m_i, k, s, N)  << " " << log(geneTreeProbability(m_i, k, s, N)) <<  endl;
    //   cout << "prob_of_history = " << setprecision(10) << prob_val << " " << log(prob_val) << endl;
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

int getNumberOfTaxa(int & arg_counter, char* argv[], Node*& newnode)
{
  std::string str="";
  ifstream fin(argv[arg_counter]);
  std::getline(fin, str, '\n');
  fin.close();
  ++arg_counter;
  removeSpaces(str);
  std::stack <Node *> stkST;
  int lbl = 0;
  double tempdist = 0.;
  countParentheses(str);

  pushNodes(lbl, stkST , str);
  newnode = stkST.top();
  newnode -> distFrRoot = 0.;

  distFromRoot(newnode);

  isUltrametric(newnode, tempdist);
  deleteStack(stkST);
  return lbl;
}

void speciesTreeProcessing(Node* newnode, int & N, double* s_times, double * s, int** ar_y)
{
  Node ** ar = new Node * [N];
  int tail = 0;
  pushToArray(newnode, tail, ar);

  getRanks(newnode, tail, ar);

  getS(newnode, s_times);
  for(int j = 2; j < N; ++j)
  {
    s[j-2] = s_times[j-2] - s_times[j-1];
  }
  int itemp = 0;
  double * arDistFrRoot = new double [N];
  saveDistFrRoot(newnode, arDistFrRoot, itemp);
  //printIsUltrametric(arDistFrRoot, N);

  //    ofstream stfile("STtopo.txt");
  //    writeRankedTreeTopology(newnode, N, stfile);
  //    stfile.close();

  size_t maxQue = (size_t) N*(1+N)/2;
  makeBeadedTree(newnode, maxQue);
  makeBeadedArray(newnode, ar_y, maxQue);

  delete[] ar;
  delete[] arDistFrRoot;
}




double calcRankedProb(int & arg_counter, char* argv[], int & N, Node * newnode, double* s, int** ar_y)
{
  string strGT="";
  string strTops="";
  ifstream finGT(argv[arg_counter]); //gtuniqtrees.txt
  ++arg_counter;
  ifstream fintops(argv[arg_counter]); //gtuniqtops.txt
  ofstream finprobGT("probForEachGT.txt");
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
  while(getline(finGT, strGT, ';'))
  {
    getline(fintops, strTops, '\n');
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
    one_gt_prob = getGeneTreeProb(N, s,  newnode, newnodeGT, ar_y, array_invcoal, arMrca, ar_rankns, k);
    finprobGT  << one_gt_prob <<  '\t' << strTops << '\n';
    prec_total_prob += one_gt_prob;

    delete [] arGT;

    deleteTree(newnodeGT);
    deleteStack(stkGT);

  }

  //printf("TOTAL: %.*f\n", 15, prec_total_prob);
  //cout << setprecision(15)  << prec_total_prob << endl;

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
  finprobGT.close();
  fintops.close();

  return prec_total_prob;

}



void calcProbsRankedGtInput(int &  arg_counter, char* argv[])
{
  Node * newnode;
  double res = 0.;
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

  res = calcRankedProb(arg_counter, argv, N, newnode, s, ar_y);
  cout << "Total: " << res << endl;
  for(int i = 0; i < N; ++i)
  {
    delete[] ar_y[i];
  }
  delete[] ar_y;
  delete[] s_times;
  delete[] s;
  deleteTree(newnode);
}



