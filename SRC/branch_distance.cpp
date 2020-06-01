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
#include "probs_calculation.h"
#include "branch_distance.h"


using namespace std;

void outputCoalIntervals(int &  arg_counter, char* argv[])
{
  Node * newnode;
  double res = 0.;
  int N = getNumberOfTaxa(arg_counter, argv, newnode);

  double * s_times = new double [N-1];
  double * s = new double [N-2];
  vector <Node *> vlabels;
  int ** ar_y = new int * [N];
  for (int i = 0; i < N; i++) ar_y[i] = new int [N];

  for (int i = 0; i < N; ++i)
    for (int j = 0; j < N; ++j)
    {
      ar_y[i][j] = 0;
    }

  Node ** ar = new Node * [N];
  int tail = 0;
  pushToArray(newnode, tail, ar);
  getRanks(newnode, tail, ar);
 
  if(vlabels.size() != 0)
  {
    while (!vlabels.empty())
    {
         vlabels.pop_back();
    }
  }

  getS(newnode, s_times, vlabels);
 
  for(int j = 2; j < N; ++j)
  {
    s[j-2] = s_times[j-2] - s_times[j-1];
    cout << setprecision(6) << s[j-2] << ' ';
  }
  cout << endl;
  delete[] ar;

  for(int i = 0; i < N; ++i)
  {
    delete[] ar_y[i];
  }
  delete[] ar_y;
  delete[] s_times;
  delete[] s;
  deleteTree(newnode);
}


void writeBranchLengths(unordered_map <string, double> & clade_and_rank, string str, int & N)
{
	Node * newnode;
	removeSpaces(str);
	std::stack <Node *> stkST;
	double tempdist = 0.;
	countParentheses(str);
	pushNodes(N, stkST , str);
	newnode = stkST.top();
	newnode -> distFrRoot = 0.;

	distFromRoot(newnode);
	isUltrametric(newnode, tempdist);
	deleteStack(stkST);

	Node ** ar = new Node * [N];
	int tail = 0;
	pushToArray(newnode, tail, ar);
	getRanks(newnode, tail, ar);

	string * tempstr = new string[N];
	for(int i = 2; i < N; ++i)
	{
		stack <Node *> allTaxa;
		int tmpcount = 0;
        Node * cur_node = getNodeFromRank(newnode, i);
		getTaxa(cur_node, allTaxa);
		while(!allTaxa.empty())
		{
			tempstr[tmpcount] = allTaxa.top() -> label;
			allTaxa.pop();
			tmpcount++;
		}
		sortString(N, tempstr);

		string combinedstr = "";    
		for(int j = 0; j < N; ++j)
			combinedstr += tempstr[j];

		pair <string, double> temp (combinedstr, cur_node -> length);
		clade_and_rank.insert(temp);

		for(int t = 0; t < N; ++t)
			tempstr[t] = "";
	}
	delete [] tempstr;
	deleteTree(newnode);
	delete[] ar;
}



void getMSEBranchDistance(int &  arg_counter, char* argv[])
{
	string str_first = "";
	string str_second = "";
	ifstream sp_file_first(argv[arg_counter]);
	arg_counter++;
	ifstream sp_file_second(argv[arg_counter]);
	ofstream dist_file("outBranchDist.txt");



	while(getline(sp_file_first, str_first, ';'))
	{
		getline(sp_file_second, str_second, ';');
    		if(str_first.size() < 3 || str_second.size() < 3) break;
		int N = 0;
		double sum_res = 0.;
		unordered_map <string, double> ar_first;
		writeBranchLengths(ar_first, str_first, N);
		unordered_map <string, double> ar_second;
		N = 0;
		writeBranchLengths(ar_second, str_second, N);

	/*	    for (auto it = ar_first.begin(); it != ar_first.end(); ++it) 
		        cout << it->first << " " << it->second  << endl;
            cout << endl;
		    
            for (auto it = ar_second.begin(); it != ar_second.end(); ++it) 
		        cout << it->first << " " << it->second  << endl;
            cout << "--------------- " << endl;
      */       
		for (auto it = ar_first.begin(); it != ar_first.end(); ++it) 
		{
		    cout << it->first << " " << it -> second << " ";
            
			if(ar_second.find(it->first) != ar_second.end())
			{
				sum_res += pow(ar_second.find(it->first) -> second - it->second, 2);
				cout << ar_second.find(it -> first) -> second << " "  << fabs(ar_second.find(it->first) -> second - it->second) << endl;
			}
			else
			{
				sum_res = 1000.;
			}
		}
		dist_file << (double) sum_res/(N-2) << endl;
       // cout << endl;
       // dist_file << endl;
		//cout << (double) sum_res/(N-2) << endl;
	}

	sp_file_first.close();
	sp_file_second.close();
    dist_file.close();

}






