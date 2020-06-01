#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cmath>
#include <algorithm>
#include "string_manipulations.h"
#include "probs_calculation.h"
#include "rank_distance.h"

using namespace std;

void writeTopoToMap(unordered_map <string, int> & clade_and_rank, string str)
{

	Node * newnode;
	removeSpaces(str);
	std::stack <Node *> stkST;
	int N = 0;
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
		getTaxa(getNodeFromRank(newnode, i), allTaxa);
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

		pair <string, int> temp (combinedstr, i);
		clade_and_rank.insert(temp);

		for(int t = 0; t < N; ++t)
			tempstr[t] = "";
	}
	delete [] tempstr;
	deleteTree(newnode);
	delete[] ar;
}



void getRankDistance(int &  arg_counter, char* argv[])
{
	string str_first = "";
	string str_second = "";
	ifstream sp_file_first(argv[arg_counter]);
	arg_counter++;
	ifstream sp_file_second(argv[arg_counter]);
	ofstream dist_file("outRankDist.txt");



	while(getline(sp_file_first, str_first, ';'))
	{
		getline(sp_file_second, str_second, ';');
    		if(str_first.size() < 3 || str_second.size() < 3) break;
		int N = 0;
		int sum_res = 0;
		unordered_map <string, int> ar_first;
		writeTopoToMap(ar_first, str_first);
		unordered_map <string, int> ar_second;
		writeTopoToMap(ar_second, str_second);

		//    for (auto it = ar_first.begin(); it != ar_first.end(); ++it) 
		//        cout << it->first << " " << it->second  << endl;

		//    for (auto it = ar_second.begin(); it != ar_second.end(); ++it) 
		//        cout << it->first << " " << it->second  << endl;


		for (auto it = ar_first.begin(); it != ar_first.end(); ++it) 
		{
			if(ar_second.find(it->first) != ar_second.end())
			{
				sum_res += pow(ar_second.find(it->first) -> second - it->second, 2);
			}
			else
			{
				sum_res = 357;
			}
		}
		dist_file << sum_res << endl;
	}
	sp_file_first.close();
	sp_file_second.close();
        dist_file.close();

}






