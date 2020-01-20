#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cmath>
#include <algorithm>
#include "string_manipulations.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"
#include "get_unranked_topology.h"

using namespace std;
void getClades (Node * p, std::stack <string> & allTaxa)
{
    if(p -> label.size() == 0)
    { 
        getClades(p -> right, allTaxa);
        getClades(p -> left, allTaxa);
        allTaxa.push(p -> desctaxa);
    }
}

int countNumberTaxa(int & indicator, string & str)
{
    int i = 0;
    int lbl = 0;
    while(i < (int) str.size())
    {
        if((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= 'A' && str[i] <= 'Z'))
        {
            while(((i+1) < (int) str.size()) && str[i+1] != '(' && str[i+1] != ',' && str[i+1] != ')')
            {
                ++i;
            }
            ++lbl;
        }

        if(str[i] == ':') ++indicator;
        ++i;
    }

    return lbl;
}



void finSort(int n, string* str)
{
    string temp="";
    string hash="|";
    for(int i = 0; i < n-1; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            if(str[i].size() == str[j].size())
            {
                int t = 0;
                string itemp = "" ;
                string jtemp = "" ;

                while(hash.compare(&str[i][t]) != 0)
                {
                    itemp += str[i][t];
                    jtemp += str[j][t];
                    ++t;
                }

                if(itemp > jtemp)
                {
                    temp = str[i];
                    str[i] = str[j];
                    str[j] = temp; 

                }
            }
        }
    }   
}


void getUnrDescTaxa(Node * p, int N)
{
    string hash = "|";
    if(p -> label.size() == 0)
    {
        stack <Node *> allTaxa;
        getTaxa(p, allTaxa);

        int tmpcount = 0;
        int cur_size = allTaxa.size();
        string tempstr[cur_size];
        string temp;
        while(!allTaxa.empty())
        {
            tempstr[tmpcount] = allTaxa.top() -> label;
            temp += allTaxa.top() -> label;
            allTaxa.pop();
            tmpcount++;  
        } 
        sortString(cur_size, tempstr);
        temp = "";

        for(int j = 0; j < cur_size; ++j)
        {
            temp +=  tempstr[j];
            temp += hash;
        }
        p -> desctaxa = temp;
        getUnrDescTaxa(p -> right, N);
        getUnrDescTaxa(p -> left, N);
    }
    else
    {
        p -> desctaxa = p -> label;
    }
}








void outputUnrankedTopology(int &  arg_counter, char* argv[])
{
    string str;
    string arg;
    ifstream gt_file(argv[arg_counter]);

    getline(gt_file, str, ';');
    int indicator = 0;
    int lbl = countNumberTaxa(indicator, str);
    int N = lbl;
    int dummy = 0;
    gt_file.close();

    gt_file.open(argv[arg_counter]);
    ofstream topo_file("outUnrTopos.txt");

    if(indicator == 0)
    {
        while(getline(gt_file, str, ';'))
        {
            if(str.size() < 3) break;
            removeSpaces(str);
            std::stack <Node *> stkST;
            string ar_strlbl[N-1];
            pushNodesUnrankedGT(dummy, stkST, str, ar_strlbl);
            Node * newnode = stkST.top();
            countParentheses(str); 

            stack <string> allTaxa;
            getUnrDescTaxa(newnode, N);
            getClades(newnode, allTaxa);

            int tmpcount = 0;
            string tempstr[N];

            while(!allTaxa.empty())
            {
                tempstr[tmpcount] = allTaxa.top();
                allTaxa.pop();
                tmpcount++;  
            }

            sortStringsBySize(tempstr, N-1);
            finSort(N-1, tempstr);

            for(int j = 0; j < N-1; ++j)
            {
                topo_file << tempstr[j] << "-";
            }

            topo_file << endl;
            deleteStack(stkST);
            deleteTree(newnode);
        }
    }
    else
    {
        while(getline(gt_file, str, ';'))
        {
            if(str.size() < 3) break;

            removeSpaces(str);
            stack <Node *> stkST;
            countParentheses(str); 
            pushNodes(lbl, stkST, str);
            Node * newnode = stkST.top();

            stack <string> allTaxa;
            getUnrDescTaxa(newnode, N);
            getClades(newnode, allTaxa);

            int tmpcount = 0;
            string tempstr[N];

            while(!allTaxa.empty())
            {
                tempstr[tmpcount] = allTaxa.top();
                allTaxa.pop();
                tmpcount++;  
            }

            sortStringsBySize(tempstr, N-1);
            finSort(N-1, tempstr);

            for(int j = 0; j < N-1; ++j)
            {
                topo_file << tempstr[j] << "-";
            }

            topo_file << endl;
            deleteStack(stkST);
            deleteTree(newnode);
        }
    }
    ofstream freq_file("outUnrFreqs.txt");
    topo_file.close();

    ifstream infile;
    infile.open("outUnrTopos.txt");
    if (infile.is_open()) countUniqueStrings(infile, freq_file);
    freq_file.close();
    gt_file.close();
    infile.close();

}






/*Jan13
  void outputUnrankedTopology(int &  arg_counter, char* argv[])
  {
  string str;
  string arg;
  ifstream gt_file(argv[arg_counter]);

  getline(gt_file, str, ';');
  int indicator = 0;
  int lbl = countNumberTaxa(indicator, str);
  int N = lbl;
  int dummy = 0;
  gt_file.close();

  gt_file.open(argv[arg_counter]);
  ofstream topo_file("utopos.txt");

  if(indicator == 0)
  {
  while(getline(gt_file, str, ';'))
  {
  if(str.size() < 3) break;
  removeSpaces(str);
  std::stack <Node *> stkST;
  string ar_strlbl[N-1];
  pushNodesUnrankedGT(dummy, stkST, str, ar_strlbl);
  Node * newnode = stkST.top();
  countParentheses(str); 

  stack <string> allTaxa;
  getUnrDescTaxa(newnode, N);
  getClades(newnode, allTaxa);

  int tmpcount = 0;
  string tempstr[N];

  while(!allTaxa.empty())
  {
  tempstr[tmpcount] = allTaxa.top();
  allTaxa.pop();
  tmpcount++;  
  }

  sortStringsBySize(tempstr, N-1);
  finSort(N-1, tempstr);

  for(int j = 0; j < N-1; ++j)
  {
  topo_file << tempstr[j] << "-";
  }
  topo_file << endl;
  deleteStack(stkST);
  deleteTree(newnode);
  }
  }
  else
  {
  while(getline(gt_file, str, ';'))
  {
  if(str.size() < 3) break;

  removeSpaces(str);
  stack <Node *> stkST;
  countParentheses(str); 
  pushNodes(lbl, stkST, str);
  Node * newnode = stkST.top();

  stack <string> allTaxa;
  getUnrDescTaxa(newnode, N);
  getClades(newnode, allTaxa);

int tmpcount = 0;
string tempstr[N];

while(!allTaxa.empty())
{
    tempstr[tmpcount] = allTaxa.top();
    allTaxa.pop();
    tmpcount++;  
}

sortStringsBySize(tempstr, N-1);
finSort(N-1, tempstr);

for(int j = 0; j < N-1; ++j)
{
    topo_file << tempstr[j] << "-";
}
topo_file << endl;
deleteStack(stkST);
deleteTree(newnode);
}
}
ofstream freq_file("ufreqs.txt");
topo_file.close();

ifstream infile;
infile.open("utopos.txt");
if (infile.is_open()) countUniqueStrings(infile, freq_file);
freq_file.close();
gt_file.close();
infile.close();

}
*/

