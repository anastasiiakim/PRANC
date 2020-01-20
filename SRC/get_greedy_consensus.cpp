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
#include "write_ranked_tree.h"
#include "get_greedy_consensus.h"
#include <sstream>

using namespace std;

string idToStr(int id)
{
    stringstream ss;
    ss << 't' << id;
    string str;
    ss >> str;
    return str;
}

int strToId(string str)
{
    cout << str << ' ';
    stringstream ss;
    str[0] = ' ';
    ss << str;
    int id;
    ss >> id;
    cout << id << ' ';
    cout << idToStr(id) << endl;
    return id;
}
/*
struct BagOfNodes
{
    vector<Node*> m;
    // Constructor

    BagOfNodes(){};
    BagOfNodes(vector<string> str)
    {
        for(unsigned i = 0; i < str.size(); ++i)
        {
            m.push_back(new Node(str[i]));
        }
    };
    BagOfNodes(vector<string> str, const BagOfNodes & b)
    {
        for(unsigned i = 0; i < str.size(); ++i)
        {
            for(unsigned j = 0; j < b.m.size(); ++j)
            {
                if (b.m[j]->label == str[i])
                {
                    m.push_back(b.m[j]);
                    break;
                } 
            }
        }
    };
    void clip(Node * rule)
    {
        auto it = m.end();
        it = find(m.begin(), m.end(), rule->left);
        if (it != m.end()) m.erase(it);
        it = find(m.begin(), m.end(), rule->right);
        if (it != m.end())
        {
            m.erase(it);
            m.push_back(rule);
        }
    };
};
*/

Node* makeTree(vector<vector<string>> strs)
{  
    BagOfNodes tmp[strs.size()];
    tmp[strs.size()-1] = BagOfNodes(strs[strs.size()-1]);
    for(int i = 0; i < strs.size()-1; ++i)
    {
        tmp[i] = BagOfNodes(strs[i], tmp[strs.size()-1]);
    }
    for(int i = 0; i < strs.size()-1; ++i)
    {
        if(tmp[i].m.size() > 2) cout << "ERROR!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        Node* p = new Node(tmp[i].m[0], tmp[i].m[1]);
        for(int j = i+1; j < strs.size(); ++j)
        {
            tmp[j].clip(p);
        }
    }

    int i = strs.size()-1;
    if(tmp[i].m.size() > 2) cout << "ERROR!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
    return new Node(tmp[i].m[0], tmp[i].m[1]);
}


inline bool cmpt(const int & a, const int & b)
{
    return ((a | b) == b) || ((a | b) == a) || !(a & b);
}

bool sortMapByVal(const pair<string, int> &a, const pair<string, int> &b) 
{ 
    return a.second > b.second; 
} 

bool sortVectBySize(const pair<string, int> & a, const pair<string, int> &b) //ToDo 
{
    return a.first.size() > b.first.size();
}






ret_t visit(const int & row, const int & col, const unsigned long & freq, unsigned nsize, unsigned long a[][2], vector<int> path)
{
    // printf("%d  %d  %3lu \n", row, col, freq);
    ret_t ret1;//r
    ret_t ret2;//d
    if(col < nsize-1)
    {
        ret1 = visit(row, col+1, freq, nsize, a, path);
        bool cmp = true;
        for (unsigned i = 0; i < path.size(); ++i)
            if(!cmpt(a[path[i]][0], a[col][0]))
            {
                cmp = false;
                break;
            }
        if (cmp)
        {
            path.push_back(col);
            ret2 = visit(col, col+1, freq + a[col][1], nsize, a, path);

            if (ret1.first < ret2.first) 
            {
                return ret2;
            }
        }

        return ret1;
    }
    else
    {
        return ret_t(freq, path);
    }
}







void getConsensusTree(int &  arg_counter, char* argv[])
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
    ofstream topo_file("outGreedyCons.txt");

    std::unordered_map <string, int> clpairs;
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


            int itemp = 0;
            for(int j = 0; j < N-1; ++j)
            {
                string str = tempstr[j];
                if(clpairs.count(str))
                {
                    clpairs[str]++;
                }  
                else
                {
                    pair <string, int> temp (str, 1);
                    clpairs.insert(temp);
                }

            }


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
            int itemp = 0;
            for(int j = 0; j < N-1; ++j)
            {
                string str = tempstr[j];
                if(clpairs.count(str))
                {
                    clpairs[str]++;
                }  
                else
                {
                    pair <string, int> temp (str, 1);
                    clpairs.insert(temp);
                }


            }
            deleteStack(stkST);
            deleteTree(newnode);
        }
    }
    gt_file.close();


    vector<pair<string, int>> vec;
    std::map<string, int> ordered(clpairs.begin(), clpairs.end());
    for(auto it = ordered.begin(); it != ordered.end(); ++it)
        vec.push_back(make_pair(it->first, it->second));

    sort(vec.begin(), vec.end(), sortMapByVal); 
    vector <int> uniq;
    uniq.push_back(vec[0].second);
    int uj = 0;

    for (int i = 0; i < vec.size(); i++)
    {
        // cout << vec[i].first << ": " << vec[i].second << endl;
        // cout << vec[i].first.size() << endl;
        if(vec[i].second != uniq[uj]) 
        {
            uniq.push_back(vec[i].second);
            ++uj;
        }
    }

    int j = 0;
    int updstart = 0;
    for(int t = 0; t < uniq.size(); ++t)
    {
        while( (j < vec.size()) && (vec[j].second == uniq[t]) )
        {
            ++j;
        }
        sort(vec.begin()+updstart, vec.begin()+j, sortVectBySize);
        updstart = j;
    }

    //vector of sorted pairs by their occurence
/* 
    for(int i = 0; i < vec.size(); i++)
    {
        cout << vec[i].first << ": " << vec[i].second << endl;
    }
*/

    //search for a tree with the highest score
    unsigned long int a [vec.size()][2];
    for(int i = 0; i < vec.size(); ++i)
    {
        a[i][0] = 0;
        a[i][1] = vec[i].second;

        stringstream ss;
        int id;
        for(int j = 0; j < vec[i].first.size(); ++j)
        {
            char c = vec[i].first[j];
            if(c == 't') continue;
            if(c == '|')
            {
                ss >> id;
                ss.str("");
                ss.clear();
                a[i][0] = a[i][0] | (1 << id);
            } 
            else
            { 
                ss << c - 'A'; //if A,B,C,..not t
            }   
           
//            if(vec[i].first[j] == '|') continue;
//            a[i][0] = a[i][0] | (1 << (int) (vec[i].first[j] - 'A'));

        }
    }

    bool b[vec.size()][vec.size()];
    for(int i = 0; i < vec.size(); ++i)
    {
        for(int j = 0; j < vec.size(); ++j)
        { 
            b[i][j] = cmpt(a[i][0],a[j][0]);
      //      cout << b[i][j] << ' ';
        }
      //  cout << " freq: " << a[i][1] << endl;
    }

    int row = 0;
    int col = 1;
    unsigned long int freq = a[0][1];
    vector<int> path;
    path.push_back(0);
    ret_t visit_ret =  visit(row, col, freq, vec.size(), a, path);

    string m[visit_ret.second.size()];
  //  printf("Freq: %lu \n", visit_ret.first);
    for (unsigned i = 0; i < visit_ret.second.size(); ++i) 
    {
        m[i] = vec[visit_ret.second[i]].first.c_str();
 //       printf("%s %d\n",vec[visit_ret.second[i]].first.c_str(),vec[visit_ret.second[i]].second);
    }


    string temp = "";
    sortStringsBySize(m, visit_ret.second.size());
    unordered_map <int, int> mfreqs;
    for (int i = 0; i < visit_ret.second.size(); ++i) 
    {
        int j = 0;
        int jcount = 0;
        string msep = "|";
        string mtemp = m[i];
        while(j < (int) mtemp.size())
        {
            if(mtemp[j] == '|')
            {
                ++jcount;
            }

            ++j;
        }
        if(mfreqs.count(jcount))
        {
            mfreqs[jcount]++;
        }  
        else
        {
            pair <int, int> temp (jcount, 1);
            mfreqs.insert(temp);
        }

    }
    /*
       stack <int> mstack;
       for (auto it = mfreqs.begin(); it != mfreqs.end(); ++it) 
       {
       cout << it->first << " " << it->second  << endl;
       mstack.push(it->second);
       }

*/


    vector<vector<string>> strs;
    for (int i = 0; i < visit_ret.second.size(); ++i) 
    {
        strs.push_back(vector<string>());
        temp = m[i];
        int j = 0;

        while(j < (int) temp.size())
        {
            if(((temp[j] >= 'a' && temp[j] <= 'z') || (temp[j] >= 'A' && temp[j] <= 'Z')) && ((j+1) < (int) temp.size()))
            {

                string strtemp;
                strtemp += temp[j];
                while(((j+1) < (int) temp.size()) &&  temp[j+1] != '|')
                {
                    strtemp += temp[j+1];
                    ++j;
                }
          //      cout << m[i] << " ------ " << strtemp << endl;
                strs[i].push_back(strtemp);
                ++j;

            }
            else
            {
                ++j;
            }

        }

    }

 /*   for (int i = 0; i < strs.size(); ++i) 
    {
        for (int j = 0; j < strs[i].size(); ++j)
        {
            cout << strs[i][j] << " "; 
        } 
    cout << endl;
    
    }
    */

    Node * tree = makeTree(strs);
  //  writeTree(tree);
    string finstr = "";
    writeUnrankTreeIntoStr(tree, N, finstr);
    deleteTree(tree);
    topo_file << finstr << ";" << endl;
    topo_file.close();
}



