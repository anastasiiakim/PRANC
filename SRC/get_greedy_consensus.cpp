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
#include <chrono> 

using namespace std;
using namespace chrono; 

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

Node* makeTree(vector<vector<string>> strs, int & correct_tree)
{ 
    BagOfNodes tmp[strs.size()];
    tmp[strs.size()-1] = BagOfNodes(strs[strs.size()-1]);
    vector<Node*> history = tmp[strs.size()-1].m;

    for(int i = 0; i < strs.size()-1; ++i)
    {
        tmp[i] = BagOfNodes(strs[i], tmp[strs.size()-1]);
    }

    /*  
        for(int i = 0; i < strs.size()-1; ++i) 
        {
        for (int j = 0; j < tmp[i].m.size(); ++j)
        {
        cout<< tmp[i].m[j]->label << " ";
        }
        cout << endl;
        }

*/

    for(int i = 0; i < strs.size()-1; ++i)
    {
        // clip only two nodes, because tree is binary
        if(tmp[i].m.size() > 2) 
        {
            cout << i << "--ERROR!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
            correct_tree = 0;
            // for(int j = 0; j < tmp[i].m.size(); ++j) 
            //     cout << tmp[i].m[j]->label << endl;
            for(auto it = history.begin(); it != history.end(); ++it)
            {
                delete *it;
            }
            return NULL;
        }
        Node* p = new Node(tmp[i].m[0], tmp[i].m[1]);
        history.push_back(p);
        for(int j = i+1; j < strs.size(); ++j)
        {
            tmp[j].clip(p);
        }
    }

    int i = strs.size()-1;
    if(tmp[i].m.size() > 2) 
    {
        cout << i << " ERROR!!!!!!!!!!!!!!!!!!!!!!!!\n\n\n";
        correct_tree = 0;
        for(auto it = history.begin(); it != history.end(); ++it)
        {
            delete *it;
        }

        return NULL;
    }
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



void getConsensusTree(int &  arg_counter, char* argv[])
{
    string str;
    string arg;
    ifstream gt_file(argv[arg_counter]);

    getline(gt_file, str, ';');
    int indicator = 0;
    int lbl = countNumberTaxa(indicator, str);
   // cout << "indicator: " << indicator << endl;
    int N = lbl;
    int dummy = 0;
    gt_file.close();

    gt_file.open(argv[arg_counter]);
    arg_counter++;

    std::unordered_map <string, int> clpairs;
    if(indicator == 0) //no branch lengths
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
        //cout << vec[i].first << ": " << vec[i].second << endl;
        //cout << vec[i].first.size() << endl;
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

   /* 
    //vector of sorted pairs by their occurence
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
                ss << c - 'A';//if A,B,C,..not t. IMPORTANT: if numtaxa < 10 then both t1,t2,.. and A,B,.. for taxa names ok; if numtaxa >=10 then names that contain more than one char like t1,t2 are not accepted and it is slower
                //ss << c - 'A';//if A,B,C,..not t
            }   

        }
    }


    bool b[vec.size()][vec.size()];
    for(int i = 0; i < vec.size(); ++i)
    {
        for(int j = 0; j < vec.size(); ++j)
        { 
            b[i][j] = cmpt(a[i][0],a[j][0]);
        }
    }

    vector <vector<string>> res_vec;
    vector <string> init_res_vec;
    vector <int> m_index;
    m_index.push_back(0);
    m_index.push_back(1);
    int idx = 0; 
    int max_freq = a[0][1] + a[1][1];
    init_res_vec.push_back(vec[0].first);
    init_res_vec.push_back(vec[1].first);
    for(int i = 2; i < vec.size(); ++i)
    {
        idx = 0;
        for(int j = 1; j < m_index.size(); ++j)
        { 
            if(b[m_index[j]][i] == 0) idx++;
        }
        if(idx == 0) 
        {
            m_index.push_back(i);
            max_freq += a[i][1];

            init_res_vec.push_back(vec[i].first);
        }
    }

    //cout << "MAX_FREQ: " << max_freq << endl;
    res_vec.push_back(init_res_vec);

    int deep_max = 0;
    int counter = 0;
    string str_top="";
    string old_top="";
    ofstream topo_file("outGreedyCons.txt");
    ofstream freq_file("outGreedyConsFreqs.txt");

    string s_temp = " ";
    vector <string> vunique;
    vector <int> res_vec_index;
    for(int j = 0; j < res_vec[0].size(); ++j)
    {
        s_temp += res_vec[0][j];
        s_temp += "-";
    }
    vunique.push_back(s_temp);
    res_vec_index.push_back(0);

    for (int global_i = 0; global_i < res_vec_index.size(); ++global_i) 
    {
        

        int msize = res_vec[res_vec_index[global_i]].size();
        string m[msize];

        for(int j = 0; j < msize; ++j)
        {   
            m[j] = res_vec[res_vec_index[global_i]][j];
        }



        string temp = "";
        sortStringsBySize(m, msize);
        unordered_map <int, int> mfreqs;
        for (int i = 0; i < msize; ++i) 
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


        vector<vector<string>> strs;
        for (int i = 0; i < msize; ++i) 
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
                    strs[i].push_back(strtemp);
                    ++j;

                }
                else
                {
                    ++j;
                }

            }

        }
 
        int correct_tree = 1;
        Node * tree = makeTree(strs, correct_tree);
        string finstr = "";
        if(correct_tree == 1)
        {
            writeUnrankTreeIntoStr(tree, N, finstr);
            deleteTree(tree);

            topo_file << finstr << ";" << endl;
            freq_file << max_freq << " " << finstr << ";" << endl;
          //  cout << max_freq << " " << finstr << ";" << endl;
            }
} 
topo_file.close();
freq_file.close();

}

BagOfNodes::BagOfNodes(vector<string> str)
{
    for(unsigned i = 0; i < str.size(); ++i)
    {
        m.push_back(new Node(str[i]));
    }
}
BagOfNodes::BagOfNodes(vector<string> str, const BagOfNodes & b)
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
}
void BagOfNodes::clip(Node * rule)
{
    auto it = m.end();
    it = find(m.begin(), m.end(), rule->left);
    if (it != m.end())
    { 
        m.erase(it);
    }
    else
    {
        return;
    }
    it = find(m.begin(), m.end(), rule->right);
    if (it != m.end())
    {
        m.erase(it);
        m.push_back(rule);
    }
}
