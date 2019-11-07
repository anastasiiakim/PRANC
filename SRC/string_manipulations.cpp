#include <string>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include "string_manipulations.h"
#include <unordered_map>

using namespace std;

void removeSpaces(string & str)
{
    for (int i = 0; i < (int) str.length(); i++)
    {
        if (str[i] == ' ')
        {
            str.erase((size_t) i, 1);
            i--;
        }
    }
}

double getBranchLengths (int  i, string & str)
{
    double temp = 0.;
    double tmp = 0.1;
    double tempsc = 0.;
    double decim = 0.;
    while(static_cast<int> (str[i] - '0') >= 0 && static_cast<int> (str[i] - '0') <= 9 && i < (int) (str.size()-1))
    {
        temp = temp * 10 + (double) (str[i] - '0');
        if(str[i+1] == '.')
        {
            ++i;
            while(str[i+1] != ',' && str[i+1] != ')' && str[i+1] != ' '  && str[i+1] != ';')
            {
                //added section to treat scient. notation
                if(str[i+1] == 'e' && str[i+2] == '-')
                {
                    i = i+2;
                    while(str[i+1] != ',' && str[i+1] != ')' && str[i+1] != ' '  && str[i+1] != ';')
                    {
                        tempsc = tempsc * 10 + double (str[i+1] - '0');
                        ++i;
                    }
                    return (temp + decim) * pow(10, -tempsc);
                }

                if(str[i+1] == 'e' && str[i+2] == '+')
                {
                    i = i+2;
                    while(str[i+1] != ',' && str[i+1] != ')' && str[i+1] != ' '  && str[i+1] != ';')
                    {
                        tempsc = tempsc * 10 + double (str[i+1] - '0');
                        ++i;
                    }
                    return (temp + decim) * pow(10, tempsc);
                }
                decim +=  (double) (str[i+1] - '0') * tmp;
                tmp *= 0.1;
                ++i;
                if (i >= (int) (str.size()-1)) return temp+decim;
            }
        }
        ++i;
        if (i >= (int) (str.size()-1)) return temp+decim;
    }
    return  temp + decim;	
}


void countParentheses (std::string & str)
{
    int lpcount = 0;
    int rpcount = 0;
    int i = 0;

    while(i < (int) str.size())
    {
        if(str[i] == '(') lpcount++;
        else if(str[i] == ')') rpcount++;
        ++i; 
    }	

    if(lpcount != rpcount)
    {
        cout << "Error: The numbers of opening and closing parentheses do not match" << endl;
    }
}


void sortString(int n, string* str)
{
    string temp;
    for(int i = 0; i < n-1; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            if(str[i] > str[j])
            {
                temp = str[i];
                str[i] = str[j];
                str[j] = temp;
            }
        }
    }   
}

void sortTwoColumns(int n, vector<string> & str, vector<int> & val)
{
    string temp;
    int itemp;
    for(int i = 0; i < n-1; ++i)
    {
        for(int j = i+1; j < n; ++j)
        {
            if(val[i] < val[j])
            {
                itemp = val[i];
                temp = str[i];
                val[i] = val[j];
                str[i] = str[j];
                val[j] = itemp;
                str[j] = temp;

            }
        }
    }   
}
  

void sortStringsBySize(string *s, int N) 
{ 
    for (int i = 1; i < N; i++) 
    { 
        string temp = s[i]; 
        int j = i - 1; 
        while (j >= 0 && temp.size() < s[j].size()) 
        { 
            s[j+1] = s[j]; 
            j--; 
        } 
        s[j+1] = temp; 
    } 
} 
/*
void countUniqueStrings(istream & infile, ofstream & outfile)
{
    vector <string> vstr;
    string str = "";
    while(getline(infile, str, '\n'))
    {
        if(str.size() < 3) break;
        vstr.push_back(str);
    }
 
    vector <int> vcounts;
    for(int i = 0; i < vstr.size(); ++i)
    {
        int count = 1;
        for(int j = i + 1; j < vstr.size(); ++j)
        {
            if(vstr[i].compare(vstr[j]) == 0)
            {
                vstr.erase(vstr.begin() + j);
                ++count;
                --j;
            }
        }
        vcounts.push_back(count);
    }

    sortTwoColumns(vcounts.size(), vstr, vcounts);
    for(int i = 0; i < vcounts.size(); ++i)
    {
        outfile << vstr[i] << '\t' << vcounts[i] << endl;

    }
}



*/

bool compareLess(pair <string, int> p1, pair <string, int> p2)
{
    return p1.second > p2.second;
}

void countUniqueStrings(istream & infile, ofstream & outfile)
{
    std::unordered_map <string, int> upairs;
    string str = "";
    while(getline(infile, str, '\n'))
    {
        if(str.size() < 3) break;
        if(upairs.count(str))
        {
            upairs[str]++;
        }  
        else
        {
            pair <string, int> temp (str, 1);
            upairs.insert(temp);
        }
    }
 
    
    vector <pair<string, int>> v(upairs.begin(), upairs.end());
    sort(v.begin(), v.end(), compareLess);

    for(auto i = v.begin(); i != v.end(); ++i)
    {
        //outfile << vunique[i] << '\t' << vcounts[i] << endl;
        outfile << i->second << '\t' << i->first << endl;

    }
}

