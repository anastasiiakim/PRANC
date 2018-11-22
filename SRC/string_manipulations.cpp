#include <string>
#include <iostream>
#include <cmath>
#include "string_manipulations.h"

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
    //int scientific = 0;
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
  


