#ifndef SM_H
#define SM_H
#include <string>
#include <vector>
#include <unordered_map>
using namespace std;

void removeSpaces(string & str);
double getBranchLengths (int  i, string & str);
void countParentheses (std::string & str);
void sortString(int n, string * str);
void sortStringsBySize(string *s, int N);
void countUniqueStrings(istream & infile, ofstream & outfile);
void sortTwoColumns(int n, vector<string> & str, vector<int> & val);
bool less(pair <string, int> p1, pair <string, int> p2);
#endif //SM_H
