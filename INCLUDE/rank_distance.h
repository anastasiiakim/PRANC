#ifndef RD_H
#define RD_H
#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <map>
#include <cmath>
#include <algorithm>

void writeTopoToMap(unordered_map <string, int> & clade_and_rank, string str);
void getRankDistance(int &  arg_counter, char* argv[]);

#endif // RD_H
