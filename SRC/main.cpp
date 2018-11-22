#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector> 
#include <iomanip> 

#include "node.h"
#include "string_manipulations.h"
#include "queue.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"

using namespace std;

int Node::st_count = 0;

int main(int argc, char * argv[])
{   
    int arg_counter = 1;
    while(arg_counter < argc)
    {
        if(strcmp(argv[arg_counter],"-rprob") == 0)
        {
            ++arg_counter;
            calcProbsRankedGtInput(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-uprob") == 0)
        {
            ++arg_counter;
            calcProbsUnrankedGtInput(arg_counter, argv);
        }
        else return 0;
    }
    return 0;
}    
