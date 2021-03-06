#include <string.h>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <vector>
#include <iomanip>
#define eps 1e-8

#include "node.h"
#include "string_manipulations.h"
#include "queue.h"
#include "probs_calculation.h"
#include "ranking_unr_gt.h"
#include "min_ancient_coal.h"
#include "write_ranked_tree.h"
#include "maxlike.h"
#include "maxlike_brent.h"
#include "symbolic.h"
#include "get_unranked_topology.h"
#include "get_ranked_topology.h"
#include "get_greedy_consensus.h"
#include "rank_distance.h"
#include "branch_distance.h"

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
        else if(strcmp(argv[arg_counter],"-sym") == 0)
        {
            ++arg_counter;
            symbolicProbsRankedGtInput(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-like_nonni") == 0)
        {
            ++arg_counter;
            calcLikeNoNNI(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-like_nni") == 0)
        {
            ++arg_counter;
            calcLikeWithNNI(arg_counter, argv);
        }     
        else if(strcmp(argv[arg_counter],"-like_nni_brent") == 0)
        {
            ++arg_counter;
            brent_calcLikeWithNNI(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-mac") == 0)
        {
            ++arg_counter;
            searchCandidateSpTreeTopology(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-write") == 0)
        {
            ++arg_counter;
            writeRankedTree(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-rank_trees") == 0)
        {
            ++arg_counter;
            rankUnrTrees(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-utopo") == 0)
        {
            ++arg_counter;
            outputUnrankedTopology(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-rtopo") == 0)
        {
            ++arg_counter;
            outputRankedTopology(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-cons") == 0)
        {
            ++arg_counter;
            getConsensusTree(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-rank_dist") == 0)
        {
            ++arg_counter;
            getRankDistance(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-coal_ints") == 0)
        {
            ++arg_counter;
            outputCoalIntervals(arg_counter, argv);
        }
        else if(strcmp(argv[arg_counter],"-br_len") == 0)
        {
            ++arg_counter;
            getMSEBranchDistance(arg_counter, argv);
        }
        else return 0;
    }
    return 0;
}
