#ifndef NODE_H
#define NODE_H
#include <string>
using namespace std;

struct Node
{
    string label;
    string desctaxa;
    Node * left;
    Node * right;
    Node * parent;
    Node * original_parent;
    Node ()
    {
        debug_name = st_count;
        st_count++;
        right = NULL;
        left = NULL;
        parent = NULL;
        original_parent = NULL;
        label = "";
        desctaxa = "";
        rank = 0;
        length = 0.;
        time = 0.;
        distFrRoot = 0.;
        leavesnum = 0.;
        outdegree = 0.;
        //  printf("nodeWithName %d created \n", debug_name);
    };
    Node (string _label)
    {
        debug_name = st_count;
        st_count++;
        right = NULL;
        left = NULL;
        parent = NULL;
        original_parent = NULL;
        label = _label;
        desctaxa = "";
        rank = 0;
        length = 0.;
        time = 0.;
        distFrRoot = 0.;
        leavesnum = 0.;
        outdegree = 0.;
        int instance_count = 0;

        //  printf("nodeWithName %d created \n", debug_name);
    };
    Node (Node *pL, Node *pR)
    {
        debug_name = st_count;
        st_count++;
        parent = NULL;
        original_parent = NULL;
        label = "";
        desctaxa = "";
        rank = 0;
        length = 0.;
        time = 0.;
        distFrRoot = 0.;
        leavesnum = 0.;
        outdegree = 0.;
        left = pL;
        right = pR;
        pL->parent = this;
        pR->parent = this;
        int instance_count = 0; 
        //  printf("nodeWithName %d created \n", debug_name);
    };

    static int st_count;
    int debug_name;
    double length;
    double time;
    int rank;
    double distFrRoot;
    bool isRight;
    int population; // from L to R
    int outdegree;
    int leavesnum;
    int leaf_id;
    int instance_count;
    
};


#endif //NODE_H
