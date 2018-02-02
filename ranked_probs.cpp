#include <string>
#include <fstream>
#include <iostream>
#include <stack>
#include <cmath>
#include <algorithm>
#include <iomanip> // output more decimals
#define eps 0.0e-15
//#define MaxQ 105
using namespace std;

struct Node
{
    string label;
    string desctaxa;
    Node * left;
    Node * right;
    Node * parent;
    Node ()
    {
        debug_name = st_count;
        st_count++;
        right = NULL;
        left = NULL;
        parent = NULL;
        label = "";
        desctaxa = "";
        rank = 0;
        length = 0.;
        time = 0.;
        distFrRoot = 0.;
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
};

int Node::st_count = 0;

inline void removeSpaces(string & str)
{
    for (int i = 0; i < str.length(); i++)
    {
        if (str[i] == ' ')
        {
            str.erase(i, 1);
            i--;
        }
    }
}

inline double getBranchLengths (int  i, string & str)
{
    double temp = 0.;
    double tmp = 0.1;
    double tempsc = 0.;
    double decim = 0.;
    int scientific = 0;
    while((double) (str[i] - '0') >= 0 && (double) (str[i] - '0') <= 9)
    {
        temp = temp * 10 + (double) (str[i] - '0');
        if(str[i+1] == '.')
        {
            ++i;
            while(str[i+1] != ',' && str[i+1] != ')' && str[i+1] != ' '  && str[i+1] != ';')
            {
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
            }
        }
        ++i;
    }

    return  temp + decim;	

}


inline void pushToArray (Node * p, int & tail, Node ** &ar)
{
    int i = tail;
    ar[i] = p;  
    while((i > 0) && (ar[i] -> time < ar[i - 1] -> time))
    {
        Node * tempe;
        tempe = ar[i];
        ar[i] = ar[i-1];
        ar[i-1] = tempe;
        i--;
    }

    tail ++;
}

inline void pushToArrayNoSort (Node * p, int & tail, Node ** &ar)
{
    int i = tail;
    ar[i] = p;   
    tail ++;
}



inline Node * popFromArray (int & tail, Node ** &ar)
{
    tail --;
    return ar[tail];
}



inline void distFromRoot (Node * p)
{	
    if (p -> label.size() == 0)
    {	
        p -> right -> distFrRoot = p -> distFrRoot + p -> right -> length;
        p -> left -> distFrRoot = p -> distFrRoot + p -> left -> length;
        p -> right -> parent = p;
        p -> left -> parent = p;

        distFromRoot(p -> right);
        distFromRoot(p -> left);
    }
    else // for leaves
    {
        p -> outdegree = 0;
        p -> population = 1;
    }

}


void isUltrametric (Node * p, double & tempdist)
{
    if (p == NULL) cout << "Tree is empty" << endl;
    if (p -> label.size() == 0)
    {
        isUltrametric(p -> right, tempdist);
        isUltrametric(p -> left, tempdist);
    }

    else 
    {	
        tempdist = p -> distFrRoot;
    }
}


void printIsUltrametric (Node * p, double tempdist, int & check)
{
    if (p -> label.size() == 0)
    {
        printIsUltrametric(p -> right, tempdist, check);
        printIsUltrametric(p -> left, tempdist, check);
    }

    else 
    {	
        if (fabs(tempdist - p -> distFrRoot) > eps)	
        {
            check++;
        }
    }
}



inline void countParentheses (std::string & str)
{
    int lpcount = 0;
    int rpcount = 0;
    int i = 0;

    while(i < str.size())
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


inline void deleteTree (Node * p)
{
    if (p != NULL )
    {
        deleteTree(p -> right);
        deleteTree(p -> left);
        delete p;
    }
}

inline void deleteStack (std::stack <Node *> stk)
{
    while(!stk.empty())
        stk.pop();
}





inline void pushNodes (int & lbl, std::stack <Node *> & stk, std::string str)
{
    int i = 0;
    while(i < str.size())
    {
        
            if(((str[i] >= 'a' && str[i] <= 'z') || (str[i] >= 'A' && str[i] <= 'Z')) && (str[i+1] != '-' && str[i+1] != '+'))
            {
                std::string strtemp;
                strtemp += str[i]; 
                while(str[i+1] != ':')
                {
                    strtemp += str[i+1];		   
                    ++i;
                }
                Node * newnode = new Node();
                newnode -> label = strtemp;
                newnode -> time = 0.;
                ++lbl;	
                ++i;
                if((double) (str[i+1] - '0') >= 0 && (double) (str[i+1] - '0') <= 9)
                {
                    newnode -> length = getBranchLengths(i+1, str);
                }
                stk.push(newnode);

            }

            else if ( str[i] == ')' )
            {	
                Node * newnode = new Node();
                newnode -> right = stk.top();
                newnode -> right -> isRight = 1;
                stk.pop();
                newnode -> left = stk.top();
                newnode -> left -> isRight = 0;
                stk.pop();
                ++i;

                if((double) (str[i+1] - '0') >= 0 && (double) (str[i+1] - '0') <= 9)
                {

                    newnode -> length = getBranchLengths(i+1, str);
                }

                newnode -> time = newnode -> right -> time + newnode -> right -> length;
                stk.push(newnode);
            } 
            ++i;
        }
    
}	

inline void getRanks(Node* newnode, int & tail, Node ** &ar)
{
    int k = 1, k1 = 1; 
    while(newnode -> time != 0)
    {
        Node * tempnode;
        tempnode = popFromArray(tail, ar);//ar[tail]
        if(tempnode -> rank == 0) tempnode -> rank = k; //added if cond to avoid segfault
        k++;
        k1++;
        if((tail > 0) && (tempnode -> time == ar[tail - 1] -> time) && (tempnode -> time > 0))
        {
            int x;

            //    cout << "Warning: some ranks are tied. Rand assignment was applied" << endl;		
            k--;	
            x = rand()%(k1-k+1)+k;
            tempnode -> rank = x;
            if(x == k) ar[tail-1] -> rank = k1; 
            else ar[tail-1] -> rank = k;

        }

        else if((tail > 0) && (tempnode -> time != ar[tail - 1] -> time))  
        {
            k = k1;
        }

        pushToArray(tempnode -> right, tail, ar);
        pushToArray(tempnode -> left, tail, ar);

        newnode = ar[tail - 1];
    }

}	


inline Node * getNodeFromRank (Node * p, int rankValueGT)
{
    if(p == NULL) cout << "Tree is empty" << endl;
    Node * k;
    Node * t;
    if(p -> label.size() == 0)
    {		
        if(p -> rank != rankValueGT) 
        {
            t =  getNodeFromRank(p -> right, rankValueGT);
            k =  getNodeFromRank(p -> left, rankValueGT);
            if(t != NULL) return t;
            return k;
        }
        else 
        {
            return p;
        }	
    }
    else return 0;

}	




inline void getTaxa (Node * p, std::stack <Node *> & allTaxa)
{
    if(p -> label.size() == 0)
    { 
        getTaxa(p -> right, allTaxa);
        getTaxa(p -> left, allTaxa);
    }
    else 
    {
        allTaxa.push(p);
    }	
}


inline Node * mrcaST (Node * p, Node * rch, Node * lch)
{
    if(p == NULL) return NULL;
    if(p -> label == rch -> label || p -> label == lch -> label)
    {
        return p;
    }
    else
    {
        Node * r = mrcaST(p -> right, rch, lch);
        Node * l = mrcaST(p -> left, rch, lch);

        if(l!=NULL && r!=NULL)//&&l->rank==0&&r->rank==0)
        {

            return p;
        }
        else
        {
            if(l!=NULL) 
            { 
                return l;
            }
            return r;
        }
    }	
}

inline Node * mrcaST2 (Node * p, Node * tp, std::string v)
{
    if(p == NULL) return NULL;
    if( p == tp || p -> label == v) return p;
    else 
    {
        Node * l = mrcaST2(p -> left, tp, v);	
        Node * r = mrcaST2(p -> right, tp, v);

        if(l&&r) return p;
        else
        {
            if(l) return l;
            return r;
        }
    }

}



Node * getMrca(Node * p, std::stack <Node * > & stk)
{
    Node * a = stk.top();
    stk.pop();
    Node * b = stk.top();
    stk.pop();
    Node * ptemp = mrcaST(p, a, b);
    while(!stk.empty())
    {
        if(ptemp -> rank == 1) return ptemp;
        Node * t = stk.top();
        stk.pop();
        ptemp = mrcaST2(p, ptemp, t -> label);	
    }	
    return ptemp;
}	

inline void storeLabels (Node * p, int N, string * ar, int & i)
{
    if(p -> label.size() != 0)
    {
        storeLabels(p -> right, N, ar, i);
        storeLabels(p -> left, N, ar, i);
    }
    else
    {
        ar[i] = p -> label;
        ++i;
    }
}


class Queue
{
    private:
        Node ** data;
        size_t maxsize;
        int head;
        int tail;
    public:
        Queue(size_t ms)
        {
            head = 0; 
            tail = 0;
            maxsize = ms;
            data = new Node * [maxsize];
        }
        ~Queue()
        {
            delete[] data;
        }
        Node * pop();
        void push(Node * p);
        void arrange();
        Node * &operator[](int i);
        int size ();
        bool isEmpty();
};

bool Queue::isEmpty()
{
    if (head == tail) {return true;}
    else {return false;}
}	

int Queue::size ()
{
    if (head > tail) return (maxsize - head + tail);
    return (tail - head);
}

Node * Queue::pop ()
{
    if(!isEmpty())
    {
        head++;
        return data[head - 1]; 
    }
    else return 0;
}

Node * &Queue::operator[] (int i)
{
    if (head > tail) return data[i - maxsize + head];
    return data[i + head];
}

void Queue::arrange()
{
    int i = tail - 1;
    while (data[i+1] -> distFrRoot < data[i] -> distFrRoot)
    {
        Node * temp = data[i];
        data[i] = data[i+1];
        data[i+1] = temp;
        i--;
    }
}


void Queue::push (Node * p)
{
    data[tail] = p;
    if(tail != head)
    {
        arrange();
    }
    tail++;
}


inline void addClone(Node * popnode, Queue & q)
{
    for (int i = 0; i < q.size(); ++i)
    {

        Node * inode = q[i];
        Node * parnd = inode -> parent;
        Node * snode = new Node();
        snode -> length = popnode -> length;
        inode -> length -= popnode -> length;
        inode -> parent = snode;
        snode -> parent = parnd;
        snode -> rank = popnode -> rank;
        snode -> distFrRoot = popnode -> distFrRoot;
        snode -> time = popnode -> time;
        if (inode -> isRight == 0) 
        {
            snode -> left = inode;
            parnd -> left = snode;
        }
        else
        {
            snode -> right = inode;
            parnd -> right = snode;

        }

    }

}

inline void makeBeadedTree (Node * root, int MaxQ)
{
    Queue q(MaxQ);
    q.push(root);
    while(!q.isEmpty())
    {	
        Node * nd = q.pop();
        if(nd -> rank != 1 && nd -> rank != 0) addClone(nd, q);
        if(nd -> left != NULL && nd -> right != NULL)
        {
            q.push(nd -> left);
            q.push(nd -> right);
        }	
    }

}       	


inline void makeBeadedArray (Node * root, int ** y, std::stack <Node *> & stk, int MaxQ)
{
    Queue q(MaxQ);
    q.push (root);
    int j = 1;
    int t = 1;
    root -> population = 1;
    root -> outdegree = 2;
    stk.push(root);
    y[root -> rank][1] = 2;
    while (!q.isEmpty())
    {
        Node * nd = q.pop();
        if (j >= nd -> rank && t != nd -> rank) j = 1;
        if (nd -> rank != 0 && nd -> rank != 1) 
        {
            if (nd -> left != NULL && nd -> right != NULL)	
            {
                y[nd -> rank][j] = 2;
                nd -> outdegree = 2;
                nd -> population = j;
                stk.push(nd);
            }
            else if (nd -> left != NULL || nd -> right != NULL) 
            {
                y[nd -> rank][j] = 1;
                nd -> outdegree = 1;
                nd -> population = j;
            }
            j++;
            t = nd -> rank;
        }

        if (nd -> left != NULL) q.push(nd -> left);
        if (nd -> right != NULL) q.push(nd -> right); 
    }	
}


inline void getPath (Node * p, std::stack <Node *> & stk)
{
    stk.push(p);	
    if(p -> rank != 1)
    {	
        p = p -> parent;
        getPath(p, stk);
    }
}

inline void getPathMatchRank(Node * &  p, int rank_val) 
{
    if (p -> rank != rank_val && p -> rank != 1) 
    {	
        p = p -> parent;
        getPathMatchRank(p, rank_val);
    }

}	


inline int * getNextHistory (int * ar, int * ar_orig, int n, int idx, int & numb)
{
    if (*(ar + n - 2) == 1) return ar;
    for(int t = n - 2; t >=0; t--)	
    {
        if(*(ar + t - 1) < *(ar + t))
        {
            *(ar + t) -= 1;
            for(int j = t + 1; j <= n - 2; ++j)
            {	
                *(ar + j) = *(ar_orig + j);
            }
            numb++;
            if (numb == idx) return ar;
            else
            {
                ar = getNextHistory(ar, ar_orig, n, idx, numb);
            }
        }
    }
    return ar;
}



inline void arrayFreq(int * arr, int * res, int n)
{
    for (int j =0; j<n; j++)
        arr[j] = arr[j]-1;
    for (int i=0; i<n; i++)
        arr[arr[i]%n] = arr[arr[i]%n] + n;
    for (int i =0; i<n; i++)
        res[i] = arr[i]/n;	   	 
}	


inline void get_hseq(int * arseq, Node * nodeGT, Node * nodeST, int n)
{
    arseq[0] = 1;
    for(int i = 2; i < n; ++i)
    {	
        std::stack <Node *> allTaxa;
        getTaxa(getNodeFromRank(nodeGT, i), allTaxa);
        Node * mrcanode  = getMrca(nodeST, allTaxa);
        if (mrcanode -> rank < i) 
        {    
            arseq[i-1] = mrcanode -> rank;
        }
        else arseq[i-1] = i;
    }
}

inline void reverseArraySort(int * arseq, int n)
{
    for(int i = n - 2; i > 0; --i)
    {
        if(arseq[i-1] > arseq[i])
            arseq[i-1] = arseq[i]; 
    }
}


inline void getDescTaxa(Node * node, int N)
{
    if(node -> label.size() == 0)
    {
        stack <Node *> allTaxa;
        getTaxa(node, allTaxa);
        string temp;
        while(!allTaxa.empty())
        {
            temp += allTaxa.top() -> label;
            allTaxa.pop();
        }
        node -> desctaxa = temp;
        getDescTaxa(node -> right, N);
        getDescTaxa(node -> left, N);

    }
}

inline void returnMrca(Node * nodeST, Node * nodeGT, int n, Node ** ar)
{
    for(int i = 1; i < n; ++i)
    {
        std::stack <Node *> allTaxa;
        getTaxa(getNodeFromRank(nodeGT, i), allTaxa);
        Node * mrcanode  = getMrca(nodeST, allTaxa);
        ar[i-1] = mrcanode;
    }   
}

inline void get_node(int * rank_hist, int ** temp, int n, Node ** ar)
{
    int tmp = 0;
    int j = 0;
    Node * mrcanode;
    for(int i = 1; i < n; ++i)
    {
        mrcanode = ar[i-1];
        getPathMatchRank(mrcanode, rank_hist[i-1]);
        if (tmp != mrcanode -> rank) j = 1;
        temp[mrcanode->rank-1][j-1] = mrcanode -> population;
        tmp = mrcanode -> rank;
        ++j;
    }
}

inline void calc_k_ijz (int ***k, int n, int * m, int ** y, int ** ar_r)
{
    for(int i = 1; i <= n; ++i) 
    {
        k[n][0][i] = 1;

    } 	

    for(int i = n - 1; i >= 1; --i)
    {
        int tempz = 1;
        for(int z = 1; z <= i; ++ z)
        {
            for(int t = 1; t <= y[i][z]; ++t) //outdegree 
            {
                k[i][m[i-1]][z] += k[i+1][0][tempz]; // (9)
                ++tempz;
            }
            for(int j = m[i-1] - 1; j >= 0; j--)
            {
                for(int h = j; h >=0; h--)
                {
                    if (ar_r[i-1][h] == z)
                        k[i][h][z] = k[i][h+1][z] - 1;
                    else
                        k[i][h][z] = k[i][h+1][z];
                }

            }
        }
    }
}       	


inline void writeTree(Node* root)
{
    cout << "Name " << root->debug_name;
    if (root->left) cout << " left child: " << root->left->debug_name;
    if (root->right) cout << " right child: " << root->right->debug_name;
    if (root->parent) cout << " parent " << root->parent->debug_name;
    cout << " rank " << root->rank;
    cout << " label " << root -> label;
    cout << " desctaxa " << root -> desctaxa;
    cout << endl;
    //	cout << " time " << root -> time;
    //	cout << " population " << root -> population << endl;
    //	cout << " outdegree " << root -> outdegree << endl;


    if (root->left) writeTree(root->left);
    if (root->right) writeTree(root->right);
}

inline int factorial (int num)
{
    long long accum = 1;
    while (num > 0)
    {
        accum *= num;
        num--;
    }
    return accum;
}


inline int calcBinomCoef (int n, int k)
{
    double prod = 1.;
    for(int i = 1; i <= k; ++i)
        prod *= (double) (n + 1 - i)/i;
    return prod;
}



inline void sortString(int n, string* str)
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


inline void storeLabels(int n, int & i, string* str, Node * p)
{
    if(p -> label.size() == 0)
    {
        storeLabels(n, i, str, p -> right);
        storeLabels(n, i, str, p -> left);
    }
    else
    {
        str[i] = p -> label;
        ++i;
    }
}

inline void getS(Node * p, double * s)
{
    if(p -> label.size() == 0)
    {
        s[p -> rank - 1] = p -> time;
        getS(p -> right, s);
        getS(p -> left, s);
    }
}

inline int lambda_ij (int i, int j, int *** k)
{
    int res = 0.;
    for(int z = 1; z < i + 1; ++z)
        res += calcBinomCoef(k[i][j][z], 2);
    return res;
}

inline double invPartialCoal(int j)
{
    return pow(2, j-1) / factorial(j) / factorial(j-1);
}

inline double conditionalProbability(int i, int * m_i, int *** array_k, double * s) // min i = 2 since l[i-2] exists
{
    double res = 0.;
    int prod;
    double m = m_i[i-1];
    double no_coal_prod = 1.;
    for(int j = 0; j < m + 1; ++j)
    {
        prod = 1;
        for(int k = 0; k < m + 1; ++k)
        {
            if(k != j)
            {
                prod *= lambda_ij(i, k, array_k) - lambda_ij(i, j, array_k);
            }
        }
        res += exp(-lambda_ij(i, j, array_k)*(s[i-2] - s[i-1]))/prod;
    }

    return res;
}

inline int max(int & x, int & y)
{
    if(x > y) return x;
    else return y;
}



inline void getInvPartialCoal(double * ar, int N)
{
    for(int i = 2; i < N + 1; ++ i)
    {
        ar[i-2] = pow(2, i - 1) / factorial(i) / factorial(i - 1);
    }
}


inline double geneTreeProbability(int * m, int *** k, double * s, double * coal, int n)
{
    double res = 1.;
    for(int i = 2; i < n; ++i)
    {
        res *= conditionalProbability(i, m, k, s);
    }
    return res*coal[m[0] - 1]; 

}

int depth(Node *node)
{
    int d = -1;
    while (node)
    {
        ++d;
        node = node->parent;
    }
    return d;
}

// To find LCA of nodes n1 and n2 in Binary Tree
Node *LCA(Node *n1, Node *n2)
{
    int d1 = depth(n1), d2 = depth(n2);
    int diff = d1 - d2;

    if (diff < 0)
    {
        Node * temp = n1;
        n1 = n2;
        n2 = temp;
        diff = -diff;
    }

    while (diff--)
        n1 = n1->parent;

    while (n1 && n2)
    {
        if (n1 == n2)
            return n1;
        n1 = n1->parent;
        n2 = n2->parent;
    }

    return NULL;
}



int main()
{ 

    //ST
    std::string str;
    ifstream fin("ST.txt");
    std::getline(fin, str, '\n');	
    fin.close();   
    removeSpaces(str);
 //   cout << str << endl;
    std::stack <Node *> stkST;

    int lbl = 0;
    double tempdist = 0.;	
    int labelcount = 0;
    countParentheses(str); 

    pushNodes(lbl, stkST , str);
    Node * newnode;
    newnode = stkST.top();
    newnode -> distFrRoot = 0.;
    int tail = 0;

    distFromRoot(newnode);
    isUltrametric(newnode, tempdist);
    /*
       int check = 0;
       printIsUltrametric(newnode, tempdist, check);
       if(check != 0)
       {
       cout << "ERROR: The species tree is NOT ultrametric" << endl;
       return 0;
       }
       */	
    int N = lbl;
    Node ** ar = new Node * [N];
    int ** ar_y = new int * [N];
    for (int i = 0; i < N; i++) ar_y[i] = new int [N];

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
        {
            ar_y[i][j] = 0;
        }

    pushToArray(newnode, tail, ar);	
    getRanks(newnode, tail, ar);


   // writeTree(newnode);

    int lblscounter = 0;
    string strLabels[N];

    double s[N-1];
    getS(newnode, s);
    //   cout << "time s_i: " << endl;
    //  for(int i = 0; i < N-1; ++i) cout << s[i] << endl;

    int maxQue = (double) (1+N)/2 * N;
    std::stack <Node *> origRanked;
    makeBeadedTree(newnode, maxQue);
    makeBeadedArray(newnode, ar_y, origRanked, maxQue);

    //GT

    string strGT;
    string strTops;
    ifstream finGT("gtuniqtrees.txt");
    ifstream fintops("gtuniqtops.txt");
    ofstream finprobGT("probForEachGT.txt");

    double total_prob = 0;
    int treecounter = 0;
  
    Node ** arMrca = new Node * [N-1];
    int ** ar_rankns = new int * [N-1];
    for (int i = 0; i < N-1; i++) ar_rankns[i] = new int [N-1];
    for(int i = 0; i < N-1; ++i)
    {
        for(int j = 0; j < N - 1; ++j)
        {
            ar_rankns[i][j] = 0; 
        }
    }	

    int ***k;
    k = new int **[N+1];
    for(int i = 0; i < N+1; i++)
    {
        k[i] = new int *[N];
        for(int j = 0; j < N; ++j)
            k[i][j] = new int[N+1];
    }


    double array_invcoal[N-1];
    getInvPartialCoal(array_invcoal, N);
    while(getline(finGT, strGT, ';'))
    {
        treecounter++;
        getline(fintops, strTops, '\n');
        if(strGT.size() < 3) break;
        removeSpaces(strGT);
        std::stack <Node *> stkGT;

        int lblGT = 0;
        double tempdistGT = 0.;
        int labelcountGT = 0;
        //  countParentheses(strGT);

        Node **arGT = new Node * [N];
        pushNodes(lblGT, stkGT , strGT);
        Node * newnodeGT;
        newnodeGT = stkGT.top();
        newnodeGT -> distFrRoot = 0.;
        int tailGT = 0;
        distFromRoot(newnodeGT);

        pushToArray(newnodeGT, tailGT, arGT);	
        getRanks(newnodeGT, tailGT, arGT);

        getDescTaxa(newnodeGT, N);
        //        writeTree(newnodeGT);
        returnMrca(newnode, newnodeGT, N, arMrca);



        int * arRankHistory = new int [N]; 
        get_hseq(arRankHistory, newnodeGT, newnode, N);
        reverseArraySort(arRankHistory, N);
        /*   cout << "Max Ranked History" << endl;
             for(int i = 0; i < N-1; ++i) cout << arRankHistory[i] << " ";
             cout << endl;	
             */
        int arRankHistoryCopy[N];
        for(int i = 0; i < N; ++i)	
            arRankHistoryCopy[i] = arRankHistory[i];

        int * m_i = new int[N-1];
        int max_ar_irows = *max_element(arRankHistory, arRankHistory+N-1); 

        arrayFreq(arRankHistoryCopy, m_i, N-1);//1st arg changes after running this func
        int max_m_i = *max_element(m_i, m_i+N-1); 

        get_node(arRankHistory, ar_rankns, N, arMrca);
        for(int i = 0; i < N+1; i++)
            for(int j = 0; j < N; ++j)
                for(int z = 0; z < N+1 ; ++z)
                {	
                    k[i][j][z] = 0;
                } 
        calc_k_ijz(k, N, m_i, ar_y, ar_rankns);//*, m_i);

        double prob_val = geneTreeProbability(m_i, k, s, array_invcoal, N);
        int arRHcopy[N-1];
        for(int i = 0; i < N - 1; ++i)
            arRHcopy[i] = arRankHistory[i];


        int history_counter = 1;
        int numb = 0;
        int * next_history;
        while(*(arRHcopy + N - 2) != 1)
        {
            next_history = getNextHistory (arRHcopy, arRankHistory, N, history_counter, numb);
            /*
               cout << "---------- Next Label history ---------- " << endl;
               for(int i = 0; i < N - 1; ++i)
               cout << next_history[i] << " ";
               cout << endl;
               */
            //repeat       
            for(int i = 0; i < N; ++i)	
                arRankHistoryCopy[i] = next_history[i];
            max_ar_irows = *max_element(next_history, next_history+N-1); 
            arrayFreq(arRankHistoryCopy, m_i, N-1);
            max_m_i = *max_element(m_i, m_i+N-1);

            get_node(next_history, ar_rankns, N, arMrca);

            for(int i = 0; i < N+1; i++)
            for(int j = 0; j <= max_m_i ; ++j)
                for(int z = 0; z < N+1 ; ++z)
                {	
                    k[i][j][z] = 0;
                }
            calc_k_ijz(k, N, m_i, ar_y, ar_rankns);//*, m_i);
            prob_val += geneTreeProbability(m_i, k, s, array_invcoal, N);
            //     cout << "get_prob = " << geneTreeProbability(m_i, k, s, N)  << " " << log(geneTreeProbability(m_i, k, s, N)) <<  endl;
            //   cout << "prob_of_history = " << setprecision(10) << prob_val << " " << log(prob_val) << endl;
            history_counter++;
        
        }
        total_prob += prob_val;
        finprobGT << setprecision(15)  << prob_val <<  '\t' << strTops << '\n'; 

        //	cout << "Prob of GT = " << setprecision(15) << prob_val << endl;

       //  cout << "Total prob = " << total_prob << endl;
        delete[] arRankHistory;
        delete[] m_i;

        delete [] arGT;

        deleteStack(stkGT);

    }
  
    finGT.close();
    finprobGT.close();

    for(int i = 0; i < N-1; ++i)
        delete[] ar_rankns[i];
    delete [] ar_rankns;

    for(int i = 0; i < N; ++i)
    {
        delete[] ar_y[i];	
    }
    delete[] ar_y;	

    delete[] ar;
    delete[] arMrca;
    for(int i = 1; i < N+1; ++i)
    {
        for(int j = 0; j < N; ++j)
            delete [] k[i][j];
        delete [] k[i];
    }
    delete [] k;

    fintops.close();
    deleteStack(stkST);

    return 0;
    } 

