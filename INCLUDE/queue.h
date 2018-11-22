#ifndef QUE_H
#define QUE_H

#include "node.h"
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
#endif //QUE_H

