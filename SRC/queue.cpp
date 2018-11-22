#include "queue.h"

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
