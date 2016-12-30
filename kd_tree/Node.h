#ifndef NODE_H
#define NODE_H

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <queue>
#include <utility>
#include <sstream>
#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <cassert>
#include "QueryBox.h"

class queueItem;

// helper function for our priority queue
class pairComparator
{
    public:
    bool operator() ( pair <void *, float> p1, pair <void *, float> p2 )
    {
        return p1.second < p2.second;
    }
};

struct ResultSet
{
    long countHits;
    long countComparisons;
    long countDistComputations;
    long transformedSpaceDistComp;
    long countVisitedNode;
  // long nodeAccesses;
    //vector <long> hitList;
};

/** Function that calculates the L-infinity distance between two
 * points (arrays of floats)
 * @param v1 First point
 * @param v2 Second point
 * @dimension Number of dimensions
 * @return L2 distance.
 */
float getLinfDistance(float * v1, float* v2, long dimension);

/** Function that calculates the L2 distance between two
 * points (arrays of floats)
 * @param v1 First point
 * @param v2 Second point
 * @dimension Number of dimensions
 * @return L2 distance.
 */
float getL2Distance(float * v1, float* v2, long dimension);

float getBoxL2Distance(float * v1, float* box, long dimension);

float getBoxLinfDistance(float * v1, float* box, long dimension);

// This class defines a kd-tree node.
class Node
{
    private:
        // Data Record
        float *record;
        // Additional data pointer
        long dataPtr;
        // Left and right children
        Node *left;
        Node *right;
        // For the transformation approach, one may wish to
        // keep the original data vector
        float * originalDataVector;
    public:
        // Number of dimensions
        static long dim;
        // Number of dimensions in the original data
        static long originalDim;

        /**
          * Conctructor that builds the node using 
          * an array of floating point values
          */
        Node(float * __dataArray, long __dataPtr=0);

        float *getRecord() {return originalDataVector;}

        /**
          * Conctructor that builds the node using 
          * an array of floating point values.
          * In addition, it also keeps track of original data vector.
          */
        Node(float* __dataArray, float * __originalData, long __dataPtr=0);

        /**
          * Conctructor that builds the node using 
          * a string containing comma separated floating point values
          */
//        Node(string    __dataString, long __dataPtr=0);

        /**
          * Insert a child node.
          * @param newNode The new node to insert
          * @param height Height of the current node.
          * @param desc Descriminator dimension
          * @return Height at which the node was inserted.
          */
        long insert(Node *newNode, long height, long desc=0);

        /**
          * Destructor
          */
        ~Node();

        /** In order traversal of the kd-tree
          */
        void inorder();

        /**
          * Copy constructor
          * @param N The node to be copied
          */
        Node(const Node&N);

        /**
          * Assignment operator. Performs a deep copy.
          */
        Node & operator=(const Node &second);

        /**
         * Returns L2 distance of this node from other point
         */
        float getDistance(float *vec);

        long getDataPointer()
        {
            return dataPtr;
        }

        /**
         * Range query in KD-Tree node. Note that we do not
         * really do a range query. Instead we are going to do a
         * box query in transformed space.
         * @param qbox Query box
         * @param rs ResultSet object which record query execution stats
         * @param desc Descriminator dimension
         */
        void rangeQuery(QueryBox *qbox, ResultSet &rs, long desc
                , float *center, float radius // May be used to eliminate false positives.
                );

        // Returns distance of "TRUE" dimensions of this node
        // from the center.
        void knnQuery(float *center, float &radius, ResultSet &rs, priority_queue  <  queueItem > *nodeQ, float * boundingBox, const long K, priority_queue < pair<void *, float > , vector < pair < void *, float > > , class pairComparator >  *resultQ, long desc, float *distArray);
};

class queueItem
{
    public:
    Node * nodePtr;
    float * box;
    float distance;
    long desc;
    
    queueItem()
    {
        box = new float[2*Node::dim];
    }
    
    ~queueItem()
    {
        delete[] box;
    }

    queueItem(const queueItem &q2)
    {
        box = new float[2*Node::dim];
        nodePtr = q2.nodePtr;
        distance = q2.distance;
        desc = q2.desc;
        for (long i=0;i<2*Node::dim;i++)
        {
            box[i] = q2.box[i];
        }
    }

    queueItem & operator=(const queueItem q2)
    {
        box = new float[2*Node::dim];
        nodePtr = q2.nodePtr;
        distance = q2.distance;
        desc = q2.desc;
        for (long i=0;i<2*Node::dim;i++)
        {
            box[i] = q2.box[i];
        }
        return *this;
    }

    bool operator<(const queueItem q2) const 
    {
        return distance < q2.distance;
    }
};

#endif
