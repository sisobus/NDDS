#ifndef KDTREE_H
#define KDTREE_H

#include "Node.h"
#include <queue>


class KDTree
{
    private:
        Node * root;
        long dim;
        long originalDim;
        long height;
    public:
        /**
          * Constructor that builds an empty
          * kd tree with given dimensions.
          */
        KDTree(long __dim, long __originalDim);

        /**
          * Insert a point in the kdtree
          */
        void insert(float * point, float * originalPoint=NULL, long dataPtr=0);

        /**
         * Box query in the KD Tree
         */
        ResultSet boxQuery(QueryBox qbox);

        /**
         * K-NN query in the KD Tree
         * @param center Center of the query
         * @param K number of neighbors
         * @return ResultSet object containing query results
         */
        ResultSet knnQuery(float* center, long K, float *distArray);

        /**
         * Range query in the KD Tree
         * @param center Center of the query
         * @param radius Radius of the query
         * @return ResultSet object containing query results
         */
        ResultSet rangeQuery(QueryBox *qbox, float *center, float radius);


        /**
          * Print the in-order traversal of the kd tree.
          */
        void printInOrder();

        /**
         * Returns height of the tree
         */
        long getHeight() {return height;};
        
};

#endif
