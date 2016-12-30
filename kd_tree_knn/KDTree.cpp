#include "KDTree.h"

KDTree::KDTree(long __dim, long __originalDim)
{
    dim = __dim;
    originalDim = __originalDim;
    Node::dim = QueryBox::dim = dim;
    Node::originalDim = originalDim;
    root = NULL;
    height = 0;
}

void KDTree::insert(float *point, float * originalPoint, long dataPtr)
{
    Node * newNode = new Node(point, originalPoint, dataPtr);
    if (root == NULL)
    {
        root = newNode;
        height = 1;
    }
    else
    {
        long insertHeight = root->insert(newNode, 1);
        if (insertHeight > height)
        {
            height = insertHeight;
        }
    }
}

ResultSet KDTree::boxQuery(QueryBox qbox)
{
}

void KDTree::printInOrder()
{
    if (root==NULL)
    {
        cout<<"The tree is empty"<<endl;
    }
    else
    {
        cout<<"Inorder traversal of the tree is"<<endl;
        root->inorder();
    }
}

ResultSet KDTree::rangeQuery(QueryBox *qbox, float *center, float radius)
{

/*    for (long uu=0;uu<originalDim;uu++)
    {
        cout<<center[uu]<<" , ";
    }
    cout<<endl;
    for (long uu=0;uu<dim;uu++)
    {
        cout<<qbox->cntBox[2*uu]<<":"<<qbox->cntBox[2*uu+1]<<"  ";
    }
    cout<<endl;
    exit(1);*/
    ResultSet rs;
    rs.countHits = 0;
    rs.countComparisons = 0;
    //rs.nodeAccesses = 0;
    rs.countDistComputations = 0;
    rs.countVisitedNode = 0;

    if (root!=NULL)
    {
        root->rangeQuery(qbox, rs, 0
        , center, radius
        );
    }
    else
    {
        cerr<<"KD Tree is empty. Query cannot be run"<<endl;
    }
    return rs;
}

ResultSet KDTree::knnQuery(float * center, long K, float *distArray)
{
    ResultSet rs;
    rs.countHits = 0;
    rs.countComparisons = 0;
    //rs.nodeAccesses = 0;
    rs.countDistComputations = 0;
    rs.transformedSpaceDistComp = 0;

    priority_queue < pair <void *, float> ,  vector <  pair < void *, float> > , pairComparator > resultQueue;
    float currentRadius = 9999999.;

    priority_queue < queueItem > nodeQueue;

    if (root!=NULL)
    {
        // COnstruct queue item for root
        queueItem rootItem;
        rootItem.nodePtr = root;
        rootItem.distance = 0.;
        rootItem.desc = 0;
        for(long i=0;i<dim;i++)
        {
            rootItem.box[2*i] = 0.;
            // Ideally this should be sqrt (dim) but any value
            // greater than that is also equally good.
            rootItem.box[2*i+1] = dim;
        }
        nodeQueue.push(rootItem);

        // Now consume each item from the queue
        while (!nodeQueue.empty())
        {
            queueItem next(nodeQueue.top());
            // Remove the first element.
            nodeQueue.pop();
            // Do the actual distance computation
            if (next.distance <= currentRadius )
                next.nodePtr->knnQuery(center, currentRadius, rs, &nodeQueue, next.box, K, &resultQueue, next.desc, distArray);
/*            else
            {
                for (long dd=0;dd<originalDim;dd++)
                {
                    cout<<center[dd]<<",";
                }cout<<endl;
                for (long dd=0;dd<originalDim;dd++)
                {
                    cout<<next.box[2*dd]<<" "<<next.box[2*dd+1]<<",";
                }
                cout<<endl;
                cout<<" Dist "<<next.distance<<"   Radius "<<currentRadius<<endl;
            } */
        }

        // Print results
        while (!resultQueue.empty())
        {
            pair <void *, float>  it = resultQueue.top();
            cout<<((Node *)(it.first))->getDataPointer()<<" : "<<it.second<<"   : "<<((Node *)(it.first))->getDistance(center)<<"  >> "<<getL2Distance(((Node *)(it.first))-> getRecord(), center, originalDim)<<endl;
            resultQueue.pop();
        }
    }
    else
    {
        cerr<<"KD Tree is empty. Query cannot be run"<<endl;
    }
    return rs;
}

