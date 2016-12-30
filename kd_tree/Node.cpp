#include "Node.h"

long Node::dim = 2;
long Node::originalDim = 2;


// Implementation of the distance function.
float getLinfDistance(float * v1, float* v2, long dimension)
{
    float distance = 0.;
    for (long i=0;i<dimension;i++)
    {
        if (fabs(v1[i] - v2[i]) > distance)
        {
            distance = fabs(v1[i] - v2[i]);
        }
    }
    return distance;
}

// Implementation of the distance function.
float getL2Distance(float * v1, float* v2, long dimension)
{
    float distance = 0.;
    for (long i=0;i<dimension;i++)
    {
        distance += (v1[i] - v2[i])*(v1[i] - v2[i]);
    }
    distance = sqrt(distance);
    return distance;
}

float getBoxLinfDistance(float *v1, float *box, long dimension)
{
    float distance = 0.;
    for (long i=0;i<dimension;i++)
    {
        float d = 0.;
        if (v1[i] < box[2*i])
        {
            d = (box[2*i] - v1[i]);
        }
        else if (v1[i] > box[2*i+1])
        {
            d = (v1[i] - box[2*i+1]);
        }
        if (d > distance)
        {
            distance = d;
        }
    }
    return distance;
}

float getBoxL2Distance(float *v1, float *box, long dimension)
{
    float distance = 0.;
    for (long i=0;i<dimension;i++)
    {
        if (v1[i] < box[2*i])
        {
            distance += (box[2*i] - v1[i]) * (box[2*i] - v1[i]);
        }
        else if (v1[i] > box[2*i+1])
        {
            distance += (v1[i] - box[2*i+1])*(v1[i] - box[2*i+1]);
        }
    }
    distance = sqrt(distance);
    return distance;
}

// Basic constructor
Node::Node(float *__dataArray, long __dataPtr)
{
    record = new float[dim];
    for (long i=0;i<dim;i++)
    {
        record[i] = __dataArray[i];
    }
    dataPtr = __dataPtr;
    left = NULL;
    right = NULL;
    originalDataVector = NULL;
}

// Basic constructor which also tracks the original data
Node::Node(float *__dataArray, float * __originalData, long __dataPtr)
{
    record = new float[dim];
    for (long i=0;i<dim;i++)
    {
        record[i] = __dataArray[i];
    }
    dataPtr = __dataPtr;
    left = NULL;
    right = NULL;
    if (__originalData != NULL)
    {
        originalDataVector = new float[originalDim];
        for (long i=0;i<originalDim;i++)
        {
            originalDataVector[i] = __originalData[i];
        }
    }
    else
    {
        originalDataVector = NULL;
    }
}

// Costructor using a string for input data
/*Node::Node(string __dataString, long __dataPtr)
{
    stringstream sin(__dataString);
    record = new float[dim];
    for (long i=0;i<dim;i++)
    {
        sin>>record[i];
        if (!sin.good())
        {
            cerr<<"Error processing the line:\n"<<__dataString<<endl;
            exit(1);
        }
    }
    dataPtr = __dataPtr;
    left = NULL;
    right = NULL;
    originalDataVector = NULL;
}*/

// Destructor
Node::~Node()
{
    delete[] record;
    if (originalDataVector != NULL)
    {
        delete[] originalDataVector;
    }
}

long Node::insert(Node *newNode, long height, long desc)
{
    if (newNode->record[desc] <= record[desc])
    {
        // Left subtree
        if (left == NULL)
        {
            left = newNode;
            return height+1;
        }
        else
        {
            return left->insert(newNode, height + 1, (desc +1)%dim);
        }
    }
    else
    {
        // Right subtree
        if (right == NULL)
        {
            right = newNode;
            return height+1;
        }
        else
        {
            return right->insert(newNode, height+1, (desc +1)%dim);
        }
    }
}

void Node::inorder()
{
    if (left != NULL)
    {
        left->inorder();
    }
    cout<<dataPtr <<"  : ";
    for (long i=0;i<originalDim;i++)
    {
        cout<<originalDataVector[i]<<" ";
    }
    cout<<" <==> ";
    for (long i=0;i<dim;i++)
    {
        cout<<record[i]<<" ";
    }
    cout<<endl;
    if (right != NULL)
    {
        right->inorder();
    }
}

float Node::getDistance(float * vec)
{
    return getL2Distance(originalDataVector, vec, originalDim);
}

// Copy constructor
Node::Node(const Node&N)
{
    record = new float[dim];
    for (long i=0;i<dim;i++)
    {
        record[i] = N.record[i];
    }
    dataPtr = N.dataPtr;
    left = N.left;
    right = N.right;
    if (N.originalDataVector != NULL)
    {
        originalDataVector = new float[originalDim];
        for (long i=0;i<originalDim;i++)
        {
            originalDataVector[i] = N.originalDataVector[i];
        }
    }
    else
    {
        originalDataVector = NULL;
    }
}

// Assignment operator
Node & Node::operator=(const Node &second)
{
    for (long i=0;i<dim;i++)
    {
        record[i] = second.record[i];
    }
    dataPtr = second.dataPtr;
    left = second.left;
    right = second.right;
    if (second.originalDataVector != NULL)
    {
        if (originalDataVector == NULL)
        {
            originalDataVector = new float[originalDim];
        }
        for (long i=0;i<originalDim;i++)
        {
            originalDataVector[i] = second.originalDataVector[i];
        }
    }
    else
    {
        if (originalDataVector != NULL)
        {
            // Second is null but first is not
            delete [] originalDataVector;
            originalDataVector = NULL;
        }
    }
    return *this;
}

void printOriginalDataVector(float *originalDataVector,long dim) {
    for ( int i = 0 ; i < dim ; i++ ) 
        printf("%d ",(int)originalDataVector[i]);
    puts("");
}
void Node::rangeQuery(QueryBox *qbox, ResultSet &rs, long desc
, float *center, float radius
)
{
    bool inside = true;
    rs.countVisitedNode ++;
    for(long d=0;d<dim && inside;d++)
    {
        rs.countComparisons+=2;
        if (record[d] < qbox->cntBox[2*d] ||
            record[d] > qbox->cntBox[2*d+1])
        {
            // Not inside
            inside = false;
        }
    }
    if (inside)
    {
        // Confirmatory test.
        // Calculate actual distance
         rs.countDistComputations++;
         /*
                2016 11 02 getL2Distance > getHammingDistance
            */
         //float dist = getL2Distance(center, originalDataVector, originalDim);
         float dist = getL2Distance(center, originalDataVector, originalDim);
        if (dist <= radius)
        {
            rs.countHits++;
    //        printOriginalDataVector(originalDataVector,dim);
//            rs.hitList.push_back(dataPtr);
        //    cout<<dataPtr<<endl;
        }
    }
    rs.countComparisons+=2;

    // Decide whether to access left and right subtree
    if (left != NULL)
    {
        // Modify the querybox
//        QueryBox left_box(qbox);
//        left_box.cntBox[2*desc + 1] = record[desc];
        
        // check the decriminator dimension
        // at least some part of the box must lie in the left
        // subspace of current node. This is possible only if
        // value of this node on desc dimension is greater than lower
        // bound of the box
        if (record[desc] >= qbox->cntBox[2*desc])
        {
            // search left subtree
            left->rangeQuery(qbox, rs, ((desc + 1)%dim)
            , center, radius
            );
        }
    }
    if (right != NULL)
    {
        // Modify the querybox
//        QueryBox right_box(qbox);
//        right_box.cntBox[2*desc] = record[desc];

        if (record[desc] <= qbox->cntBox[2*desc+1])
        {
            // search right subtree
            right->rangeQuery(qbox, rs, ((desc+1)%dim)
            , center, radius
            );
        }
    }
}

void Node::knnQuery(float *center, float &radius, ResultSet &rs, priority_queue  <  queueItem > *nodeQ, float * boundingBox, const long K, priority_queue < pair<void *, float > , vector < pair < void *, float > > , class pairComparator >  *resultQ, long desc, float *distArray)
{
    // First do the box query
    bool insideBox = true;
    for (long i=0;i<dim && insideBox ;i++)
    {
        rs.countComparisons+=2;
        if (record[i] < distArray[i] - radius ||
            record[i] > distArray[i] + radius )
        {
            insideBox = false;
        }
    }
    if (insideBox)
    {
        // Actual distance computation
        //rs.countDistComputations += 1;
        rs.transformedSpaceDistComp += 1;
       // float dist = getL2Distance(center, originalDataVector, originalDim);
        float dist = getL2Distance(distArray, record, dim);
        if (dist < radius)
        {
            // Replace one of the existing K nearest neighbors
            pair <Node *, float> nextNode(this, dist);
            resultQ->push(nextNode);
            if (resultQ->size() > K)
            {
                resultQ->pop();
                float origRadius = radius;
                radius =(resultQ->top()).second;
                assert (radius <= origRadius);
            }
        }
    }



    // Now lets figure out if we can prune the left and right
    // subtrees.
    if (left != NULL)
    {
        queueItem leftItem;
        leftItem.nodePtr = left;
        leftItem.desc = (desc+1)%dim;
        for (long i=0;i<2*dim;i++)
        {
            leftItem.box[i] = boundingBox[i];
        }
        leftItem.box[2*desc+1] = record[desc];

/*        for (long i=0;i<dim;i++)
        {
            cout<<leftItem.box[2*i]<<" : "<<leftItem.box[2*i + 1]<<", ";
        }cout<<endl;
        for (long i=0;i<dim;i++)
        {
            cout<<center[i]<<",";
        }
        cout<<endl;
*/
        rs.transformedSpaceDistComp += 1;
        leftItem.distance = getBoxL2Distance(distArray, leftItem.box, dim);
//        cout<<"Dist "<<leftItem.distance<<endl;
        nodeQ->push(leftItem);
    }

    if (right != NULL)
    {
        queueItem rightItem;
        rightItem.nodePtr = right;
        rightItem.desc = (desc+1)%dim;
        for (long i=0;i<2*dim;i++)
        {
            rightItem.box[i] = boundingBox[i];
        }
        rightItem.box[2*desc] = record[desc];
        rs.transformedSpaceDistComp += 1;
        rightItem.distance = getBoxL2Distance(distArray, rightItem.box, dim);
        nodeQ->push(rightItem);
    }
}

