#include "UserInterface.h"

void displayHelp()
{
    cout<<"Below are the important options:"<<endl;
    cout<<OPT_LOAD_FILE<<" : Specifies the data file to load"<<endl;
    cout<<OPT_ORIG_DIM<<" : Specifies the number of dimensions in the data"<<endl;
    cout<<OPT_DIM<<" : Specifies the number of reference points"<<endl;
    cout<<OPT_ORIG_DIM<<" : Specifies the number of data dimensions"<<endl;
    cout<<OPT_RQFILE<<" : Speciies the name of the file containing the range queries"<<endl;
    cout<<OPT_RANGE<<" : Specifies the query radius"<<endl;
    cout<<OPT_KNNFILE" : Specifies the name of the Knn query file"<<endl;
    cout<<OPT_KNN<<" : Specifies the number of neighbors (K) for knn query."<<endl;
    cout<<OPT_REFMETHOD<<" : Specifies the method to be used for reference point selection"<<endl;
    cout<<"\t\tValid values are :"<<endl;
    cout<<"\t\t"<<VELTKAMP<<"\n\t\t"<<ANGLE_DEV<<"\n\t\t"<<RANDOM<<"\n\t\t"<<RANDOMCORNER<<"\n\t\t"<<MVP<<endl;
    cout<<OPT_COUNT<<" : Specifies the number of data points to index"<<endl;
    cout<<OPT_BEGIN<<" : Specifies the number of lines in the data file to skip"<<endl;
    cout<<OPT_HELP" : Displays this help"<<endl;
}

// Helper function - calculates std deviation of spacing
float spacingVariance(float* distances, long numDist)
{
    vector <float> distanceVector (numDist);
    for (long i=0; i<numDist; i++)
    {
        distanceVector[i] = distances[i];
    }
    sort(distanceVector.begin(), distanceVector.end());
    accumulator_set < float , stats <tag::variance(lazy) > > spacingAcc;
    for (long i=1;i<numDist;i++)
    {
        // Add the new spacing value in the accumulator
        spacingAcc(distanceVector[i] - distanceVector[i-1]);
    }
    return variance(spacingAcc);
}

// Calculates correlation of two arrays and returns false if
// it is greater than the threshold. Returns true otherwise.
bool correlationNotOk(float * distVector1, float* distVector2, long length)
{
    // First we need to calculate linear correlation coefficient.
    float d1 = 0., d2 =0., d1d2= 0., d1square = 0., d2square = 0.;

    for (long i=0;i<length;i++)
    {
        d1 += (distVector1[i]);
        d2 += (distVector2[i]);
        d1square += (distVector1[i]*distVector1[i]);
        d2square += (distVector2[i]*distVector2[i]);
        d1d2 += (distVector1[i]*distVector2[i]);
    }

    float corr = (length*(d1d2) - d1*d2)/
        (  sqrt( length*d1square - d1*d1 ) *
           sqrt( length*d2square - d2*d2 ) );
    if (fabs(corr) > CORR_THR)
    {
        // Correlation is NOT OK
        return true;
    }
    else
    {
        // Correlation is OK
        return false;
    }

}

// Returns variance of pairwise distances
float getDistVariance(float **refPoints, long numRef, long dim)
{
    accumulator_set < float , stats <tag::variance(lazy) > > allDistances;

    for (long i=0;i<numRef;i++)
    {
        for(long j=i+1;j<numRef;j++)
        {
            float dist = getL2Distance(refPoints[i], refPoints[j], dim);
            allDistances(dist);
        }
    }

    return variance(allDistances);
}

// prints mean and variance of pairwise distances among reference
// points
void printStats(float **refPoints, long numRef, long dim)
{
    accumulator_set < float , stats <tag::variance(lazy) > > allDistances;
    long *countPairs = new long[dim+1];
    vector < pair< long, long> > at4;
    vector < pair <long, long> > at5;

    for (long d=0;d<=dim;d++)
        countPairs[d] = 0;
    for (long i=0;i<numRef;i++)
    {
        for(long j=i+1;j<numRef;j++)
        {
            float dist = getL2Distance(refPoints[i], refPoints[j], dim);
            countPairs[lround(dist*dist)] += 1;
            if (lround(dist*dist) == 4)
                at4.push_back(pair <long, long> (i, j));
            if (lround(dist*dist) == 5)
                at5.push_back(pair <long, long> (i, j));
            allDistances(dist);
        }
    }/*
    cout<<"AT 4 "<<endl;
    for (long i=0;i<at4.size();i++)
    {
        cout<<at4[i].first<<", "<<at4[i].second<<endl;
    }
    cout<<endl<<"AT 5 "<<endl;
    for (long j=0;j<at5.size();j++)
    {
        cout<<at5[j].first<<", "<<at5[j].second<<endl;
    }*/
    cout<<"Pairwise distance mean "<<mean(allDistances)<<endl;
    cout<<"Pairwise distance variance "<<variance(allDistances)<<endl;
    cout<<"Number of pairs at given distance "<<endl;
    for(long d=0;d<=dim;d++)
    {
        cout<<d<<" : "<<countPairs[d]<<endl;
    }
    cout<<endl;
}

// Generate reference points using Veltkamp's method
void generateVeltkampRefPoints(float **refPoints, long numRef, float ** dataPoints, long numData, long dim)
{
    // Matrix containing distances from vantage points
    float **distanceMatrix = new float*[numRef];
    for (long rf=0;rf < numRef;rf++)
    {
        distanceMatrix[rf] = new float[numData];
    }

    // Initial candidates are random data points
    generateRandomRefPoints(refPoints, numRef, dataPoints, numData, dim);

    // Calculate distance matrix
    for (long rf=0;rf < numRef;rf++)
    {
        for (long i=0;i<numData;i++)
        {
            distanceMatrix[rf][i] = getL2Distance(dataPoints[i], refPoints[rf], dim);
        }
    }

    // Spacing Heuristic
    for (long rf=0;rf<numRef;rf++)
    {
        float spacingVar = spacingVariance(distanceMatrix[rf], numData);
        while (spacingVar > SPACING_THR)
        {
            // Select a new point as the reference point
            long next = rand()%numData;
            // Copy the "next" point
            for (long d = 0;d < dim;d++)
            {
                refPoints[rf][d] = dataPoints[next][d];
            }
            // Update the distance matrix
            for (long i=0;i<numData;i++)
            {
                distanceMatrix[rf][i] = getL2Distance(dataPoints[i], refPoints[rf], dim);
            }
           spacingVar = spacingVariance(distanceMatrix[rf], numData);
        }   
    }

    // Correlation Heuristic
    bool change = true;
    while (change)
    {
        change = false;
        // No point running the loops if anything changes.
        // We need to restart in such a case
        for (long i=0;i<numRef && !change ;i++)
        {
            for (long j=i+1;j<numRef && !change ;j++)
            {
                // Calculate correlation of i-th and
                // j-th distance vector.
                if (correlationNotOk(distanceMatrix[i], distanceMatrix[j], numData))
                {
                    float V1, V2;
                    V1 = spacingVariance(distanceMatrix[i], numData);
                    V2 = spacingVariance(distanceMatrix[j], numData);
                    long rfToChange = i;
                    if (V2 > V1)
                    {
                        // j-th point has higher variance
                        rfToChange = j;
                    }
                    bool spacingNotOk = true;
                    while (spacingNotOk)
                    {
                        // Select a new point as the reference point
                        long next = rand()%numData;
                        // Copy the "next" point
                        for (long d = 0;d < dim;d++)
                        {
                            refPoints[rfToChange][d] = dataPoints[next][d];
                        }
                        // Update the distance matrix
                        for (long i=0;i<numData;i++)
                        {
                            distanceMatrix[rfToChange][i] = getL2Distance(dataPoints[i], refPoints[rfToChange], dim);
                        }
                        if (spacingVariance(distanceMatrix[rfToChange], numData) < SPACING_THR)
                        {
                            // Spacing is ok
                            spacingNotOk = false;
                        }
                    }
                    change = true;
                }
            }
        }
    }

    // Delete Distance matrix
    for (long rf=0; rf < numRef;rf++)
    {
        delete[] distanceMatrix[rf];
    }
    delete[] distanceMatrix;
}

// Generate random data points as reference points
void generateRandomRefPoints(float **refPoints, long numRef, float ** dataPoints, long numData, long dim)
{
    srand(time(NULL));
    // Map of points already used
    map <long, bool> alreadySelected;
    for (long i=0;i<numRef;i++)
    {
        long next = rand() % numData;
        while (alreadySelected.find(next) != alreadySelected.end())
        {
            next = rand() % numData;
        }
        alreadySelected[next] = true;
        // Copy the "next" point
        for (long d = 0;d < dim;d++)
        {
            refPoints[i][d] = dataPoints[next][d];
        }
    }
}

// Another helper function. First calculates all possible
// angles among the reference points. Then calculates
// the variance of those angles.
float getAngleVariance(float **setOfPoints, long numRef, long dim)
{
    accumulator_set < float , stats <tag::variance(lazy) > > allAngles;
    //float fancyMean = 0.;//*M_PI/180.;
    for (long i=0;i<numRef;i++)
    {
        for(long j=i+1;j<numRef;j++)
        {
            for(long k=0;k<numRef;k++)
            {
                if (k != i && k != j)
                {
                    // Calculate angle. Put it in accumulator
                    float D_ij = getL2Distance(setOfPoints[i], setOfPoints[j], dim);
                    float D_jk = getL2Distance(setOfPoints[j], setOfPoints[k], dim);
                    float D_ki = getL2Distance(setOfPoints[k], setOfPoints[i], dim);
                    // Use Cosine rule to get the angle.
                    float cosAlpha = (D_jk*D_jk + D_ki*D_ki - D_ij*D_ij)/(2.*D_jk*D_ki);
                    //cout<<D_ij<<"  :  "<<D_jk<<"  :  "<<D_ki<<"  :  "<<cosAlpha<<endl;
                    allAngles(acos(cosAlpha));
                }
            }
        }
    }
    return variance(allAngles);
}

// Generate reference points using data independent method
// Uses repeated random samples to get the final vantage point set.
void generateDataIndependentRefPointsRandomized(float **refPoints, long numRef, long dim)
{
    // Temporary set of candidate reference points.
    float ** candidateRefPoints;
    candidateRefPoints = new float* [numRef];
    for (long i=0;i<numRef;i++)
    {
        candidateRefPoints[i] = new float[dim];
    }

    // We have to minimize variance of angles
    // among the reference points.
    float minAngleVar = 999999.;
    for (long iter = 0;iter < ITERATIONS;iter++)
    {
        // Generate a new sample of candidate reference
        // points
        generateRandomCornerRefPoints(candidateRefPoints, numRef, dim);
        // Calculate angle variance
        float angleVar = getAngleVariance(candidateRefPoints, numRef, dim);
        if (angleVar < minAngleVar)
        {
            minAngleVar = angleVar;
            // Copy the candidate ref points to the main set
            for(long i=0;i<numRef;i++)
            {
                for (long j=0;j<dim;j++)
                {
                    refPoints[i][j] = candidateRefPoints[i][j];
                }
            }
        }
    }

    // Delete the temporary set of candidate referenec points.
    for (long i=0;i<numRef;i++)
    {
        delete [] candidateRefPoints[i];
    }
    delete [] candidateRefPoints;

}

// A recursive function to generate strings of 0,1
// used as an intermediate function for reference point generator.
vector <string> recursiveGenerator(long dim)
{
    // Number of dimensions MUST be a power of 2.
    if (dim == 2)
    {
        vector <string> result;
        result.push_back("00");
        result.push_back("01");
        return result;
    }
    vector <string> tempResult = recursiveGenerator(dim/2);
    vector <string> result;
    for (long i=0;i<dim/2;i++)
    {
        result.push_back(tempResult[i] + tempResult[i]);
        string temp = tempResult[i];
        for (long j=0;j< dim/2;j++)
        {
            if (temp[j] == '0')
            {
                temp[j] = '1';
            }
            else
            {
                temp[j] = '0';
            }
        }
        result.push_back(temp + tempResult[i]);
    }
    return result;

}


bool pairComparator (pair <long, long> p1, pair <long, long > p2)
{
    return p1.second < p2.second;
}

// Generate reference points using data independent method
void generateDataIndependentRefPointsRecursive(float **refPoints, long numRef, long dim)
{
    if (numRef > dim)
    {
        cerr<<"Cannot generate more than "<<dim<<" points."<<endl;
        exit(1);
    }
    // Find the nearest power of 2
    float log2 = log(float(dim))/log(float(2.));
    long dimPower2 = dim;

    if ( fabs(log2 - int(log2)) > 0.00001 )
    {
        // The number is not a power of two.
        dimPower2 = long(pow(2., 1+int(log2)));
    }


    float ** superRefPoints = new float * [dimPower2];
    for (long i=0;i<dimPower2;i++)
    {
        superRefPoints[i] = new float [dimPower2];
    }
    
    vector <string> pointString = recursiveGenerator(dimPower2);




    for (long i=0;i<dimPower2;i++)
    {
        for (long j=0;j<dimPower2;j++)
        {
            if (pointString[i][j] == '0')
            {
                superRefPoints[i][j] = 0.;
            }
            else
            {
                superRefPoints[i][j] = 1.;
            }
        }
    }
    // Now eliminate columns
    bool even = true;
    bool* eliminated = new bool [dimPower2];
    for (long c=0;c<dimPower2;c++)
    {
        eliminated[c] = false;
    }
    float ** candidate = new float * [dimPower2];
    for (long col=0;col < dimPower2 - dim ; col++)
    {
        // Eliminate a column from second half
        long begin = dimPower2/2;
        if (!even)
        {
            // Eliminate a column from first half
            begin = 0;
        }

        // Allocate space for candidates
        for (long rf=0;rf<dimPower2;rf++)
        {
            candidate[rf] = new float[dimPower2 - col -1];
        }

        float minVariance = 99999.;
        long nextColumn = -1;
        for (long c = begin;c < (begin + dimPower2/2);c++)
        {
            if (!eliminated[c])
            {
            // Now copy stuff
            for (long rf=0;rf<dimPower2;rf++)
            {
                for (long i=0, nextC =0 ;i<dimPower2;i++)
                {
                    if (i!= c && !eliminated[i])
                    {
                        // This column is neither eliminated nor the current elimination candidate
                        candidate[rf][nextC++] = superRefPoints[rf][i];
                    }
                }
            }
            // Ahh... Finally we have the next candidate set
            // Get the distance variance
            float distVar = getDistVariance(candidate, dimPower2, (dimPower2 - col -1));
            if (distVar < minVariance )
            {
                minVariance = distVar;
                nextColumn = c;
            }
            }
        }

        // Deallocate space
        for (long rf=0;rf<dimPower2;rf++)
        {
            delete [] candidate[rf];
        }
        // Eliminate next column
        assert(nextColumn >= 0);
        eliminated[nextColumn] = true;
        even = !even;
    }
    float ** level1RefPoints = new float* [dimPower2];
    for (long i=0;i<dimPower2;i++)
    {
        level1RefPoints[i] =  new float [dim];
    }

    for (long i=0;i<dimPower2;i++)
    {
        for (long j =0, next=0;j<dimPower2;j++)
        {
            if (!eliminated[j])
            {
                level1RefPoints[i][next++] = superRefPoints[i][j];
            }
        }
    }

    bool * repeated = new bool[dimPower2];
    for (long i=0;i<dimPower2;i++)
    {
        repeated[i] = false;
    }


    // Now see if there are any rows that are repeated.
    for (long i=0;i<dimPower2;i++)
    {
        for (long j =i+1;j<dimPower2;j++)
        {
            // If distance between two points is too small, one of them
            // is flagged as "repeated"
            if (getL2Distance (level1RefPoints[i], level1RefPoints[j], dim) < 0.0001 )
            {
                repeated[j] = true;
            }
        }
    }

    vector < pair <long, long> > countFlips;
    long totalNonRepeatedRef = 0;
    for (long i=0;i<dim;i++)
    {
        if (!repeated[i])
        {
            pair <long, long> foo (i, 0);
            countFlips.push_back(foo);
            totalNonRepeatedRef++;
        }
    }
    for (long i=0;i<totalNonRepeatedRef;i++)
    {
        long flips = 0;
        for (long j=1;j<dim;j++)
        {
            if (level1RefPoints[countFlips[i].first][j] != level1RefPoints[countFlips[i].first][j-1])
            {
                flips += 1;
            }
        }
        countFlips[i].second = flips;
    }

    // Now sort the dimensions in asceneding
    // order of number of flips
    sort(countFlips.begin(), countFlips.end(), pairComparator);
    if (numRef > totalNonRepeatedRef)
    {
        cerr<<"Unable to generate required number of unique points"<<endl;
        exit(1);
    }

    for (long i=0;i< numRef;i++)
    {
        for (long j=0;j<dim;j++)
        {
            refPoints[i][j] = level1RefPoints[countFlips[totalNonRepeatedRef - 1 -i].first ][j];
            //refPoints[i][j] = level1RefPoints[countFlips[i].first ][j];
        }
    }
/*
    for (long i=0;i< numRef;i++)
    {
        for (long j=0;j<dim;j++)
        {
            refPoints[i][j] = level1RefPoints[i][j];
        }
    }*/
    for (long i=0;i<dimPower2;i++)
    {
        delete [] level1RefPoints[i];
    }
    delete [] level1RefPoints;
    
    // Free up "super" refpoints
    for (long i=0;i<dimPower2;i++)
    {
        delete [] superRefPoints[i];
    }
    delete [] superRefPoints;

    delete [] eliminated;
    delete [] repeated;
}

// Generate reference points using data independent method
void generateDataIndependentRefPoints(float **refPoints, long numRef, long dim)
{
    // We have to minimize variance of distances
    // among the reference points.

    // Genearate all corner points
    long numCorners = pow(2, dim);
    float ** allCorners = new float*[numCorners];
    for (long i=0;i<numCorners;i++)
    {
        // Convert i to binary and express that as an array of floats
        float *nextPoint = new float[dim];
        long mask = 1;
        for (long bit = 0;bit < dim;bit++, mask = mask << 1)
        {
            nextPoint[bit] = (i & mask)? 1.0 : 0.0;
        }
        allCorners[i] = nextPoint;
    }
    map <long, bool> selectedPoints;
    // Select first two points randomly
    long point = rand() % numCorners;
    point = 0;
    for (long d=0;d<dim;d++)
    {
        refPoints[0][d] = allCorners[point][d];
        // We copy the same point in next location as well
        // because the next point shares 50% of the bits.
        if (numRef > 1)
        {
            refPoints[1][d] = allCorners[point][d];
        }
    }
    selectedPoints[point] = true;

    if (numRef > 1)
    {
        // Next point is obtained by flipping 50% of the bits
        for(long d = 0;d<(dim)/2 ; d++)
        {
            // Flip the bit
            refPoints[1][d] = 1. - refPoints[1][d];
        }

        // Calculate corresponding point
        point = 0;
        long adder = 1;
        for (long d=0;d<dim;d++)
        {
            if (refPoints[1][d] > 0.)
            {
                point += adder;
            }
            adder += adder;
        }

        selectedPoints[point] = true;
    }
    // Now the first two points are "seeded" properly
    // Now select rest of the point using greedy algorithm
    for (long nextRef = 2;nextRef< numRef;nextRef++)
    {
        float minDistVar = 999999.;
//        float maxAngleVar = -1.;
        long next = -1;
        for (long i=0;i<numCorners;i++)
        {
            // continue only if this point is not already included.
            if (selectedPoints.find(i) == selectedPoints.end())
            {
                // Temporarily add this point
                for (long d=0;d<dim;d++)
                {
                    refPoints[nextRef][d] = allCorners[i][d];
                }
                float distVar = getDistVariance(refPoints, nextRef+1, dim);
                if (distVar < minDistVar)
                {
                    // Finalize this point
                    next = i;
                    minDistVar = distVar;
                }

            }
        }
        // Finalize "next" point
        for (long d=0;d<dim;d++)
        {
            refPoints[nextRef][d] = allCorners[next][d];
        }
    }

    // Delete all the corners
    for (long i=0;i<numCorners;i++)
    {
        delete[] allCorners[i];
    }
    delete [] allCorners;

}

// Returns true if the set of point contains a repeated point.
bool containsRepeatedPoint(float **refPoints, long numRef, long dim)
{
    // The basic idea is very simple. If there are repeated points then
    // the minimum pairwise distance will be zero. To make matters easy,
    // we can use L1 distance instead of L2
    float minDist = 99999.;
    for (long i=0;i<numRef;i++)
    {
        for (long j=i+1;j<numRef;j++)
        {
            float dist = 0.;
            for (long d=0;d<dim;d++)
            {
                dist += fabs(refPoints[i][d] - refPoints[j][d]);
            }
            if (dist < minDist)
            {
                minDist = dist;
            }
        }
    }
    if (minDist < 0.0001)
    {
        return true;
    }
    else
    {
        return false;
    }
}


// Generate reference points using data independent method
void generateDataIndependentRefPointsHighDim(float **refPoints, long numRef, long dim)
{
    // We have to minimize variance of distances
    // among the reference points.

    long long numCorners = pow((long double)2, dim);
    cout<<"Num Corners "<<numCorners;
    // First point is origin
    for (long d=0;d<dim;d++)
    {
        refPoints[0][d] = 0.;
        // The next point shares 50% of the bits.
        if (d < (dim+1)/2)
        {
            refPoints[1][d] = 1.;
        }
        else
        {
            refPoints[1][d] = 0.;
        }
    }
    float *next = new float[dim];
    // Now the first two points are "seeded" properly
    // Now select rest of the point using greedy algorithm
    for (long nextRef = 2;nextRef< numRef;nextRef++)
    {
        float minDistVar = 999999.;
        for (long long i=0;i<numCorners;i++)
        {
                // Temporarily add this point
                long long mask = 1;
                for (long d=0;d<dim;d++, mask = mask<<1)
                {
                    if (mask & i)
                    {
                        refPoints[nextRef][d] = 1.;
                    }
                    else
                    {
                        refPoints[nextRef][d] = 0.;
                    }
                }

            if (! containsRepeatedPoint(refPoints, numRef, dim))
            {
                // Evaluate angleVariance of new set
                float distVar = getDistVariance(refPoints, nextRef+1, dim);
                if (distVar < minDistVar)
                {
                    // Finalize this point
                    for (long d=0;d<dim;d++)
                    {
                        next[d] = refPoints[nextRef][d];
                    }
                    minDistVar = distVar;
                }

            }
        }
        // Finalize "next" point
        for (long d=0;d<dim;d++)
        {
            refPoints[nextRef][d] = next[d];
        }
    }
    delete[] next;
}

// Generate reference points using data independent method
void generateDataIndependentRefPoints2(float **refPoints, long numRef, long dim)
{
    // Temporary set of candidate reference points.
    float ** candidateRefPoints;
    candidateRefPoints = new float* [numRef];
    for (long i=0;i<numRef;i++)
    {
        candidateRefPoints[i] = new float[dim];
    }

    // We have to minimize variance of angles
    // among the reference points.

    // Genearate all corner points
    long numCorners = pow(2, dim);
    float ** allCorners = new float*[numCorners];
    for (long i=0;i<numCorners;i++)
    {
        // Convert i to binary and express that as an array of floats
        float *nextPoint = new float[dim];
        long mask = 1;
        for (long bit = 0;bit < dim;bit++, mask = mask << 1)
        {
            nextPoint[bit] = (i & mask)? 1.0 : 0.0;
        }
        allCorners[i] = nextPoint;
    }
    map <long, bool> selectedPoints;
    // Select first three points using exhaustive search
    float minAngleVar = 99999.;
    long ii = -1 ,jj = -1,kk = -1;
    for (long i=0;i<numCorners;i++)
    {
        for (long j = i+1;j<numCorners;j++)
        {
            for(long k = j+1;k<numCorners;k++)
            {
                // Copy i-th, j-th and k-th point as candidates
                for (long d=0;d<dim;d++)
                {
                    candidateRefPoints[0][d] = allCorners[i][d];
                    candidateRefPoints[1][d] = allCorners[j][d];
                    candidateRefPoints[2][d] = allCorners[k][d];
                }
                float angleVar = getAngleVariance(candidateRefPoints, 3, dim);
                if (angleVar < minAngleVar ||
                    (fabs(angleVar - minAngleVar) < 0.001 && 
                    // Accept this point with some randomness
                    float(rand())/float(RAND_MAX) > 0.5 ) )
                {
                    minAngleVar = angleVar;
                    for (long d=0;d<dim;d++)
                    {
                        refPoints[0][d] = candidateRefPoints[0][d];
                        refPoints[1][d] = candidateRefPoints[1][d];
                        refPoints[2][d] = candidateRefPoints[2][d];
                    }
                    ii = i; jj = j; kk = k;
                }
            }
        }
    }

    selectedPoints[ii] = true;
    selectedPoints[jj] = true;
    selectedPoints[kk] = true;

    // Now select rest of the point using greedy algorithm
    for (long nextRef = 3;nextRef< numRef;nextRef++)
    {
        minAngleVar = 999999.;
        long next = -1;
        for (long i=0;i<numCorners;i++)
        {
            // continue only if this point is not already included.
            if (selectedPoints.find(i) == selectedPoints.end())
            {
                // Temporarily add this point
                for (long d=0;d<dim;d++)
                {
                    refPoints[nextRef][d] = allCorners[i][d];
                }
                // Evaluate angleVariance of new set
                float angleVar = getAngleVariance(refPoints, nextRef + 1, dim);
                if (angleVar <minAngleVar)
                {
                    // Finalize this point
                    next = i;
                    minAngleVar = angleVar;
                }

            }
        }
        // Finalize "next" point
        for (long d=0;d<dim;d++)
        {
            refPoints[nextRef][d] = allCorners[next][d];
        }
    }

    // Delete the temporary set of candidate referenec points.
    for (long i=0;i<numRef;i++)
    {
        delete [] candidateRefPoints[i];
    }
    delete [] candidateRefPoints;

    // Delete all the corners
    for (long i=0;i<numCorners;i++)
    {
        delete[] allCorners[i];
    }
    delete [] allCorners;

}

// Generate random corner points as reference points
void generateRandomCornerRefPoints(float **refPoints, long numRef, long dim)
{
    // Select a random corner point.
    // Map of points already used
    map <long, bool> alreadySelected;
    for (long rf=0; rf < numRef;rf++)
    {
        bool nextPointNotFound = true;
        while (nextPointNotFound)
        {
            long next = 0;
            // adder is 2^d
            long adder = 1;
            for (long d=0;d<dim;d++)
            {
                float random = (float)rand()/(float)RAND_MAX;
                if (random >= 0.5)
                {
                    refPoints[rf][d] = 1.;
                    next += adder;
                }
                else
                {
                    refPoints[rf][d] = 0.;
                }
                adder += adder;
            }
            if (alreadySelected.find(next) == alreadySelected.end())
            {
                // Not found. Hence a new point
                alreadySelected[next] = true;
                nextPointNotFound = false;
            }
        }
    }
}

// Function to read command line options
void getOptions(int argCount, char **argVector, options * opt)
{
    string optDim(OPT_DIM);
    string optOrigDim(OPT_ORIG_DIM);
    string optLoadFile(OPT_LOAD_FILE);
    string optBqfile(OPT_BQFILE);
    string optRqfile(OPT_RQFILE);
    string optReffile(OPT_REF_FILE);
    string optKnnfile(OPT_KNNFILE);
    string optKnn(OPT_KNN);
    string optRange(OPT_RANGE);
    string optSkip(OPT_BEGIN);
    string optCount(OPT_COUNT);
    string optStats(OPT_STATS);
    string optHelp(OPT_HELP);
    string optRefMethod(OPT_REFMETHOD);

    for(long i=0;i<argCount;i++)
    {
        if(!optBqfile.compare(argVector[i]))
        {
            opt->bqfile = argVector[i+1];
            i++;
        }
        if(!optRqfile.compare(argVector[i]))
        {
            opt->rqfile = argVector[i+1];
            i++;
        }
        if(!optReffile.compare(argVector[i]))
        {
            opt->reffile = argVector[i+1];
            i++;
        }
        if(!optKnnfile.compare(argVector[i]))
        {
            opt->knnfile = argVector[i+1];
            i++;
        }
        if(!optKnn.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>(opt->knn);
        }
        if(!optLoadFile.compare(argVector[i]))
        {
            opt->datafile = argVector[i+1];
            i++;
        }
        if(!optOrigDim.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>opt->origDim;
        }
        if(!optDim.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>opt->dim;
        }
        if(!optRange.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>opt->range;
        }
        if(!optCount.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>opt->count;
        }
        if(!optSkip.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>opt->skip;
        }
        if (!optRefMethod.compare(argVector[i]))
        {
            opt->refMethod = argVector[i+1];
            i++;
        }
        if(!optStats.compare(argVector[i]))
        {
            opt->stats = true;
        }
        if(!optHelp.compare(argVector[i]))
        {
            opt->help = true;
        }
    }
}

int main(int argc, char ** argv)
{
    srand ( time(NULL) );
    options opt;
    //Initialize options to Undef values
    opt.datafile = UNDEF_STR;
    opt.bqfile = UNDEF_STR;
    opt.rqfile = UNDEF_STR;
    opt.knnfile = UNDEF_STR;
    opt.reffile = UNDEF_STR;
    // Default method for reference point generation is
    // to use random data points as reference points.
    opt.refMethod = RANDOM;
    //opt.knnfile = UNDEF_STR;
    //opt.knn = UNDEF_LONG;
    opt.range = 0.;
    opt.origDim = UNDEF_LONG;
    opt.dim = UNDEF_LONG;
    opt.skip = 0;
    opt.knn=0;
    opt.count = 99999999;
    opt.help = false;
    opt.stats = false;
    getOptions(argc, argv, &opt);

    if (opt.help)
    {
        displayHelp();
        return 0;
    }

    if (argc < 2)
    {
        cerr<<"Use -help option to get usage help"<<endl;
        return 1331;
    }

    if (!opt.datafile.compare(UNDEF_STR))
    {
        cerr<<"Specify data file to load using the mandatory "<<OPT_LOAD_FILE<<" option"<<endl;
        return 1331;
    }

/*    if (opt.origDim > 63)
    {
        cerr<<"This program does not work for dimensions > 63"<<endl;
        return 1331;
    }*/
/*
    if (opt.origDim > 24)
    {
        cerr<<"The number of dimensions is too high. I cannot even guarantee correctness !!!"<<endl;
        return 1729;
    }
*/
    // Create an empty KD tree with given dimensions.
    float constructionTime = 0.;
    clock_t start = clock();
    KDTree kdTree(opt.dim, opt.origDim);
    clock_t end = clock();
    constructionTime += end-start;

    ifstream inFile;
    inFile.open(opt.datafile.c_str(), ios::in);

    string line;
    // now read the actual data.

    if(opt.skip > 0)
    {
        // Skip intial lines
        for(long i=0;i<opt.skip && !inFile.eof();i++)
        {
            getline(inFile, line);
        }
    }
    long countInserted = opt.skip;
    float *distanceBuffer = new float[opt.dim];
    // We are going to read all the points
    if (opt.count > 1000000)
    {
        cerr<<" The number of points is too big. There may be problems with the memory management"<<endl;
    }

    // TODO Add exception handling to protect against possible
    // allocation issues.
    float ** allPoints= new (nothrow) float*[opt.count];
    if (!allPoints)
    {
        cerr<<"Could not allocate space"<<endl;
        exit(1);
    }
    for (int i=0;i<opt.count;i++)
    {
        // Allocate
        allPoints[i] = new (nothrow) float[opt.origDim];
        if (!allPoints[i])
        {
            cerr<<"Could not allocate space"<<endl;
            exit(1);
        }
    }

    for(int i=0;i< opt.count && !inFile.eof(); i++)
    {
        getline(inFile, line);
        if(line[0] != '#' && line[0] != '\n')
        {
            // This not the first non-comment line so it contains vectors
            if(line.length() > 3) // vector must have at least 3 characters
            {
                stringstream input(line);
                for(long d =0;d<opt.origDim;d++)
                {
                    input>>allPoints[i][d];
                }
            }
        }
    }
    // File reading complete

    // Create reference points
    float ** refPoints = new (nothrow) float *[opt.dim];
    for (long rf=0; rf < opt.dim;rf++)
    {
        refPoints[rf] = new float[opt.origDim];
    }

    // Generate reference point. Although the program allows
    // reading reference points from a file (mainly for
    // historical reasons), almost always, we will generate
    // them on the fly using heuristics.

    clock_t ref_points_start_time = clock();
    if (opt.reffile.compare(UNDEF_STR))
    {
        ifstream refFile;
        // User did pass some value for ref file
        refFile.open(opt.reffile.c_str(), ios::in);
        for (long rf=0; rf < opt.dim;rf++)
        {
            getline(refFile, line);
            stringstream input(line);
            for(long d =0;d<opt.origDim;d++)
            {
                input>>refPoints[rf][d];
            }
        }
        refFile.close();
    }
    else
    {
        // User did not pass a value for ref file
        // Use the heuristics

        if (!opt.refMethod.compare(VELTKAMP))
        {
            generateVeltkampRefPoints(refPoints, opt.dim, allPoints, opt.count, opt.origDim);
        }
        else if (!opt.refMethod.compare(RANDOMCORNER))
        {
            generateRandomCornerRefPoints(refPoints, opt.dim, opt.origDim);
        }
        else if (!opt.refMethod.compare(RANDOM))
        {
            generateRandomRefPoints(refPoints, opt.dim, allPoints, opt.count, opt.origDim);
        }
        else if (!opt.refMethod.compare(ANGLE_DEV))
        {/*
            if (opt.origDim <33 && opt.origDim != 16 && opt.origDim != 8 && opt.origDim != 32)
            {
                generateDataIndependentRefPoints(refPoints, opt.dim, opt.origDim);
            }
            else
            {
                if (opt.origDim == 8 || opt.origDim == 16 || opt.origDim == 32 || opt.origDim == 64 || opt.origDim == 128)
                {*/
                   // generateDataIndependentRefPointsRecursive(refPoints, opt.dim, opt.origDim);
                    generateDataIndependentRefPoints(refPoints, opt.dim, opt.origDim);
               /* }
                else
                {
                    generateDataIndependentRefPointsRandomized(refPoints, opt.dim, opt.origDim);
                }
            }*/
        }
        else if (!opt.refMethod.compare(MVP))
        {
            ///////////// MVP Tree Heuristic /////////////
            for (long rf=0; rf < opt.dim;rf++)
            {
                int randomIndex = rand() % opt.count;
                float maxDist = 0.;
                long maxIndex = -1;
                for (long j=0;j<opt.count;j++)
                {
                    // get distance between elements at randomIndex and at j
                    float dist = getL2Distance(allPoints[randomIndex], allPoints[j], opt.origDim);
                    if (dist > maxDist)
                    {
                        maxDist = dist;
                        maxIndex = j;
                    }
                }

                for (long d=0;d<opt.origDim;d++)
                {
                    refPoints[rf][d] = allPoints[maxIndex][d];
                }
            }
            ///////////////////////////////////////////
        }
        else
        {
            cerr<<"Incorrect reference point method spcified !"<<endl;
            exit(2);
        }
    }

    // We just need to collect statistics
/*    cout<<"Method  : "<<opt.refMethod<<endl;
    cout<<"Point                    #Ones"<<endl;
    for(long i=0;i<opt.dim;i++)
    {
        long countOnes = 0;
        for (long j=0;j<opt.origDim;j++)
        {
            cout<<lround(refPoints[i][j])<<",";
            countOnes += lround(refPoints[i][j]);
        }
        cout<<"        "<<countOnes<<endl;
    }
    printStats(refPoints, opt.dim, opt.origDim);
*/
    clock_t ref_points_end_time = clock();
    float refPointCreationTime = float(ref_points_end_time - ref_points_start_time)/float(CLOCKS_PER_SEC);
//    cout<<"DEBUG  REF_PT_TIME,"<<refPointCreationTime<<endl;
//    exit(1);
    // Reference point reading complete
        /*for (long rf=0; rf < opt.dim;rf++)
        {
            for(long d=0;d<opt.origDim;d++)
            {
                cout<<refPoints[rf][d]<<",";
            }
            cout<<endl;
        }*/
    //cout<<"Angle variance "<<getAngleVariance(refPoints, opt.dim, opt.origDim)<<endl<<endl<<endl;
    // Insertion begins
    /*
     *  kd_tree build sisobus
     */
    start = clock();
    for(int i=0;i< opt.count ; i++)
    {
        // Calculate distances from reference points
     /*
        for (long rf=0;rf < opt.dim;rf++)
        {
            distanceBuffer[rf] = getL2Distance(refPoints[rf], allPoints[i], opt.origDim);
        }
    */
        // Create a KD-tree node.
        kdTree.insert(allPoints[i], allPoints[i], countInserted);
        //kdTree.insert(distanceBuffer, allPoints[i], countInserted);
        countInserted++;
    }
    end = clock();
    constructionTime += (end -start);
    inFile.close();

    // Free the memory
    for (int i=0;i<opt.count;i++)
    {
        delete[] allPoints[i];
    }
    delete[] allPoints;
    // Verify KD-tree by in-order traversal.
    constructionTime /= CLOCKS_PER_SEC;
//    kdTree.printInOrder();
//    cout<<"DEBUG  Insertion complete"<<endl;

    // Range query
    // We do not execute the range query directly.
    // We first convert it to a distance based
    // box query.

    /*
     * range query sisobus
     */
    if (opt.rqfile.compare(UNDEF_STR) && opt.range > 0.)
    {
        // RQ file is defined and the radius is defined
        ifstream rqFile;
        rqFile.open(opt.rqfile.c_str(), ios::in);

        float avgHits = 0.;
        float avgComparisons = 0.;
        float avgDistCalculations = 0.;
        long queryCnt = 0;
        float avgVisitedNode = 0.;
        clock_t totalTime = 0;

        vector <float *> rangeQueries;
        vector <QueryBox *> tRangeQueries;

        getline(rqFile, line);
        while (!rqFile.eof())
        {
            // This not the first non-comment line so it contains vectors
            if(line.length() > 3) // vector must have at least 3 characters
            {
                stringstream input(line);
                float *buffer= new float[opt.origDim];

                for(long d =0;d<opt.origDim;d++)
                {
                    input>>buffer[d];
                }
                // Calculate distances from reference points
                for (long rf=0;rf < opt.dim;rf++)
                {
            //        distanceBuffer[rf] = getL2Distance(refPoints[rf], buffer, opt.origDim);
                      distanceBuffer[rf] = buffer[rf];
                }
                // Convert range query to box query
                QueryBox *qbox = new QueryBox;
                for (long rf=0;rf < opt.dim;rf++)
                {
                    qbox->cntBox[2*rf] = distanceBuffer[rf] - opt.range;
                    qbox->cntBox[2*rf+1] = distanceBuffer[rf] + opt.range;
                }

/*                if (rangeQueries.size() == 0)
                {
                    for (long rf=0;rf < opt.origDim;rf++)
                    {
                        cout<<buffer[rf]<<" , ";
                    }
                    cout<<endl;
                    for (long rf=0;rf < opt.dim;rf++)
                    {
                        cout<<qbox->cntBox[2*rf]<<" "<<qbox->cntBox[2*rf+1]<<",";
                    }
                    cout<<endl;
                }
*/
                // Append the range queries and transformed queries in the list
                // TODO FIX memory leak here... I am not going to free up memory
                // which is a sin.
                rangeQueries.push_back(buffer);
                tRangeQueries.push_back(qbox);
                getline(rqFile, line);
            }
        }

        assert(rangeQueries.size() == tRangeQueries.size());
        start = clock();
        for (long qNum=0;qNum<rangeQueries.size();qNum++)
        {
        //    printf("#%ld\n",qNum);
                // Execute "box query". We pass the original
                // query parameters to be able to eliminate FP.
        /*        if (qNum == 0)
                {
                    for (long rf=0;rf < opt.origDim;rf++)
                    {
                        cout<<rangeQueries[qNum][rf]<<" , ";
                    }
                    cout<<endl;
                    for (long rf=0;rf < opt.dim;rf++)
                    {
                        cout<<tRangeQueries[qNum]->cntBox[2*rf]<<" "<<tRangeQueries[qNum]->cntBox[2*rf+1]<<",";
                    }
                    cout<<endl;
                }
*/
                ResultSet rs = kdTree.rangeQuery(tRangeQueries[qNum], rangeQueries[qNum], opt.range);
                avgHits += rs.countHits;
                avgComparisons += rs.countComparisons;
                avgDistCalculations += rs.countDistComputations;
                queryCnt ++;
                avgVisitedNode += rs.countVisitedNode;

                //cout<<"HITS,"<<rs.countHits<<",COMP,"<<rs.countComparisons<<",DIST,"<<rs.countDistComputations<<endl;
        }
        end = clock();
        totalTime += (end -start);

        rqFile.close();

        avgHits /= queryCnt;
        avgComparisons /= queryCnt;
        avgDistCalculations /= queryCnt;
        avgVisitedNode /= queryCnt;
        printf("AVG_VISITED_NODE%lf\n",avgVisitedNode);
        float avgTime = (float) totalTime / float(CLOCKS_PER_SEC*queryCnt);
        //cout<<"CONST_TIME,"<<constructionTime<<",AVG_HITS,"<<avgHits<<",AVG_COMP,"<<avgComparisons<<",AVG_DIST,"<<avgDistCalculations<<",AVG_TIME,"<<avgTime<<",NumQ,"<<queryCnt<<endl;
        cout<<"REF_PT_TIME,"<<refPointCreationTime<<",INSRT_TIME,"<<constructionTime<<",AVG_QUERY_TIME,"<<avgTime<<",AVG_HITS,"<<avgHits<<",AVG_COMP,"<<avgComparisons<<",AVG_DIST,"<<avgDistCalculations<<",END"<<endl;
    }

                float *buffer= new float[opt.origDim];
//    cout<<"DEBUG About to do knn"<<endl;
    if (opt.knnfile.compare(UNDEF_STR) && opt.knn > 0)
    {
        // KNN file is defined and K is defined
        ifstream knnFile;
        knnFile.open(opt.knnfile.c_str(), ios::in);

        float avgHits = 0.;
        float avgComparisons = 0.;
        float avgDistCalculations = 0.;
        float avgTransformedSpaceDistCalculations = 0.;
        long queryCnt = 0;
        clock_t totalTime = 0;

        getline(knnFile, line);
        start = clock();
        while (!knnFile.eof())
        {
            // This not the first non-comment line so it contains vectors
            if(line.length() > 3) // vector must have at least 3 characters
            {
                stringstream input(line);
                for(long d =0;d<opt.origDim;d++)
                {
                    input>>buffer[d];
                }
                // Calculate distances from reference points
                for (long rf=0;rf < opt.dim;rf++)
                {
                    distanceBuffer[rf] = getL2Distance(refPoints[rf], buffer, opt.origDim);
                }
                ResultSet rs = kdTree.knnQuery(buffer, BLOWUP* opt.knn, distanceBuffer);
                avgHits += rs.countHits;
                avgComparisons += rs.countComparisons;
                avgDistCalculations += rs.countDistComputations;
                avgTransformedSpaceDistCalculations  += rs.transformedSpaceDistComp;
                queryCnt ++;
                cout<<"__________________"<<endl;

                //cout<<"HITS,"<<rs.countHits<<",COMP,"<<rs.countComparisons<<",DIST,"<<rs.countDistComputations<<endl;
            }
            getline(knnFile, line);
        }
        end = clock();
        totalTime += (end -start);

        knnFile.close();

        avgHits /= queryCnt;
        avgComparisons /= queryCnt;
        avgDistCalculations /= queryCnt;
        avgTransformedSpaceDistCalculations /= queryCnt;
        float avgTime = (float) totalTime / float(CLOCKS_PER_SEC*queryCnt);
        //cout<<"CONST_TIME,"<<constructionTime<<",AVG_HITS,"<<avgHits<<",AVG_COMP,"<<avgComparisons<<",AVG_DIST,"<<avgDistCalculations<<",AVG_TIME,"<<avgTime<<",NumQ,"<<queryCnt<<endl;
        cout<<"REF_PT_TIME,"<<refPointCreationTime<<",INSRT_TIME,"<<constructionTime<<",AVG_QUERY_TIME,"<<avgTime<<",AVG_HITS,"<<avgHits<<",AVG_COMP,"<<avgComparisons<<",AVG_DIST,"<<avgDistCalculations<<",AVG_TR_DIST,"<<avgTransformedSpaceDistCalculations<<",END"<<endl;
    }


    //cerr<<" sleeping so that memory snapshot is captured"<<endl;
    //sleep(60);
    delete[] buffer;
    delete[] distanceBuffer;
}

