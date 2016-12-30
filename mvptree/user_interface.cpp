/*

    MVPTree c library 
    Copyright (C) 2008-2009 by D. Grant Starkweather
    All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    D Grant Starkweather - starkd88@gmail.com

*/

/*
MVP tree implementation for in-memory index searching by Alok Watve
This program uses MVP tree library. For questions regarding the library,
contact the original author (D Grant Starkweather). For questions
about this file contact Alok at watvealo@cse.msu.edu
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include "mvptree.h"
#include "user_interface.h"

#define MVP_BRANCHFACTOR 3
#define MVP_LEAFCAP      80
#define MVP_PATHLENGTH   5 

using namespace std;

void testmvpfree(void *ptr){
    free(ptr);
}

// Number of distance computations.
static unsigned long long nbcalcs = 0;


void displayHelp()
{
    cout<<"Sorry. Not yet implemented :-)\nYou are on your own."<<endl;
}

// Function to read command line options
void getOptions(int argCount, char **argVector, options * opt)
{
    string optDim(OPT_DIM);
    string optLoadFile(OPT_LOAD_FILE);
    string optBqfile(OPT_BQFILE);
    string optRqfile(OPT_RQFILE);
    string optRange(OPT_RANGE);
    string optSkip(OPT_BEGIN);
    string optCount(OPT_COUNT);
    string optHelp(OPT_HELP);

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
        if(!optLoadFile.compare(argVector[i]))
        {
            opt->datafile = argVector[i+1];
            i++;
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
        if(!optHelp.compare(argVector[i]))
        {
            opt->help = true;
        }
    }
}

float point_hamming_distance(MVPDP *pointA, MVPDP *pointB){
    if (!pointA || !pointB)
    {
        /*
        cerr<<" Error in distance computation. One of the points is NULL"<<endl;
        // print both the points
        if (pointA)
        {
            cerr<<" A :";
            for (long i=0;i<pointA->datalen;i++)
            {
                cerr<<( (float *)pointA->data)[i]<<" ";
            }
            cerr<<endl;
        }
        if (pointB)
        {
            cerr<<"B :";
            for (long i=0;i<pointB->datalen;i++)
            {
                cerr<<((float *)pointB->data)[i]<<" ";
            }
            cerr<<endl;
        }
        */
        return 10.0;
    } else if (pointA->datalen != pointB->datalen)
    {
        fprintf(stdout," dataA %d, dataB %d\n", pointA->datalen, pointB->datalen);
        return -2.0f;
    }
    int i;
    float sum = 0;
    float *data1 = (float*)pointA->data;
    float *data2 = (float*)pointB->data;
    for (i=0;i<pointA->datalen;i++)
    {
        int d1 = (int)data1[i];
        int d2 = (int)data2[i];
        sum += (d1 != d2);
    }
    nbcalcs++;
    // We do not want to normalize the distance.
    // Hence Alok removed the denominator.
    return sum;
}

float point_l2_distance(MVPDP *pointA, MVPDP *pointB){
    /* L2 distance */
    if (!pointA || !pointB)
    {
        cerr<<" Error in distance computation. One of the points is NULL"<<endl;
        // print both the points
        if (pointA)
        {
            cerr<<" A :";
            for (long i=0;i<pointA->datalen;i++)
            {
                cerr<<( (float *)pointA->data)[i]<<" ";
            }
            cerr<<endl;
        }
        if (pointB)
        {
            cerr<<"B :";
            for (long i=0;i<pointB->datalen;i++)
            {
                cerr<<((float *)pointB->data)[i]<<" ";
            }
            cerr<<endl;
        }
        return 10.0;
    } else if (pointA->datalen != pointB->datalen)
    {
        fprintf(stdout," dataA %d, dataB %d\n", pointA->datalen, pointB->datalen);
        return -2.0f;
    }
    int i;
    float sum = 0;
    float *data1 = (float*)pointA->data;
    float *data2 = (float*)pointB->data;
    for (i=0;i<pointA->datalen;i++)
    {
        float d1 = (float)data1[i];
        float d2 = (float)data2[i];
        float diff = d1 - d2;
        sum += diff*diff;
    }
    nbcalcs++;
    // We do not want to normalize the distance.
    // Hence Alok removed the denominator.
    return sqrt((float)sum);//(float)pointA->datalen;
}


/** Sort of a copy constructor
 * Returns a copy of the point passed as the argument.
 */
MVPDP* mvpdp_dup(MVPDP *dp)
{
    MVPDP *newpnt = dp_alloc(dp->type);
    if (!newpnt) return NULL;
    newpnt->datalen = dp->datalen;
    newpnt->data = malloc(dp->datalen*dp->type);
    memcpy(newpnt->data, dp->data, dp->datalen*dp->type);
//    newpnt->id = strdup(dp->id);
    newpnt->id = dp->id;
    return newpnt;
}

/**
 * Convert a float array in to a MVP data point.
 * This function assumes that sizeof(float) = 4
 * dp_length is the number of dimensions
 */
MVPDP* get_data_point(float * point, unsigned int dp_length, long point_id)
{
    assert(sizeof(float) == 4);
    // Each item in the data is a float which has 32 bits
    MVPDP *newpnt = dp_alloc(UINT32ARRAY);

    if (newpnt == NULL)
    {
        cerr<<"Unable to allocate space"<<endl;
        exit(1);
    }

    newpnt->datalen = dp_length;
    newpnt->data = malloc(dp_length*sizeof(float));

    if (newpnt->data == NULL) {
        free (newpnt);
        cerr<<"Unable to allocate  space for the data point"<<endl;
        exit(1);
    }
    float *row = (float *)newpnt->data;

    for (unsigned int i=0;i<dp_length;i++)
    {
        row[i] = point[i];
    }
    // Copy the ID
    newpnt->id = point_id;
    return newpnt;
}

/**
 * Function that reads all the points from a file.
 * @param nbpoints Maximum number of points to read
 * @param dp_length Number of dimensions
 * @param filename Name of the data file
 */
MVPDP** read_all_points(const unsigned int nbpoints, const unsigned int dp_length, const char * filename)
{
    MVPDP **pointlist = (MVPDP**)malloc(nbpoints*sizeof(MVPDP*));
    if (pointlist == NULL)
    {
        cerr<<"Unable to allocate space for points"<<endl;
        exit(1);
    }

    ifstream inFile;
    inFile.open(filename, ios::in);
    string line;
    long point_id = 0;
    for (unsigned int i=0;i < nbpoints;i++)
    {
        getline(inFile, line);
        if (inFile.eof())
        {
            cerr<<"End of file"<<endl;
            cerr<<"I could read only "<<i-1<<" points"<<endl;
            break;
        }
        stringstream input(line);
        float *point = new float[dp_length];
        for (long d=0;d<dp_length;d++)
        {
            input>>point[d];
        }
        // Convert this float array into a MVP point instance.
        pointlist[i] = get_data_point(point, dp_length, point_id);
        point_id++;
        delete[] point;
    }
    inFile.close();
    return pointlist;
}

int main(int argc, char ** argv)
{
    options opt;
    //Initialize options to Undef values
    opt.datafile = UNDEF_STR;
    opt.bqfile = UNDEF_STR;
    opt.rqfile = UNDEF_STR;
    opt.range = 0.;
    opt.dim = UNDEF_LONG;
    opt.skip = 0;
    opt.count = 99999999;
    opt.help = false;
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
    if (opt.skip > 0)
    {
        cerr<<"Skip option is not yet supported. Will use 0"<<endl;
    }

    // We are going to use L2 distance measure
    // Hence set the distance function accordingly
    /*
        sisobus modi
       */
//    CmpFunc distance_func = point_l2_distance;
    CmpFunc distance_func = point_hamming_distance;
    nbcalcs = 0;

    float constructionTime = 0.;
    clock_t start = clock();
    // Read all the points.
    MVPDP ** pointlist = read_all_points(opt.count, (unsigned int)opt.dim, opt.datafile.c_str());
    
    /* create tree */
    MVPTree *tree = mvptree_alloc(NULL, distance_func, MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP);
    assert(tree);

    printf("Add list of %d points to tree.\n", opt.count);
    /* add points to tree */
    MVPError err = mvptree_add(tree, pointlist, opt.count);
    clock_t end = clock();
    printf("[END]Add list of %d points to tree.\n", opt.count);
    constructionTime = end-start;
    constructionTime /= CLOCKS_PER_SEC;
    unsigned long long constructionDistCalc = nbcalcs;
    if (err != MVP_SUCCESS)
    {
        cerr<<"Error "<<endl<<mvp_errstr(err)<<endl;
    }
    //assert(err == MVP_SUCCESS);


    //FILE *outfile = fopen("mvptree_out.txt", "w");
    //mvptree_print(outfile, tree);
    //fclose(outfile);
    
    /*printf("Write tree out to file, %s.\n", testfile);
    //write out the tree
    err = mvptree_write(tree, testfile, 00755);
    assert(err == MVP_SUCCESS);
    
    mvptree_clear(tree, free);
    free(tree);


    printf("Read tree from %s.\n", testfile);
    // read it back in again
    tree = mvptree_read(testfile, distance_func, MVP_BRANCHFACTOR, MVP_PATHLENGTH, MVP_LEAFCAP, &err);
    assert(tree);


    //add some points to a tree one at a time
    MVPDP *querypoint = NULL;
    int count = 0, total = nbpoints/10;
    do {
    MVPDP *pnt = generate_point(dplength);
    if (count == 0) querypoint = mvpdp_dup(pnt);

    printf("Add point, %s to tree.\n", pnt->id);
    err = mvptree_add(tree, &pnt, 1);
    assert(err == MVP_SUCCESS);
    } while (++count < total);


    printf("Retrieve point, %s.\n", cluster_center->id);

     // retrieve the cluster center
    nbcalcs = 0;
    unsigned int nbresults;
    MVPDP **results = mvptree_retrieve(tree, cluster_center , knearest, radius, &nbresults, &err);
    assert(results);
    assert(err == MVP_SUCCESS);

    unsigned int i;
    for (i = 0;i < nbresults;i++){
    fprintf(stdout,"  FOUND --> (%d) %s\n", i, results[i]->id);
    }
    free(results);
    */
    if (opt.rqfile.compare(UNDEF_STR) && opt.range > 0.)
    {
        // RQ file is defined and the radius is defined
        ifstream rqFile;
        string line;
        rqFile.open(opt.rqfile.c_str(), ios::in);
        getline(rqFile, line);

        float avgHits = 0.;
        float avgDistCalculations = 0.;
        long queryCnt = 0;
        // For a range query the upper bound on number of
        // neighbors is equal to the number of points.
        // Apparently the MVP tree code has a minor bug
        // if the number of points retrieved is exactly K
        // then it throws an error. Hence we set k to be 1 more
        // than actual max value.
        unsigned int knearest = (unsigned int)opt.count + 1;
        clock_t totalTime = 0;
        float *buffer = new float[opt.dim];
        int numberOfRQ = 0;
        start = clock();
        while (!rqFile.eof())
        {
//            printf("%d\n",numberOfRQ++);
            // This not the first non-comment line so it contains vectors
            if(line.length() > 3) // vector must have at least 3 characters
            {
                stringstream input(line);
                for(long d =0;d<opt.dim;d++)
                {
                    input>>buffer[d];
                }
                unsigned int nbresults;
                // Number of distance calculations
                nbcalcs = 0;
                // Convert float array to MVP point.
                MVPDP * querypoint = get_data_point(buffer, (unsigned int) opt.dim, 0);
                /*    float * data = (float*)querypoint->data;
                    for (long dd=0;dd<opt.dim;dd++)
                    {
                        cout<<data[dd]<<" , ";
                    }
                    cout<<endl;*/
                MVPDP** results = mvptree_retrieve(tree, querypoint, knearest, opt.range, &nbresults, &err);
                assert(results);
                if(err != MVP_SUCCESS)
                {
                    cerr<<"Error "<<endl<<mvp_errstr(err)<<endl;
                    return 1;
                }

                /*for (unsigned int i=0;i<nbresults;i++)
                {
                    float * data = (float*)results[i]->data;
                    fprintf(stdout,"%d  \n", results[i]->id);
                }*/
                avgHits += nbresults;
                avgDistCalculations += nbcalcs;
                queryCnt ++;
              //  cout<<"ACC,9999"<<",HITS,"<<nbresults<<",DIST,"<<nbcalcs<<endl;
                free(results);
                dp_free(querypoint, free);
            }
            getline(rqFile, line);
        }
        end = clock();
        totalTime += end -start;

        delete[] buffer;
        rqFile.close();

        avgHits /= queryCnt;
        avgDistCalculations /= queryCnt;
        float avgTime = (float)totalTime/(CLOCKS_PER_SEC*queryCnt);
        cout<<"CONST_TIME,"<<constructionTime<<",CONST_DIST,"<<constructionDistCalc<<",AVG_HITS,"<<avgHits<<",AVG_DIST,"<<avgDistCalculations<<",AVG_TIME,"<<avgTime<<endl;
    }
    mvptree_clear(tree, free);
    free(tree);

 cleanup:
//    dp_free(cluster_center, free);
    //printf("Done.\n");
    return 0;

}
