
#include "UserInterface.h"
#include <unistd.h>

void displayHelp()
{
    cout<<"Sorry. Not yet implemented :-)\nYou are on your own."<<endl;
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
//    string optKnnfile(OPT_KNNFILE);
//    string optKnn(OPT_K);
    string optRange(OPT_RANGE);
    string optSkip(OPT_BEGIN);
    string optCount(OPT_COUNT);
    string optStats(OPT_STATS);
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
        if(!optReffile.compare(argVector[i]))
        {
            opt->reffile = argVector[i+1];
            i++;
        }
        /*if(!optKnnfile.compare(argVector[i]))
        {
            opt->knnfile = argVector[i+1];
            i++;
        }
        if(!optKnn.compare(argVector[i]))
        {
            stringstream s(argVector[i+1]);
            s>>(opt->knn);
        }*/
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
    options opt;
    //Initialize options to Undef values
    opt.datafile = UNDEF_STR;
    opt.bqfile = UNDEF_STR;
    opt.rqfile = UNDEF_STR;
    opt.reffile = UNDEF_STR;
    //opt.knnfile = UNDEF_STR;
    //opt.knn = UNDEF_LONG;
    opt.range = 0.;
    opt.origDim = UNDEF_LONG;
    opt.dim = UNDEF_LONG;
    opt.skip = 0;
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

    // Create a KD Tree

    if (opt.count > 100000)
    {
        cout<<"The KD-Tree contains more than 100000 points. It may take several minutes. Be patient."<<endl;
    }

    // Create an empty KD tree with given dimensions.
    float constructionTime = 0.;
    clock_t start = clock();
    KDTree kdTree(opt.dim, opt.origDim);
    clock_t end = clock();
    constructionTime += end-start;

    ifstream inFile;
    inFile.open(opt.datafile.c_str(), ios::in);

    // Read reference points
    float ** refPoints = new (nothrow) float *[opt.dim];
    for (long rf=0; rf < opt.dim;rf++)
    {
        refPoints[rf] = new float[opt.origDim];
    }
    ifstream refFile;
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
    if (opt.count > 100000)
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
    cout<<" Space allocated"<<endl;

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
    cout<<"File reading complete"<<endl;

    if (opt.reffile.compare(UNDEF_STR))
    {
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
        cout<<" Refrence points read"<<endl;
    }
    else
    {
        // User did not pass a value for ref file
        // Use MVP-Tree (TODS version) heuristic
        // Select a random point.
        srand ( time(NULL) );
        for (long rf=0; rf < opt.dim;rf++)
        {
            int randomIndex = rand() % opt.count;
            for (long d=0;d<opt.origDim;d++)
            {
                refPoints[rf][d] = allPoints[randomIndex][d];
            }
        }
        cout<<" Refrence points generated"<<endl;
    }
    // Reference point reading complete
        /*for (long rf=0; rf < opt.dim;rf++)
        {
            for(long d=0;d<opt.origDim;d++)
            {
                cout<<refPoints[rf][d]<<",";
            }
            cout<<endl;
        }*/
    // Insertion begins
    start = clock();
    for(int i=0;i< opt.count ; i++)
    {
        // Calculate distances from reference points
        for (long rf=0;rf < opt.dim;rf++)
        {
            distanceBuffer[rf] = getL2Distance(refPoints[rf], allPoints[i], opt.origDim);
        }
        // Create a KD-tree node.
        kdTree.insert(distanceBuffer, allPoints[i], countInserted);
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
    cout<<" Height of the tree is "<<kdTree.getHeight()<<endl;
    constructionTime /= CLOCKS_PER_SEC;
//    kdTree.printInOrder();

    // Range query
    // We do not execute the range query directly.
    // We first convert it to a distance based
    // box query.
    float *buffer= new float[opt.origDim];
    if (opt.rqfile.compare(UNDEF_STR) && opt.range > 0.)
    {
        // RQ file is defined and the radius is defined
        ifstream rqFile;
        rqFile.open(opt.rqfile.c_str(), ios::in);

        float avgHits = 0.;
        float avgComparisons = 0.;
        float avgDistCalculations = 0.;
        long queryCnt = 0;
        clock_t totalTime = 0;

        getline(rqFile, line);
        start = clock();
        while (!rqFile.eof())
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
                // Convert range query to box query
                QueryBox qbox;
                for (long rf=0;rf < opt.dim;rf++)
                {
                    qbox.cntBox[2*rf] = distanceBuffer[rf] - opt.range;
                    qbox.cntBox[2*rf+1] = distanceBuffer[rf] + opt.range;
                }
                // Execute "box query". We pass the original
                // query parameters to be able to eliminate FP.
                ResultSet rs = kdTree.rangeQuery(&qbox, buffer, opt.range);
                avgHits += rs.countHits;
                avgComparisons += rs.countComparisons;
                avgDistCalculations += rs.countDistComputations;
                queryCnt ++;

                cout<<"HITS,"<<rs.countHits<<",COMP,"<<rs.countComparisons<<",DIST,"<<rs.countDistComputations<<endl;
            }
            getline(rqFile, line);
        }
        end = clock();
        totalTime += (end -start);

        rqFile.close();

        avgHits /= queryCnt;
        avgComparisons /= queryCnt;
        avgDistCalculations /= queryCnt;
        float avgTime = (float) totalTime / (CLOCKS_PER_SEC*queryCnt);
        cout<<"CONST_TIME,"<<constructionTime<<",AVG_HITS,"<<avgHits<<",AVG_COMP,"<<avgComparisons<<",AVG_DIST,"<<avgDistCalculations<<",AVG_TIME,"<<avgTime<<",NumQ,"<<queryCnt<<endl;
    }
    //cerr<<" sleeping so that memory snapshot is captured"<<endl;
    //sleep(60);
    delete[] buffer;
    delete[] distanceBuffer;
}
