#ifndef USERINTERFACE_H
#define USERINTERFACE_H

#include "KDTree.h"
#include <ctime>
#include <cstdlib>
#include <map>
#include <fstream>
#include <unistd.h>
#include <cmath>
#include <utility>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/variance.hpp>

using namespace std;
using namespace boost::accumulators;

#define OPT_DIM "-dim"
#define OPT_ORIG_DIM "-orig_dim"
#define OPT_LOAD_FILE "-load_file"
#define OPT_BQFILE "-bqfile"
#define OPT_RQFILE "-rqfile"
#define OPT_KNNFILE "-knnfile"
#define OPT_REF_FILE "-ref_file"
#define OPT_RANGE "-range"
#define OPT_BEGIN "-skip"
#define OPT_COUNT "-count"
#define OPT_KNN "-knn"
#define OPT_STATS "-stats"
#define OPT_HELP "-help"
#define OPT_REFMETHOD "-ref_method"

#define UNDEF_STR "__UNDEF__"
#define UNDEF_LONG -999999

// Various methods for reference point generation
#define VELTKAMP "veltkamp"
#define ANGLE_DEV "angle_dev"
#define RANDOM "random"
#define RANDOMCORNER "randomcorner"
#define MVP "mvp"

#define TEST_VERSION "test_version"

// Various control parameters
///////     VELTKAMP METHOD    //////
#define SPACING_THR 0.0001
#define CORR_THR 0.95
/////// DATA INDEPENDENT METHOD//////
#define ITERATIONS 1000
#define BLOWUP 1
struct options
{
    // Name of the data file.
    string datafile;
    string bqfile;// box query file
    string rqfile;// range query file
    string reffile;// Reference point  file
    string knnfile; // K-NN query file
    string refMethod; // Reference point generation method
    float range;
    long knn;
    long dim;
    long origDim;
    long skip;
    long count;
    bool stats;
    bool help;
};

/**
 * Options parser
 * @param argCount Number of arguments
 * @param argVector Argument vector
 * @param opt object that will contain parsed options
 */
void getOptions(int argCount, char **argVector, options * opt);

// Displays the help message.
void displayHelp();

// Generate reference points using Veltkamp's method
void generateVeltkampRefPoints(float **refPoints, long numRef, float ** dataPoints, long numData, long dim);

// Generate random data points as reference points
void generateRandomRefPoints(float **refPoints, long numRef, float ** dataPoints, long numData, long dim);

// Generate random corner points as reference points
void generateRandomCornerRefPoints(float **refPoints, long numRef, long dim);

// Generate reference points using data independent method
// A slower more exhaustive search
void generateDataIndependentRefPoints2(float **refPoints, long numRef, long dim);

// Generate reference points using data independent method
void generateDataIndependentRefPointsRandomized(float **refPoints, long numRef, long dim);

// Generate reference points using data independent method - tuned to handle very high
// dimension
void generateDataIndependentRefPointsHighDim(float **refPoints, long numRef, long dim);

// Generate reference points using data independent method
void generateDataIndependentRefPoints(float **refPoints, long numRef, long dim);

// A recursive function to generate strings of 0,1
// used as an intermediate function for reference point generator.
vector <string> recursiveaGenerator(long dim);

// Generate reference points using recursive data independent method
void generateDataIndependentRefPointsRecursive(float **refPoints, long numRef, long dim);

// Helper function - calculates std deviation of spacing
float spacingVariance(float* distances, long numDist);

// Calculates correlation of two arrays and returns false if
// it is greater than the threshold. Returns true otherwise.
bool correlationNotOk(float * distVector1, float* distVector2, long length);

// Another helper function. First calculates all possible
// angles among the reference points. Then calculates
// the variance of those angles.
float getAngleVariance(float **setOfPoints, long numRef, long dim);

// Returns variance of pairwise distances
float getDistVariance(float **refPoints, long numRef, long dim);

// Returns true if the set of point contains a repeated point.
bool containsRepeatedPoint(float **refPoints, long numRef, long dim);

#endif
