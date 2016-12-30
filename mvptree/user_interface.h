using namespace std;

#define OPT_DIM "-dim"
#define OPT_LOAD_FILE "-load_file"
#define OPT_BQFILE "-bqfile"
#define OPT_RQFILE "-rqfile"
#define OPT_RANGE "-range"
#define OPT_BEGIN "-skip"
#define OPT_COUNT "-count"
#define OPT_HELP "-help"
#define UNDEF_STR "__UNDEF__"
#define UNDEF_LONG -999999

// I know this progam is in C. But using
// some C++ components lets me copy-paste
// a lot of existing code.
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


struct options
{
    // Name of the data file.
    string datafile;
    string bqfile;// box query file
    string rqfile;// range query file
    string knnfile; // K-NN query file
    float range;
    long dim;
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

/**
 * Convert a float array in to a MVP data point.
 * This function assumes that sizeof(float) = 4
 * dp_length is the number of dimensions
 */
MVPDP* get_data_point(float * point, unsigned int dp_length, long point_id);

/** Sort of a copy constructor
 * Returns a copy of the point passed as the argument.
 */
MVPDP* mvpdp_dup(MVPDP *dp);

// Function to calculate L2 distance
float point_l2_distance(MVPDP *pointA, MVPDP *pointB);

