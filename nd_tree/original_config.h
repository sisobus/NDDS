#ifndef CONFIG_H
#define CONFIG_H

#include <stack>

// This file defines all the configuration parameters for the ND Tree program

//#define LOG_SPLITTING_VECTOR_INDEX
//#ifdef LOG_SPLITTING_VECTOR_INDEX
//	long int vector_index;
//#endif

//#define LOG_SPLITTING_HEIGHT
//#ifdef LOG_SPLITTING_VECTOR_INDEX
//	int splitting_height;
//#endif

int debug_dataIndex;

enum ENVIRONMENT {UNIX, WINDOWS};
ENVIRONMENT RUNNING_ENVIRONMENT = UNIX;

const int QUERY_RESULTS_BUFFER_SIZE = 1000;

typedef int ND_tree_record;
const int MAX_LINE_IN_SOURCE_FILE = 100000; // MUST_CHANGE_FOR_TEST

const int DIM = 10; // MUST_CHANGE_FOR_TEST
const int CNTDIM = 0; 

//rest 2 are used as part of sourceData file name: sourceData32+10
const int TOTAL_DSC_VALUE = DIM;
const int TOTAL_CNT_VALUE =0;//total CNT values to be read out
const int TOTAL_BOX_QUERY_NUM = 100;
const int TOTAL_RANGE_QUERY_NUM =100;


const int ALPHA=10;
int A[] = {ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA
, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA,ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA, ALPHA}; 
const int BYTES_PER_DIM_IN_DMBR = ((ALPHA%8)==0)?(ALPHA/8):(ALPHA/8+1);

const int BOX_SIZE_START_AT = 1;
const int BOX_SIZE_STEP = 1;
const int BOX_SIZE_STOP_BEFORE = A[0];
const int RANGE_SIZE_STEP = 1;
const int RANGE_STOP_BEFORE = (DIM + CNTDIM)<10?(DIM + CNTDIM):10;


const int DISK_BLOCK_SIZE = (1<<14);
//const int DISK_BLOCK_SIZE = 8192;
//const int DISK_BLOCK_SIZE = 1024; 
const int DMBR_SIZE = DIM * BYTES_PER_DIM_IN_DMBR; 



const int TOTAL_DIM =DIM + CNTDIM; //<--only used for leaf node LEAF_NODE_SIZE calculation

int boxSize;

const int usedBytesForEntryHeader = ((DIM%8)==0)?(DIM/8):(DIM/8+1);;
//const int USE_LINK_FOR_BOX_QUERY=1;//don't change, set to 1 in all conditions
//const int USE_LINK_FOR_RANGE_QUERY=1;//don't change, set to 1 in all conditions
		
const int MAX_DIM_AND_CONTDIM = DIM; //jan 09, maximum lines for cnt and dsc in boxqueryall file

enum SPLIT_TYPE { ORIGINAL, TO_DEATH/*do not use*/,TO_DEATH_MAX_UTIL/*do not use*/,TO_DEATH_RANDOM_UTIL }; 
SPLIT_TYPE nodeSplitType  = ORIGINAL; // NEVER CHANGE, IS USED TO CALL FUNCTION split_algorithm_addPerimeter
//SPLIT_TYPE nodeSplitType  = TO_DEATH_RANDOM_UTIL; // NEVER CHANGE, IS USED TO CALL FUNCTION split_algorithm_addPerimeter

enum ENUM_TREE_TYPE {STATIC_TREE, DYNAMIC_TREE};ENUM_TREE_TYPE TREE_TYPE = STATIC_TREE;

bool RESORT_2_OLD_SPLIT = true;//if true, might resort to old split when nodeSplitType is TO_DEATH_RANDOM_UTIL

//#define LOG_MBR_SPLITTING_LEAF
//#define LOG_MBR_SPLITTING_DIR

//#define LOG_EDGE_LENGTH_LEAF
//#define LOG_EDGE_LENGTH_DIR

//#define LOG_OLD_AND_NEW_BLOCK_LEAF

////#define LOG_UTILIZATION

//#define LOG_ENTRY_ON_SPLITTING_DIM 
//#define ShorterThanSplitDim

//#define LOG_KANPSACK_VALUE_AND_WEIGHT

		int heuristics_overlap_used_leaf =0;
		int heuristics_area_used_leaf =0;
		int heuristics_overlap_used_dir =0;
		int heuristics_area_used_dir =0;
		int BestChild_covered_area =0;
		int BestChild_notcovered_overlap_enlarge =0;
		int BestChild_notcovered_area =0;
		int BestChild_notcovered_area_enlarge =0;

const int enforce_minUtil_for_exhaustSplit =1;//always 1
const int try_all_dim_with_minUtil = 1;



#endif
