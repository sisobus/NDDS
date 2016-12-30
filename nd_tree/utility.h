#ifndef UTILITY_H
#define UTILITY_H

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <time.h>
#include <limits>
#include <cassert>
#include <math.h>
#include <vector>
#include <utility>//for std::pair
#include <algorithm>
#include <stack>

using namespace std;

enum Error_code { duplicate_error, not_present, success, overflow, underflow, range_err, fail };

const int INT_INF = numeric_limits<int> :: max( );
const double DOUBLE_INF = numeric_limits<double> :: max( );
const double DOUBLE_THRESHOLD = 1E-14;

const int BITS_PER_BYTE = 8;

const unsigned char MASKS[] = {1, 2, 4, 8, 16, 32, 64, 128};

// A lookup table counts the number of bit 1 in each number between 0 and 255
const int bit_1_count_lut[] = { // 256 entries
   0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
   1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
   2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
   3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
   4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
};

const unsigned char bit_mask[] = {1,2,4,8,16,32,64,128};//set corresponding bit to 1, others to 0



/* Code to compute the above table
   int tmp_i = 0;
   for(unsigned char ch = 0; tmp_i < 256; ch++, tmp_i++){
      int tmp_sum = 0;
      for(int bit_index = 0; bit_index < 8; bit_index++)
         if(MASKS[bit_index] & ch)tmp_sum++;
      cout<< tmp_sum << ", ";
   }
   cout << endl;
*/

int string_to_int(string s);
string int_to_string(int i);

// pad_s is single char, padded to the left
string lpad(string s, unsigned int str_len, string pad_s);

// sort weight_list, sorted_index_list is also an output, ascending by default
void NDT_sort(int count, double* weight_list, int* sorted_index_list, bool asc = true);
void NDT_stable_sort(int count, double* weight_list, int* sorted_index_list, bool asc = true);

void NDT_sort(int count, unsigned int* weight_list);
void NDT_sort(unsigned int* a, int lo, int hi);
void binary2string(const unsigned char* const ptr, int lengthInBytes, string & resultString);
void string2binary(const string &origString, unsigned char* const ptr);

vector<int> Knapsack_recursive_forNDTree(
	int n, 
	int available_size, //initialy equals knapsack_capacity = block size - overhead
const int * const item_size,
const int * const item_value,
const int 	all_entries_size,
const float knapsack_min_capacity,
const int knapsack_capacity	);

#endif
