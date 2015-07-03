#ifndef DIR_NODE_H
#define DIR_NODE_H

#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include "logClass.h"
#include "Dir_entry.h"

const int DIR_NODE_OVERHEAD = 4; // int count
extern logClass logO;
extern int debug_height;


class Dir_node:public Node{
public:
   Dir_node(int* alphabet_sizes);
   int count; // curent number of entries in the node

   //int get_count();
   //void find_covering_entries(const Leaf_entry &query_data, int &num_of_covering_entries, unsigned int* covering_block_nos); // Given query data, find all the entries in node that covers the query data
   void find_covering_entries(const Leaf_entry &query_data, vector<int> &covered_block_numbers); // Given query data, find all the entries in node that covers the query data


   void find_entries_by_hamming_dist(const Leaf_entry &query_data, int range, vector<int>& block_numbers_in_range);
   void find_best_child(const Leaf_entry &query_data, int* alphabet_sizes, unsigned int& child_block_no, int& child_index);
   void set_DMBR(int entry_index, unsigned char* new_DMBR, bool& not_same); // set the DMBR of a particular entry given by entry_index, if not the same, cal and return a new DMBR
   Error_code insert_new_entry(unsigned int new_entry_block_no, unsigned char* new_entry_DMBR, double dir_min_util, int *alphabet_sizes, unsigned char* cur_DMBR, Dir_node*& new_dir_node, unsigned char* new_DMBR);

   //void find_entries_by_box_query(	const Dir_entry box_query_data, int& num_of_relevant_entries, unsigned int* relevant_block_nos);
   void find_entries_by_box_query(	const Dir_entry box_query_data, vector<int>& block_numbers_in_box);

   void set_node_count(int new_count);
   void set_node_entry(int index, Dir_entry& new_entry);
   void set_node_entry(int index, unsigned int new_entry_block_no, unsigned char* new_entry_DMBR);

   int get_node_count();
   void get_node_entry(int index, Dir_entry& cur_entry);
   void get_DMBR(unsigned char* cur_DMBR);
   void get_DMBR(unsigned char* cur_DMBR,int index);

   //virtual void read_node_dynamic(fstream& ND_file, unsigned int block_number);
   //virtual void write_node_dynamic(fstream& ND_file, unsigned int block_number);
   void read_node_dynamic(fstream& ND_file, unsigned int block_number);
   void write_node_dynamic(fstream& ND_file, unsigned int block_number);
   void read_node_static(fstream& ND_file, unsigned int block_number);
   void write_node_static(fstream& ND_file, unsigned int block_number);


   // public data members
   static const int DIR_NODE_SIZE =  (DISK_BLOCK_SIZE - DIR_NODE_OVERHEAD) / (DMBR_SIZE + sizeof(unsigned int));


void log2file( );
void logDMBR( );
void print_OneEntry_OnOneDscDim( Dir_entry*  oneEntry,int entryIndex,int dimNum,int* alphabet_sizes );
void log_OneEntry_OnOneDscDim( Dir_entry*  oneEntry,int entryIndex,int dimNum,int* alphabet_sizes );
void logSingleChildDMBR(const unsigned char* const DMBR);

private:
   //static const int max_dynamic_DIR_NODE_SIZE =  (DISK_BLOCK_SIZE - DIR_NODE_OVERHEAD) / (usedBytesForEntryHeader + sizeof(unsigned int));

   // private types
   // used in sort entries by size
   struct set_entry{
      unsigned char set[MAX_DMBR_DIM_SIZE]; // bitmap
      int element_size; // num of elements in the set
      int freq; // number of entries whose set corresponds to the set_entry
      //int entry_set[DIR_NODE_SIZE + 1]; // entries (indices) whose component set corresponds to the set_entry
      int *entry_set; // entries (indices) whose component set corresponds to the set_entry
	  int num_of_entry_set;
	  set_entry(){
		  num_of_entry_set=0;
		  entry_set = NULL;
	  };
	  ~set_entry(){delete []entry_set ;}

	  set_entry(const set_entry& orig) {

			  num_of_entry_set=orig.num_of_entry_set;
			  entry_set=new int[num_of_entry_set];

		  for(int x=0;x<num_of_entry_set;x++) entry_set[x]=orig.entry_set[x];

		  for(int x=0;x<MAX_DMBR_DIM_SIZE;x++) set[x]=orig.set[x];

		  element_size = orig.element_size;
		  freq = orig.freq;

	  }

	  set_entry& operator=(const set_entry& orig) {
		  //assert ((count+1)==sizeof(orig.entry_set));
		  if(this !=&orig)
		  {
			  if(num_of_entry_set<orig.num_of_entry_set)
			  {
				  delete []entry_set;
				  num_of_entry_set=orig.num_of_entry_set;
				  entry_set=new int[num_of_entry_set];
			  }
			  else
				  num_of_entry_set=orig.num_of_entry_set;


			  for(int x=0;x<num_of_entry_set;x++) entry_set[x]=orig.entry_set[x];

		  }

		  for(int x=0;x<MAX_DMBR_DIM_SIZE;x++) set[x]=orig.set[x];

		  element_size = orig.element_size;
		  freq = orig.freq;

		  return *this;

	  }


	  void resizeSetEntry(int numSet)
	  {
			  if(num_of_entry_set<numSet)
			  {
				  delete []entry_set;
				  num_of_entry_set=numSet;
				  entry_set=new int[num_of_entry_set];
			  }
			  else
				  num_of_entry_set=numSet;
	  }
   };

   struct aux_tree_node{
      unsigned char letters[MAX_DMBR_DIM_SIZE]; // bitmap representation
      //int sets[DIR_NODE_SIZE + 1]; // array of indices into set_array
      int *sets; // array of indices into set_array
      int set_count;
      int freq; // number of entries whose set corresponds to the set in sets
      int level; // level in the tree.  Leaf is level 1
      int children[MAX_ALPHABET_SIZE];
      int cnt; // number of children
      int parent; // index of parent
	  int numLetters;	
	  int num_of_sets;

	  aux_tree_node(){
		  num_of_sets=0;
		  sets=NULL;
	  };
	  ~aux_tree_node(){delete []sets ;}

	  aux_tree_node(const aux_tree_node& orig) {
		  num_of_sets=orig.num_of_sets;
		  sets=new int[num_of_sets];
		  for(int x=0;x<num_of_sets;x++) sets[x]=orig.sets[x];

		  for(int x=0;x<MAX_DMBR_DIM_SIZE;x++) 
		  {
			  letters[x]=orig.letters[x];
			  children[x]=orig.children[x];
		  }

		  set_count = orig.set_count;
		  freq = orig.freq;
		  level = orig.level;
		  cnt = orig.cnt;
		  parent = orig.parent;
	  }

	  aux_tree_node& operator=(const aux_tree_node& orig) {
		  if(this !=&orig)
		  {

			  if(num_of_sets<orig.num_of_sets)
			  {
				  delete []sets;
				  num_of_sets=orig.num_of_sets;
				  sets=new int[num_of_sets];
			  }
			  else
				  num_of_sets=orig.num_of_sets;

			  for(int x=0;x<num_of_sets;x++) sets[x]=orig.sets[x];

			  for(int x=0;x<MAX_DMBR_DIM_SIZE;x++) 
			  {
				  letters[x]=orig.letters[x];
				  children[x]=orig.children[x];
			  }

			  set_count = orig.set_count;
			  freq = orig.freq;
			  level = orig.level;
			  cnt = orig.cnt;
			  parent = orig.parent;
		  }
		  return *this;

	  }

	  void resizeSets(int numSet)
	  {
			  if(num_of_sets<numSet)
			  {
				  delete []sets;
				  num_of_sets=numSet;
				  sets=new int[num_of_sets];
			  }
			  else
				  num_of_sets=numSet;
	  }
   };

   struct mixed_list_entry{
      int index; // If set, index of a set_entry.  If roots, index in the aux_tree
      bool is_set; // true: sets, false: roots
      int weight; // for roots only
   };

   // private methods
   bool is_covered(const Leaf_entry &query_data, unsigned char* DMBR); // check whether the DMBR covers the query data
   bool is_within_Hamming_dist(const Leaf_entry &query_data, unsigned char* DMBR, int range);
   void enlarge_DMBR(unsigned char* DMBR, const Leaf_entry &new_data); // Enlarge the DMBR by union the new data
   void cal_DMBR(Dir_entry* dir_entries, int num_of_dir_entries, unsigned char* new_DMBR); // Create a DMBR from a set of dir_entries
   void cal_DMBR(Dir_entry* dir_entries, int* indices_of_dir_entries_used, int num_of_dir_entries_used, unsigned char* new_DMBR); // Create a DMBR from a subset of dir_entries
   void cal_DMBR(vector<Dir_entry> & dir_entries, int num_of_dir_entries, unsigned char* new_DMBR); // Create a DMBR from a set of dir_entries
   void cal_DMBR(vector<Dir_entry> & dir_entries, int* indices_of_dir_entries_used, int num_of_dir_entries_used, unsigned char* new_DMBR); // Create a DMBR from a subset of dir_entries

   void split_algorithm(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
   void split_algorithm_after_knapsack(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);

   void split_algorithm_TO_DEATH(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
   void split_TO_DEATH_RANDOM_UTIL(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
   void split_TO_DEATH_RANDOM_UTIL_v2(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
   void split_TO_DEATH_RANDOM_UTIL_v3(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);
   void split_TO_DEATH_RANDOM_UTIL_v4(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2);

   void sort_entries_by_size(int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list);
   void sort_entries_by_size_Random_Util_v2(int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list, int & rst_leftNumLetters, double & rst_overlap, int &rst_cutPoint);
   void sort_entries_by_size_Random_Util_v3(double dir_min_util, int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list, int & rst_leftNumLetters, double & rst_overlap, int &rst_cutPoint);


   int sort_entries_by_size_RANDOM_UTIL(int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list,bool forceSplitonThisDim=false);

   void decide_order_by_aux_tree(aux_tree_node* tree, int root, int letters_bitmap_byte_size, set_entry* set_entry_array, int set_entry_array_size, int* sorted_set_list, int& list_size);
   void find_minOverlap_by_aux_tree(int sort_dim, Dir_entry* sort_entries, aux_tree_node* tree, int root, int letters_bitmap_byte_size, set_entry* set_entry_array, int set_entry_array_size, vector<int> & leftSets,  double & rst_overlap);

	bool is_within_box(	const Dir_entry box_query_data, 	unsigned char* DMBR);

	int getCompressedEntriesSize(vector<Dir_entry> & tmp_entries, int *alphabet_sizes);

	int getCompressedDiskSize(vector<Dir_entry> & tmp_entries, int *alphabet_sizes);
   // Data members

   //Dir_entry entries[DIR_NODE_SIZE]; // Array of entries
	vector<Dir_entry> entriesVec;

};



Dir_node::Dir_node(int* alphabet_sizes):Node(alphabet_sizes){
   count = 0;
   entriesVec.clear();
}



//int Dir_node::get_count(){
//	assert (count==entriesVec.size());
//   return count;
//}



void Dir_node::find_covering_entries(const Leaf_entry &query_data, vector<int> &covered_block_numbers){
   int i;

   //num_of_covering_entries = 0;

   for(i = 0; i < count; i++){
      if(is_covered(query_data, entriesVec.at(i).DMBR)){ // Current entry covers the data
         //covering_block_nos[num_of_covering_entries] = entriesVec.at(i).child;
         //num_of_covering_entries++;
			covered_block_numbers.push_back(entriesVec.at(i).child);
      }
   }
} 


void Dir_node::find_entries_by_box_query(
	const Dir_entry  box_query_data,
	vector<int>& block_numbers_in_box)
{
	int i;


	for(i = 0; i < count; i++)
	{

		if(is_within_box(box_query_data, entriesVec.at(i).DMBR))
		{ 
			block_numbers_in_box.push_back(entriesVec.at(i).child);

		}

	}

}


void Dir_node::find_entries_by_hamming_dist(const Leaf_entry &query_data, int range, vector<int>& block_numbers_in_range){
   int i;



   for(i = 0; i < count; i++){
      if(is_within_Hamming_dist(query_data, entriesVec.at(i).DMBR, range)){ 
         block_numbers_in_range.push_back(entriesVec.at(i).child);
      }
   }
}



void Dir_node::find_best_child(const Leaf_entry &new_data, int* alphabet_sizes, unsigned int& child_block_no, int& child_index){
   int i, j;
	int lastHeuristicUsed;
   unsigned int best_block_no;
   int best_child_index;
   bool proper_child_found;
   double best_area, cur_area;

   proper_child_found = false;
   best_area = DOUBLE_INF;
   for(i = 0; i < count; i++){
      if(is_covered(new_data, entriesVec.at(i).DMBR)){
         proper_child_found = true;
         cur_area = cal_normed_area(entriesVec.at(i).DMBR, alphabet_sizes);
         if(cur_area < best_area){
            best_area = cur_area;
            best_block_no = entriesVec.at(i).child;
            best_child_index = i;
			lastHeuristicUsed = 1;
         }
      }
   }
 
#ifdef LOG_WHICH_HEURISTICS_CHOOSE_LEAF
   if(proper_child_found && (debug_height==2))
	   logO.log2File(" is covered ");
#endif

   if(!proper_child_found)
   {
      double orig_area, new_area;
      double best_overlap_enlarge, cur_overlap_enlarge, old_overlap, new_overlap;
      double best_enlarge, cur_enlarge;

      // sort children based on area enlargement
      //double enlargement[DIR_NODE_SIZE]; // store the enlargement of each child
      double *enlargement=new double[count]; // store the enlargement of each child
      
      // Used to indicate a group of DMBRs
      unsigned char tmp_DMBR[DMBR_SIZE];

      for(i = 0; i < count; i++){ // check the enlargement of all children
         for(j = 0; j < DMBR_SIZE; j++) // copy the DMBR
            tmp_DMBR[j] = entriesVec.at(i).DMBR[j];
         enlarge_DMBR(tmp_DMBR, new_data); //enlarge the tmp DMBR using new data

         orig_area = cal_normed_area(entriesVec.at(i).DMBR, alphabet_sizes);
         new_area = cal_normed_area(tmp_DMBR, alphabet_sizes);
         enlargement[i] = new_area - orig_area;
      }
      //int sorted_child_list[DIR_NODE_SIZE]; // store the sorted child list
      int *sorted_child_list=new int[count]; // store the sorted child list


      // Sort the child based on enlargement
      NDT_sort(count, enlargement, sorted_child_list);
      
      best_overlap_enlarge = DOUBLE_INF;
      best_enlarge = DOUBLE_INF; 
      best_area = DOUBLE_INF;
		best_child_index=-1;

      int cur_child;
      for(i = 0; i < count; i++){ // check all child in sorted order
         cur_child = sorted_child_list[i];
         for(j = 0; j < DMBR_SIZE; j++) // copy the DMBR
            tmp_DMBR[j] = entriesVec.at(cur_child).DMBR[j];
         enlarge_DMBR(tmp_DMBR, new_data); //enlarge the tmp DMBR using new data
 
		 if(nodeSplitType!=ORIGINAL)
		 {
			 bool enlargeDIMwithOneLetter=false;
			 for(int k = 0; k < DIM; k++)
			 {
				 int tmp_sum1 = 0;
				 int tmp_sum2 = 0;
				 for(j = this->DMBR_start_byte_lut[k]; j <= this->DMBR_end_byte_lut[k]; j++)
				 {
					 tmp_sum1 += bit_1_count_lut[tmp_DMBR[j]];
					 tmp_sum2 += bit_1_count_lut[entriesVec.at(cur_child).DMBR[j]];
				 }

				 if((tmp_sum1!=1)&&(tmp_sum2==1))
					 enlargeDIMwithOneLetter=true;
			 }

			 if((enlargeDIMwithOneLetter)&&(!((i==(count-1))&&(best_child_index==-1))))
				 continue;
		 }
         old_overlap = 0.0;
         new_overlap = 0.0;
         for(j = 0; j < count; j++) // cal all overlap for cur child
            if(j != cur_child){
               old_overlap += cal_normed_overlap(entriesVec.at(j).DMBR, entriesVec.at(cur_child).DMBR, alphabet_sizes);
               new_overlap += cal_normed_overlap(entriesVec.at(j).DMBR, tmp_DMBR, alphabet_sizes);
            }

         cur_overlap_enlarge = new_overlap - old_overlap;
         if((best_overlap_enlarge - cur_overlap_enlarge) > DOUBLE_THRESHOLD){ // cur_overlap_enlarge < best_overlap_enlarge
            best_block_no = entriesVec.at(cur_child).child;
            best_child_index = cur_child;
            best_overlap_enlarge = cur_overlap_enlarge;
            orig_area = cal_normed_area(entriesVec.at(cur_child).DMBR, alphabet_sizes);
            new_area = cal_normed_area(tmp_DMBR, alphabet_sizes);
            best_enlarge = new_area - orig_area;
            best_area = orig_area;
					lastHeuristicUsed = 3;
#ifdef LOG_WHICH_HEURISTICS_CHOOSE_LEAF
			if(debug_height==2)
				logO.log2File(" overlap_enlarge ");
#endif
         }else if (((best_overlap_enlarge - cur_overlap_enlarge) > 0 ? best_overlap_enlarge - cur_overlap_enlarge : cur_overlap_enlarge - best_overlap_enlarge) <= DOUBLE_THRESHOLD){ // cur_overlap_enlarge == best_overlap_enlarge 
            orig_area = cal_normed_area(entriesVec.at(cur_child).DMBR, alphabet_sizes);
            new_area = cal_normed_area(tmp_DMBR, alphabet_sizes);
            cur_enlarge = new_area - orig_area;

lastHeuristicUsed = 5;
            if((best_enlarge - cur_enlarge) > DOUBLE_THRESHOLD){ // cur_enlarge < best_enlarge
               best_block_no = entriesVec.at(cur_child).child;
               best_child_index = cur_child;
               best_enlarge = cur_enlarge;
               best_area = orig_area;
					
#ifdef LOG_WHICH_HEURISTICS_CHOOSE_LEAF
			if(debug_height==2)
				logO.log2File(" area_enlarge ");
#endif
            }else if(((best_enlarge - cur_enlarge) > 0 ? best_enlarge - cur_enlarge : cur_enlarge - best_enlarge) <= DOUBLE_THRESHOLD){ // cur_enlarge == best_enlarge
lastHeuristicUsed = 7;
               if((best_area - orig_area) > DOUBLE_THRESHOLD){ //orig_area < best_area
                  best_block_no = entriesVec.at(cur_child).child;
                  best_child_index = cur_child;
                  best_area = orig_area;
					
#ifdef LOG_WHICH_HEURISTICS_CHOOSE_LEAF
			if(debug_height==2)
				logO.log2File(" area ");
#endif
               }
            }
         }
      } // for

	  delete []enlargement;
	  delete []sorted_child_list;
   }



   child_block_no = best_block_no;
   child_index = best_child_index;

	switch(lastHeuristicUsed)
	{
	case 1:
		BestChild_covered_area++;
		break;
	case 3:
		BestChild_notcovered_overlap_enlarge++;
		break;

	case 5:
		BestChild_notcovered_area_enlarge++;
		break;
	case 7:
		BestChild_notcovered_area++;
	}









}



void Dir_node::set_DMBR(int entry_index, unsigned char* new_DMBR, bool& not_same){
   int i;

   not_same = false;
   for(i = 0; i < DMBR_SIZE; i++){
      if(entriesVec.at(entry_index).DMBR[i] != new_DMBR[i]){ // not the same
         entriesVec.at(entry_index).DMBR[i] = new_DMBR[i];
         not_same = true;
      }
   }
   if(not_same)cal_DMBR(entriesVec, count, new_DMBR);
}



Error_code Dir_node::insert_new_entry(unsigned int new_entry_block_no, unsigned char* new_entry_DMBR, double dir_min_util, int *alphabet_sizes, unsigned char* cur_DMBR, Dir_node*& new_dir_node, unsigned char* new_DMBR){
   int i;

	  Dir_entry tempEntry;
      tempEntry.child = new_entry_block_no;
      for(i = 0; i < DMBR_SIZE; i++)tempEntry.DMBR[i] = new_entry_DMBR[i];
	bool noOverFlow;

	if(TREE_TYPE == STATIC_TREE)
		noOverFlow=count < DIR_NODE_SIZE?true:false;
	else
	{
	
	  vector<Dir_entry> tmp_vec(entriesVec);
	  tmp_vec.push_back(tempEntry);

	  int compressedSize = getCompressedDiskSize(tmp_vec,alphabet_sizes);
		noOverFlow=(compressedSize<= DISK_BLOCK_SIZE)?true:false;
	}
   //if(count < DIR_NODE_SIZE){ // room avaialbe
	if(noOverFlow)
	{
      //entriesVec.at(count).child = new_entry_block_no;
      //for(i = 0; i < DMBR_SIZE; i++)entriesVec.at(count).DMBR[i] = new_entry_DMBR[i];
      //count++;
      //// create a new DMBR
      //cal_DMBR(entriesVec, count, cur_DMBR);
      //return success;


	  entriesVec.push_back(tempEntry);
      count++;
      // create a new DMBR
      cal_DMBR(entriesVec, count, cur_DMBR);
      return success;
   }else
   { // overflows, split
	   //Dir_entry tmp_entries[DIR_NODE_SIZE+1];

	   Dir_entry *tmp_entries=new Dir_entry[count+1];
      for(i = 0; i < count; i++)tmp_entries[i] = entriesVec.at(i);
	  tmp_entries[count]=tempEntry;


      //int entry_group1[DIR_NODE_SIZE], count1, entry_group2[DIR_NODE_SIZE], count2;
      int *entry_group1, count1, *entry_group2, count2;
	  entry_group1 = new int[count];
	  entry_group2 = new int[count];

		switch (nodeSplitType)
		{
		case ORIGINAL:
			split_algorithm(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
			break;
		//case TO_DEATH:
		//	split_algorithm_TO_DEATH(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
		//	break;
		//case TO_DEATH_MAX_UTIL:
		//	split_TO_DEATH_MAX_UTIL(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
		//	break;
		case TO_DEATH_RANDOM_UTIL:
			{
			//split_TO_DEATH_RANDOM_UTIL(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
		if(enforce_minUtil_for_exhaustSplit==0)
			split_TO_DEATH_RANDOM_UTIL_v2(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
		else
			if(try_all_dim_with_minUtil==0)
				split_TO_DEATH_RANDOM_UTIL_v3(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR
			else
				split_TO_DEATH_RANDOM_UTIL_v4(tmp_entries, dir_min_util, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR


			if(RESORT_2_OLD_SPLIT==true)
			{
				float split_overlap = cal_normed_overlap(new_DMBR, cur_DMBR, alphabet_sizes);
				if(split_overlap!=0)
				{
					split_algorithm_after_knapsack(tmp_entries, 0.3, alphabet_sizes, entry_group1, count1, entry_group2, count2, new_DMBR, cur_DMBR); // group1 corresponds new_DMBR. group2 corresponds cur_DMBR

				}
			}
			break;
			}
		default:
			cout<<"nodeSplitType not defined"<<endl;
			exit(0);

		}


      //assert((count1 + count2) == DIR_NODE_SIZE+1);
      assert((count1 + count2) == count+1);
		assert(count1>0);assert(count2>0);

      // create a new dir node for group1
      new_dir_node = new Dir_node(alphabet_sizes);
      new_dir_node->set_node_count(count1);
      for(i = 0; i < count1; i++)
         new_dir_node->set_node_entry(i, tmp_entries[entry_group1[i]]);

      //count = count2;
	  set_node_count(count2);

      for(i = 0; i < count2; i++)
         entriesVec.at(i) = tmp_entries[entry_group2[i]]; // insert an entry knowing it will not overflow

	delete []tmp_entries;	
	delete []entry_group1;
	delete []entry_group2;



      return overflow;
   }
}



void Dir_node::set_node_count(int new_count){
   assert(count >= 0);
   count = new_count;
   entriesVec.resize(new_count);
}



void Dir_node::set_node_entry(int index, Dir_entry& new_entry){
   assert(index >= 0 && index < count);
   entriesVec.at(index) = new_entry;
}



void Dir_node::set_node_entry(int index, unsigned int new_entry_block_no, unsigned char* new_entry_DMBR){
   assert(index >= 0 && index < count);
   entriesVec.at(index).child = new_entry_block_no;
   for(int i = 0; i < DMBR_SIZE; i++)
      entriesVec.at(index).DMBR[i] = new_entry_DMBR[i];
}



int Dir_node::get_node_count(){
		assert (count==entriesVec.size());
   return count;
}


void Dir_node::get_node_entry(int index, Dir_entry& cur_entry){
   assert(index >= 0 && index < count);
   cur_entry = entriesVec.at(index);
}



void Dir_node::get_DMBR(unsigned char* cur_DMBR){
   int i, j;

   for(i = 0; i < DMBR_SIZE; i++)
      cur_DMBR[i] = 0;
   for(i = 0; i < count; i++)
      for(j = 0; j < DMBR_SIZE; j++)
         cur_DMBR[j] |= entriesVec.at(i).DMBR[j];
}

void Dir_node::get_DMBR(unsigned char* cur_DMBR,int index)
{
	assert(index<entriesVec.size());
   int i, j;

   for(i = 0; i < DMBR_SIZE; i++)
	   cur_DMBR[i] = 0;
   for(j = 0; j < DMBR_SIZE; j++)
	   cur_DMBR[j] |= entriesVec.at(index).DMBR[j];
}




void Dir_node::read_node_dynamic(fstream& ND_file, unsigned int block_number){
	unsigned char block[DISK_BLOCK_SIZE];
	unsigned char * block_ptr=block;
	unsigned char entryHeader[usedBytesForEntryHeader];//enough bytes for entry header;

   ND_file.seekg(block_number*DISK_BLOCK_SIZE);

   ND_file.read((char*)block, DISK_BLOCK_SIZE);
assert(!ND_file.fail());

//count=*(int*)block_ptr;
	memcpy((void *)(&count), (void *)block_ptr, sizeof(int));
	block_ptr+=sizeof(int);

	entriesVec.resize(count);
	for(int e=0;e<count;e++)
	{
		memcpy((void *)(&(entriesVec.at(e).child)), (void *)block_ptr, sizeof(unsigned int));
//entriesVec.at(e).child = *(unsigned int *)block_ptr;
		block_ptr+=sizeof(unsigned int);

		memcpy((void *)entryHeader, (void *)block_ptr, usedBytesForEntryHeader);
		block_ptr+=usedBytesForEntryHeader;
		for(int i =0;i<DIM;i++)
		{
			if(entryHeader[(int)floor((float)i/8)]&bit_mask[7-i%8])//this dimension has bit 0 - not compressed
			{
				for(int t=0;t<BYTES_PER_DIM_IN_DMBR;t++)
					entriesVec.at(e).DMBR[i*BYTES_PER_DIM_IN_DMBR+t]=0x00;

				for(int t=0;t<A[i];t++)
				{
					int byte_no = this->DMBR_byte_lut[i][t];
					int bit_no = this->DMBR_bit_lut[i][t];
					entriesVec.at(e).DMBR[byte_no] |= MASKS[bit_no];
				}

			}
			else
			{
				memcpy((void *)(&(entriesVec.at(e).DMBR[i*BYTES_PER_DIM_IN_DMBR])),(void *)block_ptr, BYTES_PER_DIM_IN_DMBR);
				block_ptr+=BYTES_PER_DIM_IN_DMBR;
			}
		}


	}



}

void Dir_node::write_node_dynamic(fstream& ND_file, unsigned int block_number){

	unsigned char block[DISK_BLOCK_SIZE];
	unsigned char * block_ptr=block;
	
	unsigned char entryHeader[usedBytesForEntryHeader];//enough bytes for entry header;

	ND_file.seekp(block_number * DISK_BLOCK_SIZE);
	

/////////////debug
//if(block_number==10)
//	block_number=block_number;

////////////



	memcpy((void *)block_ptr, (void *)(&count), sizeof(int));
	block_ptr+=sizeof(int);

	for(int e = 0; e < count; e++)
	{
		memcpy((void *)block_ptr, (void *)(&(entriesVec.at(e).child)), sizeof(unsigned int));

		block_ptr+=sizeof(unsigned int);

		for(int x=0;x<usedBytesForEntryHeader;x++)
			entryHeader[x]=0x00;//INIT

		vector<int> unFullDim;

		for(int i =0;i<DIM;i++)
		{
			int tmp_sum = 0;
			for(int j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
				tmp_sum += bit_1_count_lut[entriesVec.at(e).DMBR[j]];

			if(tmp_sum==A[i])//a full dimension
			{
				entryHeader[(int)floor((float)i/8)]^=bit_mask[7-i%8];

				//compressedDMBRByteNum-=BYTES_PER_DIM_IN_DMBR;
			}
			else
				unFullDim.push_back(i);

		}


		memcpy((void *)block_ptr, (void *)entryHeader, usedBytesForEntryHeader);
		block_ptr+=usedBytesForEntryHeader;

		for(int t=0;t<unFullDim.size();t++,block_ptr+=BYTES_PER_DIM_IN_DMBR)
			memcpy((void *)block_ptr, (void *)(&(entriesVec.at(e).DMBR[unFullDim.at(t)*BYTES_PER_DIM_IN_DMBR])), BYTES_PER_DIM_IN_DMBR);
	}

	ND_file.write(reinterpret_cast<const char*>(block), DISK_BLOCK_SIZE);

	assert(!ND_file.fail());
}




//void Dir_node::write_node_dynamic(fstream& ND_file, unsigned int block_number){
//
//	unsigned char entryHeader[usedBytesForEntryHeader];//enough bytes for entry header;
//
//	ND_file.seekp(block_number * DISK_BLOCK_SIZE);
//
/////*debug*/
////cout<<"tellp"<<ND_file.tellp()<<endl;
//
//	ND_file.write(reinterpret_cast<const char*>(&count), sizeof(int));
//
//	//unsigned int child_array[DIR_NODE_SIZE]; 
//	//unsigned char DMBR_array[DIR_NODE_SIZE][DMBR_SIZE];
//
//	//for(int i = 0; i < count; i++){
//	//	child_array[i] = entries[i].child;
//	//	//for(int j =0; j < DMBR_SIZE; j++)
//	//	//   DMBR_array[i][j] = entries[i].DMBR[j];
//	//}
//
//	//ND_file.write((const char*)child_array, sizeof(child_array));
//
//	for(int e = 0; e < count; e++)
//	{
//		ND_file.write((const char*)(&(entriesVec.at(e).child)), sizeof(unsigned int));
//
//		for(int x=0;x<usedBytesForEntryHeader;x++)
//			entryHeader[x]=0x00;//INIT
//
//		vector<int> unFullDim;
//
//		for(int i =0;i<DIM;i++)
//		{
//			int tmp_sum = 0;
//			for(int j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
//				tmp_sum += bit_1_count_lut[entriesVec.at(e).DMBR[j]];
//
//			if(tmp_sum==A[i])//a full dimension
//			{
//				entryHeader[(int)floor((float)i/8)]^=bit_mask[7-i%8];
//
//				//compressedDMBRByteNum-=BYTES_PER_DIM_IN_DMBR;
//			}
//			else
//				unFullDim.push_back(i);
//
//		}
//		ND_file.write((const char*)entryHeader, usedBytesForEntryHeader);
//
//
//		for(int t=0;t<unFullDim.size();t++)
//			ND_file.write((const char*)(&(entriesVec.at(e).DMBR[unFullDim.at(t)*BYTES_PER_DIM_IN_DMBR])), BYTES_PER_DIM_IN_DMBR);
//
//	//					/////////////debug
//	//cout<<unFullDim.size()<<"skipped"<<endl;
//	//			////////////
//	}
//
//
///////////////////////////////debug
////		cout<<"is failed? "<<ND_file.fail()<<endl;
////
//////////////////////////////
//assert(!ND_file.fail());
//}
//
//
void Dir_node::read_node_static(fstream& ND_file, unsigned int block_number){
//	clock_t start, finish;
//   start = clock();

   ND_file.seekg(block_number*DISK_BLOCK_SIZE);
   ND_file.read((char*)(&count), sizeof(int));

   unsigned int child_array[DIR_NODE_SIZE]; 
   unsigned char DMBR_array[DIR_NODE_SIZE][DMBR_SIZE];

   
   ND_file.read((char*)child_array, sizeof(child_array));
   ND_file.read((char*)DMBR_array, sizeof(DMBR_array));

//   finish = clock();
//   cout << (double)(finish - start) / CLOCKS_PER_SEC << " ";
entriesVec.resize(count);
   for(int i = 0; i < count; i++){
      entriesVec.at(i).child = child_array[i];
      for(int j =0; j < DMBR_SIZE; j++)
         //entries[i].DMBR[j] = DMBR_array[i][j];
		 entriesVec.at(i).DMBR[j] = DMBR_array[i][j];
   }
}



void Dir_node::write_node_static(fstream& ND_file, unsigned int block_number){
   ND_file.seekp(block_number * DISK_BLOCK_SIZE);
   ND_file.write(reinterpret_cast<const char*>(&count), sizeof(int));

   unsigned int child_array[DIR_NODE_SIZE]; 
   unsigned char DMBR_array[DIR_NODE_SIZE][DMBR_SIZE];

   for(int i = 0; i < count; i++){
      child_array[i] = entriesVec.at(i).child;
      for(int j =0; j < DMBR_SIZE; j++)
         DMBR_array[i][j] = entriesVec.at(i).DMBR[j];
   }

   ND_file.write((const char*)child_array, sizeof(child_array));
   ND_file.write((const char*)DMBR_array, sizeof(DMBR_array));
}

bool Dir_node::is_covered(const Leaf_entry &query_data, unsigned char* DMBR){
   int i;
   int byte_no, bit_no;

   for(i = 0; i < DIM; i++){
      byte_no = this->DMBR_byte_lut[i][query_data.key[i]];
      bit_no = this->DMBR_bit_lut[i][query_data.key[i]];
      if(!(DMBR[byte_no] & MASKS[bit_no]))return false; // one dim not covered by current entry, try next
   }
   return true;
}



bool Dir_node::is_within_box(
	const Dir_entry box_query_data,	
	unsigned char* DMBR)
{
	int j,i;
	bool matchOnThisDim;


	for(j = 0; j < DMBR_SIZE; j=j+BYTES_PER_DIM_IN_DMBR)
	{
		matchOnThisDim=false;

		for(i=0;i<BYTES_PER_DIM_IN_DMBR;i++)
			if((DMBR[j+i] & box_query_data.DMBR[j+i])!= 0)
				matchOnThisDim = true;

		if (matchOnThisDim == false)
			return false;
	}


	return true;
}



bool Dir_node::is_within_Hamming_dist(const Leaf_entry &query_data, unsigned char* DMBR, int range){
   int i;
   int byte_no, bit_no;

   int cur_dist = 0;
   for(i = 0; i < DIM && cur_dist <= range; i++){
      byte_no = this->DMBR_byte_lut[i][query_data.key[i]];
      bit_no = this->DMBR_bit_lut[i][query_data.key[i]];
      if(!(DMBR[byte_no] & MASKS[bit_no]))cur_dist++; 
   }
   if(cur_dist <= range)return true;
   else return false;
}



void Dir_node::enlarge_DMBR(unsigned char* DMBR, const Leaf_entry &new_data){
   int i;
   int byte_no, bit_no;
   for(i = 0; i < DIM; i++){ // union the new data -- enlarge the DMBR
      byte_no = this->DMBR_byte_lut[i][new_data.key[i]];
      bit_no = this->DMBR_bit_lut[i][new_data.key[i]];
      DMBR[byte_no] |= MASKS[bit_no];
   }
}



void Dir_node::cal_DMBR(vector<Dir_entry> & dir_entries, int num_of_dir_entries, unsigned char* new_DMBR){
   int i, j;

   // Initialize
   for(i = 0; i < DMBR_SIZE; i++)new_DMBR[i] = dir_entries.at(0).DMBR[i];

   // create the DMBR
   for(i = 0; i < DMBR_SIZE; i++)
      for(j = 1; j < num_of_dir_entries; j++)
         new_DMBR[i] |= dir_entries.at(j).DMBR[i];
}



void Dir_node::cal_DMBR(vector<Dir_entry> &  dir_entries, int* indices_of_dir_entries_used, int num_of_dir_entries_used, unsigned char* new_DMBR){
   int i, j;
   int tmp_index;

   // Initialize
   tmp_index = indices_of_dir_entries_used[0];
   for(i = 0; i < DMBR_SIZE; i++)new_DMBR[i] = dir_entries.at(tmp_index).DMBR[i];

   // create the DMBR
   for(i = 0; i < DMBR_SIZE; i++)
      for(j = 1; j < num_of_dir_entries_used; j++){
         tmp_index = indices_of_dir_entries_used[j];
         new_DMBR[i] |= dir_entries.at(tmp_index).DMBR[i];
      }
}




void Dir_node::cal_DMBR(Dir_entry* dir_entries, int num_of_dir_entries, unsigned char* new_DMBR){
   int i, j;

   // Initialize
   for(i = 0; i < DMBR_SIZE; i++)new_DMBR[i] = dir_entries[0].DMBR[i];

   // create the DMBR
   for(i = 0; i < DMBR_SIZE; i++)
      for(j = 1; j < num_of_dir_entries; j++)
         new_DMBR[i] |= dir_entries[j].DMBR[i];
}



void Dir_node::cal_DMBR(Dir_entry* dir_entries, int* indices_of_dir_entries_used, int num_of_dir_entries_used, unsigned char* new_DMBR){
   int i, j;
   int tmp_index;

   // Initialize
   tmp_index = indices_of_dir_entries_used[0];
   for(i = 0; i < DMBR_SIZE; i++)new_DMBR[i] = dir_entries[tmp_index].DMBR[i];

   // create the DMBR
   for(i = 0; i < DMBR_SIZE; i++)
      for(j = 1; j < num_of_dir_entries_used; j++){
         tmp_index = indices_of_dir_entries_used[j];
         new_DMBR[i] |= dir_entries[tmp_index].DMBR[i];
      }
}




//derived from split_TO_DEATH_RANDOM_UTIL_v3
// sort dimensions
//split all dim with more than 1 letter
//Knapsack problem
void Dir_node::split_TO_DEATH_RANDOM_UTIL_v4(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;
   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   //int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   //if(split_start < 0) split_start = 0;
   //int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   //if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   //if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif
   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in ascending order using span
   NDT_stable_sort(DIM, span_list, sorted_dim_list, true);
	//for(i= 0; i < DIM; i++)
	//	sorted_dim_list[i]=i;


   int firstIndexWithMoreThan1Letters;

   for(i = 0; i < DIM; i++)
   {
	 if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
	 {
		firstIndexWithMoreThan1Letters=i;
		break;
	 }
   }

   if(firstIndexWithMoreThan1Letters==DIM)
		firstIndexWithMoreThan1Letters=DIM-1;


	double best_overlap=DOUBLE_INF, cur_overlap;
	int best_leftNumLetters=A[0]*10, cur_leftNumLetters;
	int cur_cutPoint;
	float best_span = DOUBLE_INF;
   for(i = firstIndexWithMoreThan1Letters; i < DIM; i++)
   {
	   int cur_dim=sorted_dim_list[i];
	   if(span_list[i]>best_span)
		   continue;
	   else
		   best_span=span_list[i];

//////////////////////debug
	   //cout<<"-----------------------\n";
	   //for(int debugp=0;debugp<count+1;debugp++)
	   //{
		  // print_OneEntry_OnOneDscDim(&split_entries[debugp ],debugp,cur_dim, A);

	   //}
///////////////////////


	   sort_entries_by_size_Random_Util_v3(dir_min_util,cur_dim,split_entries, num_of_entries, sorted_entry_list,cur_leftNumLetters,cur_overlap,cur_cutPoint);

	   bool getBetterOne=false;
	   if(cur_overlap<best_overlap)
	   {
		   getBetterOne = true;
		   best_overlap=cur_overlap;
		   best_leftNumLetters = cur_leftNumLetters;

	   }
	   //do not need else because if there's a tie on overlap the next criteria should be the shorter dim, since dim already sorted by span, no need of else
	   ////else if((cur_overlap<=DOUBLE_THRESHOLD)&&(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD))
		  //// if(cur_leftNumLetters<best_leftNumLetters)
		  //// {
			 ////  getBetterOne = true;
			 ////  best_overlap=cur_overlap;
			 ////  best_leftNumLetters = cur_leftNumLetters;

		  //// }

	   if(getBetterOne)
	   {

 		 count1 = cur_cutPoint;
         count2 = num_of_entries - count1;

		   for(j=0;j<num_of_entries; j++)
		   {
            if(j < count1)
				group1[j] = sorted_entry_list[j];
            else 
				group2[j - count1] = sorted_entry_list[j];
		   }

	   }
			
   }//end of for


   for(j=0;j<count1; j++)
	   sorted_entry_list[j]=  group1[j] ;
   for(j=0;j<count2; j++)
	   sorted_entry_list[count1+j]=  group2[j] ;

   cal_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
   cal_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2


////////////////////debug
//cout<<"=======================\n";
//for (int di=0;di<(count+1);di++)
//{
//	cout<<sorted_entry_list[di]<<endl;
//
//}
//////////////////


   delete []sorted_entry_list;


#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif

}











//derived from split_TO_DEATH_RANDOM_UTIL_v2
//do not sort dimensions
//only try to split on first dim with more than 1 letter
//Knapsack problem
void Dir_node::split_TO_DEATH_RANDOM_UTIL_v3(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;
   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   //int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   //if(split_start < 0) split_start = 0;
   //int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   //if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   //if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif
   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in descending order using span
   //NDT_stable_sort(DIM, span_list, sorted_dim_list, true);
	for(i= 0; i < DIM; i++)
		sorted_dim_list[i]=i;


   int firstIndexWithMoreThan1Letters;

   for(i = 0; i < DIM; i++)
   {
	 if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
	 {
		firstIndexWithMoreThan1Letters=i;
		break;
	 }
   }

   if(firstIndexWithMoreThan1Letters==DIM)
		firstIndexWithMoreThan1Letters=DIM-1;


	double best_overlap=DOUBLE_INF, cur_overlap;
	int best_leftNumLetters=A[0]*10, cur_leftNumLetters;
	int cur_cutPoint;
   //for(i = firstIndexWithMoreThan1Letters; i < DIM; i++)
   //{
	   int cur_dim=sorted_dim_list[firstIndexWithMoreThan1Letters];

//////////////////////debug
	   //cout<<"-----------------------\n";
	   //for(int debugp=0;debugp<count+1;debugp++)
	   //{
		  // print_OneEntry_OnOneDscDim(&split_entries[debugp ],debugp,cur_dim, A);

	   //}
///////////////////////


	   sort_entries_by_size_Random_Util_v3(dir_min_util,cur_dim,split_entries, num_of_entries, sorted_entry_list,cur_leftNumLetters,cur_overlap,cur_cutPoint);

	   bool getBetterOne=false;
	   if(cur_overlap<best_overlap)
	   {
		   getBetterOne = true;
		   best_overlap=cur_overlap;
		   best_leftNumLetters = cur_leftNumLetters;

	   }
	   ////else if((cur_overlap<=DOUBLE_THRESHOLD)&&(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD))
		  //// if(cur_leftNumLetters<best_leftNumLetters)
		  //// {
			 ////  getBetterOne = true;
			 ////  best_overlap=cur_overlap;
			 ////  best_leftNumLetters = cur_leftNumLetters;

		  //// }

	   if(getBetterOne)
	   {

 		 count1 = cur_cutPoint;
         count2 = num_of_entries - count1;

		   for(j=0;j<num_of_entries; j++)
		   {
            if(j < count1)
				group1[j] = sorted_entry_list[j];
            else 
				group2[j - count1] = sorted_entry_list[j];
		   }

	   }
			
   //}

   for(j=0;j<count1; j++)
	   sorted_entry_list[j]=  group1[j] ;
   for(j=0;j<count2; j++)
	   sorted_entry_list[count1+j]=  group2[j] ;

   cal_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
   cal_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2

   delete []sorted_entry_list;





#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif

}






void Dir_node::split_TO_DEATH_RANDOM_UTIL_v2(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;
   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   //int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   //if(split_start < 0) split_start = 0;
   //int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   //if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   //if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif
   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in descending order using span
   //NDT_sort(DIM, span_list, sorted_dim_list, false);
   NDT_stable_sort(DIM, span_list, sorted_dim_list, true);

  // for(int x=0;x<DIM;x++)
		//sorted_dim_list[x]=x;

   int firstIndexWithMoreThan1Letters;

   for(i = 0; i < DIM; i++)
   {
	 if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
	 {
		firstIndexWithMoreThan1Letters=i;
		break;
	 }
   }

   if(firstIndexWithMoreThan1Letters==DIM)
		firstIndexWithMoreThan1Letters=DIM-1;


	double best_overlap=DOUBLE_INF, cur_overlap;
	int best_leftNumLetters=A[0]*10, cur_leftNumLetters;
	int cur_cutPoint;
   for(i = firstIndexWithMoreThan1Letters; i < DIM; i++)
   {
	   int cur_dim=sorted_dim_list[i];

//////////////////////debug
	   //cout<<"-----------------------\n";
	   //for(int debugp=0;debugp<count+1;debugp++)
	   //{
		  // print_OneEntry_OnOneDscDim(&split_entries[debugp ],debugp,cur_dim, A);

	   //}
///////////////////////


	   sort_entries_by_size_Random_Util_v2(cur_dim,split_entries, num_of_entries, sorted_entry_list,cur_leftNumLetters,cur_overlap,cur_cutPoint);

	   bool getBetterOne=false;
	   if(cur_overlap<best_overlap)
	   {
		   getBetterOne = true;
		   best_overlap=cur_overlap;
		   best_leftNumLetters = cur_leftNumLetters;

	   }
	   else if((cur_overlap<=DOUBLE_THRESHOLD)&&(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD))
	   //else if((cur_overlap<=DOUBLE_THRESHOLD)&&(cur_overlap < best_overlap))
		   if(cur_leftNumLetters<best_leftNumLetters)
		   {
			   getBetterOne = true;
			   best_overlap=cur_overlap;
			   best_leftNumLetters = cur_leftNumLetters;

		   }

	   if(getBetterOne)
	   {

 		 count1 = cur_cutPoint;
         count2 = num_of_entries - count1;

		   for(j=0;j<num_of_entries; j++)
		   {
            if(j < count1)
				group1[j] = sorted_entry_list[j];
            else 
				group2[j - count1] = sorted_entry_list[j];
		   }

	   }
			
   }

   for(j=0;j<count1; j++)
	   sorted_entry_list[j]=  group1[j] ;
   for(j=0;j<count2; j++)
	   sorted_entry_list[count1+j]=  group2[j] ;

   cal_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
   cal_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2

   delete []sorted_entry_list;





#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif

}




void Dir_node::split_TO_DEATH_RANDOM_UTIL(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;
   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   if(split_start < 0) split_start = 0;
   int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif
   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in descending order using span
   //NDT_sort(DIM, span_list, sorted_dim_list, false);
   NDT_stable_sort(DIM, span_list, sorted_dim_list, true);

  // for(int x=0;x<DIM;x++)
		//sorted_dim_list[x]=x;

   int firstIndexWithMoreThan1Letters;

   for(i = 0; i < DIM; i++)
   {
	 if(span_list[i]>static_cast<double>(1.0) / alphabet_sizes[sorted_dim_list[i]])//if this dim has more than 1 letter
	 {
		firstIndexWithMoreThan1Letters=i;
		break;
	 }
   }

   if(firstIndexWithMoreThan1Letters==DIM)
		firstIndexWithMoreThan1Letters=DIM-1;





	vector<int> shortestDim,longerDim;
   shortestDim.push_back(sorted_dim_list[firstIndexWithMoreThan1Letters]);
   for(i = firstIndexWithMoreThan1Letters+1; i < DIM; i++)
   {
		if(span_list[i]==span_list[firstIndexWithMoreThan1Letters])
			shortestDim.push_back(sorted_dim_list[i]);
		else
			longerDim.push_back(sorted_dim_list[i]);
   }

   int best_util_dim;
   int most_entries=-1;
   for(i=0;i<shortestDim.size();i++)//first round, split to 1 and rest on shortest dims
   {
	  int num_left_entries=sort_entries_by_size_RANDOM_UTIL(shortestDim.at(i), split_entries, num_of_entries, sorted_entry_list);
	  if(num_left_entries>0)
	  {
		best_util_dim=shortestDim.at(i);
		most_entries=num_left_entries;
		break;
	  }
   }

   if (most_entries<=0)//first round, split to 1 and rest on longer dims
   {
	   for (int t1=0;t1<longerDim.size();t1++)
	   {
		   int num_left_entries=sort_entries_by_size_RANDOM_UTIL(longerDim.at(i), split_entries, num_of_entries, sorted_entry_list);
		   if(num_left_entries>0)
		   {
			   best_util_dim=longerDim.at(i);
			   most_entries=num_left_entries;
			   break;
		   }
	   }
   }

   if(most_entries<=0)//third round, force to split on shortest dim
   {
		most_entries=sort_entries_by_size_RANDOM_UTIL(shortestDim.at(0), split_entries, num_of_entries, sorted_entry_list,true);
		best_util_dim=shortestDim.at(0);
   }
   ////splitting dim is best_util_dim
   // //# of entries on left side is most_entries

		 //most_entries=sort_entries_by_size_RANDOM_UTIL(best_util_dim, split_entries, num_of_entries, sorted_entry_list);

		 count1 = most_entries;
         count2 = num_of_entries - count1;
         for(j = 0; j < num_of_entries; j++){
            if(j < count1)group1[j] = sorted_entry_list[j];
            else group2[j - count1] = sorted_entry_list[j];
         }
         cal_DMBR(split_entries, sorted_entry_list, count1, DMBR1); // create DMBR for grp1
         cal_DMBR(split_entries, sorted_entry_list + count1, num_of_entries - count1, DMBR2); // create DMBR for grp2

   
	delete []sorted_entry_list;
#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif

}

void Dir_node::split_algorithm_TO_DEATH(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;

   int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   if(split_start < 0) split_start = 0;
   int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif
   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in descending order using span
   //NDT_sort(DIM, span_list, sorted_dim_list, false);
   NDT_stable_sort(DIM, span_list, sorted_dim_list, true);

  // for(int x=0;x<DIM;x++)
		//sorted_dim_list[x]=x;

   double best_overlap, cur_overlap; 
   double best_span, cur_span;
   double best_balance, cur_balance; // the differenece of span of two group
   double best_area, cur_area; // the sum of area of the current best split
   int best_dim, cur_dim; // the dimension where to split
   unsigned char cur_DMBR1[DMBR_SIZE], cur_DMBR2[DMBR_SIZE];
   int tmp_sum1, tmp_sum2;

   best_overlap = DOUBLE_INF; // max possible value
   best_span = 0;
   best_balance = DOUBLE_INF; // max possible value
   best_area = DOUBLE_INF; // max possible value
   best_dim = 0; 

   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   // try every dim in order sorted by span
   for(i = 0; i < DIM; i++){

      cur_dim = sorted_dim_list[i];
      cur_span = span_list[i];
	
	  //if(cur_span*alphabet_sizes[i] ==1)
		 // continue;

      //if(best_overlap < DOUBLE_THRESHOLD && (best_span - cur_span) > DOUBLE_THRESHOLD) // best_overlap == 0 & cur_span < best_span
      //   break;
    
      // sort entries on the dimension 
	  sort_entries_by_size(cur_dim, split_entries, num_of_entries, sorted_entry_list);



	  /////////////////////////////////////////////
	  //cout<<"=========================="<<endl;
	  //for(int debugp=0;debugp<num_of_entries;debugp++)
	  //{
		 // print_OneEntry_OnOneDscDim(&split_entries[sorted_entry_list[debugp]],sorted_entry_list[debugp],	cur_dim, alphabet_sizes);

	  //}



      // try every possible group
      int best_j;
      for(j = split_start; j <= split_end; j++){
         cal_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
         cal_DMBR(split_entries, sorted_entry_list + j + 1, num_of_entries - j - 1, cur_DMBR2); // create DMBR for grp2
         cur_overlap = cal_normed_overlap(cur_DMBR1, cur_DMBR2, alphabet_sizes);
         if((best_overlap - cur_overlap) > 0){ //cur_overlap < best_overlap
            best_j = j;
            for(k = 0; k < DMBR_SIZE; k++){
               DMBR1[k] = cur_DMBR1[k];
               DMBR2[k] = cur_DMBR2[k];
            }
            best_overlap = cur_overlap;
            best_span = cur_span;
            tmp_sum1 = 0;
            tmp_sum2 = 0;
            for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
               tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
               tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
            }
            // use normalized best_balance
            best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
            if(best_balance < 0)
                best_balance = - best_balance;
            best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
            best_dim = cur_dim;
         }
		 //else if(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD){ // cur_overlap == best_overlap
		 else if(false){
            if((cur_span - best_span) > DOUBLE_THRESHOLD){ // never happen since dimensions are sorted by span
               best_j = j;
               for(k = 0; k < DMBR_SIZE; k++){
                  DMBR1[k] = cur_DMBR1[k];
                  DMBR2[k] = cur_DMBR2[k];
               }
               best_span = cur_span;
               tmp_sum1 = 0;
               tmp_sum2 = 0;
               for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
                  tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                  tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
               }
               // use normalized best_balance
               best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
               if(best_balance < 0)
                  best_balance = - best_balance;
               best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
               best_dim = cur_dim;
            }else if(((cur_span - best_span) > 0 ? (cur_span - best_span) : (best_span - cur_span)) <= DOUBLE_THRESHOLD){ // cur_span == best_span
               tmp_sum1 = 0;
               tmp_sum2 = 0;
               for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
                  tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                  tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
               }
               // use normalized balance
               cur_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
               if(cur_balance < 0)cur_balance = - cur_balance;
               if((best_balance - cur_balance) > DOUBLE_THRESHOLD){ // cur_balance < best_balance
                  best_j = j;
                  for(k = 0; k < DMBR_SIZE; k++){
                     DMBR1[k] = cur_DMBR1[k];
                     DMBR2[k] = cur_DMBR2[k];
                  }
                  best_balance = cur_balance;
                  best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                  best_dim = cur_dim;
               }else if(((cur_balance - best_balance) > 0 ? (cur_balance - best_balance) : (best_balance - cur_balance)) <= DOUBLE_THRESHOLD){ // cur_balance == best_balance
                  cur_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                  if(best_area - cur_area > DOUBLE_THRESHOLD){ // cur_area < best_area
                     best_j = j;
                     for(k = 0; k < DMBR_SIZE; k++){
                        DMBR1[k] = cur_DMBR1[k];
                        DMBR2[k] = cur_DMBR2[k];
                     }
                     best_area = cur_area;
                     best_dim = cur_dim;
                  }
               }
            }
         }
      } // for j
      // remember group1 and group 2 now
      if(best_dim == cur_dim){
         count1 = best_j + 1;
         count2 = num_of_entries - count1;
         for(j = 0; j < num_of_entries; j++){
            if(j < count1)group1[j] = sorted_entry_list[j];
            else group2[j - count1] = sorted_entry_list[j];
         }
      }
   } // for dim
#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif


#ifdef LOG_ENTRY_ON_SPLITTING_DIM
logO.log2File("dir split on dim ");logO.log2File(best_dim);logO.log2File("\n");
   for(int debugp = 0; debugp < num_of_entries;debugp++)
   {//each entry holds one line in log file
   	log_OneEntry_OnOneDscDim(&split_entries[sorted_entry_list[debugp]],sorted_entry_list[debugp],	best_dim, alphabet_sizes);
    logO.log2File("\n");
   }

#endif

delete []sorted_entry_list;
}

void Dir_node::split_algorithm_after_knapsack(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2)
{
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;

   int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   if(split_start < 0) split_start = 0;
   int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif

   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in descending order using span
   NDT_sort(DIM, span_list, sorted_dim_list, false);

   double best_overlap, cur_overlap; 
   double best_span, cur_span;
   double best_balance, cur_balance; // the differenece of span of two group
   double best_area, cur_area; // the sum of area of the current best split
   int best_dim, cur_dim; // the dimension where to split
   unsigned char cur_DMBR1[DMBR_SIZE], cur_DMBR2[DMBR_SIZE];
   int tmp_sum1, tmp_sum2;

   best_overlap = DOUBLE_INF; // max possible value
   best_span = 0;
   best_balance = DOUBLE_INF; // max possible value
   best_area = DOUBLE_INF; // max possible value
   best_dim = 0; 
	int whichHeuristicUsedLast;

   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   // try every dim in order sorted by span
   for(i = 0; i < DIM; i++){

      cur_dim = sorted_dim_list[i];
      cur_span = span_list[i];
	
	  //if(cur_span*alphabet_sizes[i] ==1)
		 // continue;

      //if(best_overlap < DOUBLE_THRESHOLD && (best_span - cur_span) > DOUBLE_THRESHOLD) // best_overlap == 0 & cur_span < best_span
      //   break;
    
      // sort entries on the dimension 
      sort_entries_by_size(cur_dim, split_entries, num_of_entries, sorted_entry_list);
      // try every possible group
      int best_j;
      for(j = split_start; j <= split_end; j++){
		 
		  if((RESORT_2_OLD_SPLIT==true)&&(nodeSplitType == TO_DEATH_RANDOM_UTIL)&&(TREE_TYPE == DYNAMIC_TREE))
		  {
			  vector<Dir_entry> left, right;

			  for(int num=0;num<j + 1;num++)
				  left.push_back(split_entries[sorted_entry_list[num]]);
			  for(int num=0;num<num_of_entries - j - 1;num++)
				  right.push_back(split_entries[sorted_entry_list[num+j + 1]]);

			  bool couldCompress;
			  couldCompress=(getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE);

			  if(!couldCompress)
				  continue;
		  }


         cal_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
         cal_DMBR(split_entries, sorted_entry_list + j + 1, num_of_entries - j - 1, cur_DMBR2); // create DMBR for grp2
         cur_overlap = cal_normed_overlap(cur_DMBR1, cur_DMBR2, alphabet_sizes);
         if((best_overlap - cur_overlap) > DOUBLE_THRESHOLD){ //cur_overlap < best_overlap
			 whichHeuristicUsedLast=1;
            best_j = j;
            for(k = 0; k < DMBR_SIZE; k++){
               DMBR1[k] = cur_DMBR1[k];
               DMBR2[k] = cur_DMBR2[k];
            }
            best_overlap = cur_overlap;
            best_span = cur_span;
            tmp_sum1 = 0;
            tmp_sum2 = 0;
            for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
               tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
               tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
            }
            // use normalized best_balance
            best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
            if(best_balance < 0)
                best_balance = - best_balance;
            best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
            best_dim = cur_dim;
         }else if(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD)
		 { // cur_overlap == best_overlap
		 //else if(false){
            //if((cur_span - best_span) > DOUBLE_THRESHOLD){ // never happen since dimensions are sorted by span
			 if(false){
               best_j = j;
               for(k = 0; k < DMBR_SIZE; k++){
                  DMBR1[k] = cur_DMBR1[k];
                  DMBR2[k] = cur_DMBR2[k];
               }
               best_span = cur_span;
               tmp_sum1 = 0;
               tmp_sum2 = 0;
               for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
                  tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                  tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
               }
               // use normalized best_balance
               best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
               if(best_balance < 0)
                  best_balance = - best_balance;
               best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
               best_dim = cur_dim;
            }
			 //else if(((cur_span - best_span) > 0 ? (cur_span - best_span) : (best_span - cur_span)) <= DOUBLE_THRESHOLD)
			 else if(true)
			 { // cur_span == best_span
               tmp_sum1 = 0;
               tmp_sum2 = 0;
               for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
                  tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                  tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
               }
               // use normalized balance
               cur_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
               if(cur_balance < 0)
				   cur_balance = - cur_balance;
               //if((best_balance - cur_balance) > DOUBLE_THRESHOLD)
			   if(false)
			   { // cur_balance < best_balance
                  best_j = j;
                  for(k = 0; k < DMBR_SIZE; k++){
                     DMBR1[k] = cur_DMBR1[k];
                     DMBR2[k] = cur_DMBR2[k];
                  }
                  best_balance = cur_balance;
                  best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                  best_dim = cur_dim;
               }
			   //else if(((cur_balance - best_balance) > 0 ? (cur_balance - best_balance) : (best_balance - cur_balance)) <= DOUBLE_THRESHOLD)
			   else if(true)
			   { // cur_balance == best_balance
                  cur_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                  if(best_area - cur_area > DOUBLE_THRESHOLD){ // cur_area < best_area
					whichHeuristicUsedLast=4;
                     best_j = j;
                     for(k = 0; k < DMBR_SIZE; k++){
                        DMBR1[k] = cur_DMBR1[k];
                        DMBR2[k] = cur_DMBR2[k];
                     }
                     best_area = cur_area;
                     best_dim = cur_dim;
                  }
               }
            }
         }
      } // for j
      // remember group1 and group 2 now
      if(best_dim == cur_dim){
         count1 = best_j + 1;
         count2 = num_of_entries - count1;
         for(j = 0; j < num_of_entries; j++){
            if(j < count1)group1[j] = sorted_entry_list[j];
            else group2[j - count1] = sorted_entry_list[j];
         }
      }
   } // for dim
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif
delete []sorted_entry_list;

		switch (whichHeuristicUsedLast)
		{
		case 1:
			heuristics_overlap_used_dir++;
			break;
		case 4:
			heuristics_area_used_dir++;
			break;
		default:
			break;
		}



}




void Dir_node::split_algorithm(Dir_entry* split_entries, double dir_min_util, int* alphabet_sizes, int* group1, int &count1, int* group2, int &count2, unsigned char* DMBR1, unsigned char* DMBR2){
   int i, j, k;

   //int num_of_entries = DIR_NODE_SIZE + 1;
   int num_of_entries = count + 1;

   int split_start = static_cast<int>(floor(DIR_NODE_SIZE * dir_min_util)) - 1; // the index in the entries where possible split point starts
   if(split_start < 0) split_start = 0;
   int split_end = static_cast<int>(ceil(DIR_NODE_SIZE * (1 - dir_min_util))); // the index in the entries where possible split point ends
   if(split_end >= num_of_entries - 1)split_end = num_of_entries - 2; // index must at least ends at num_of_entries - 1
   if(split_end < split_start)split_end = split_start;

   count1 = 0;
   count2 = 0;

   // sort the dims based on span
   // create a DMBR for the entries
   unsigned char tmp_DMBR[DMBR_SIZE];
   cal_DMBR(split_entries, num_of_entries, tmp_DMBR);
#ifdef LOG_MBR_SPLITTING_DIR

   logO.log2File("dir  before ");Node::log_DMBR(tmp_DMBR);
#endif

   // calculate the span of each dimension 
   double span_list[DIM];
   for(i = 0; i < DIM; i++){
      int tmp_sum = 0;
      for(j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
         tmp_sum += bit_1_count_lut[tmp_DMBR[j]];
      span_list[i] = static_cast<double>(tmp_sum) / alphabet_sizes[i];
   }
   int sorted_dim_list[DIM]; // hold the sorted dimension list
   // sort the dimension in descending order using span
   NDT_sort(DIM, span_list, sorted_dim_list, false);

   double best_overlap, cur_overlap; 
   double best_span, cur_span;
   double best_balance, cur_balance; // the differenece of span of two group
   double best_area, cur_area; // the sum of area of the current best split
   int best_dim, cur_dim; // the dimension where to split
   unsigned char cur_DMBR1[DMBR_SIZE], cur_DMBR2[DMBR_SIZE];
   int tmp_sum1, tmp_sum2;

   best_overlap = DOUBLE_INF; // max possible value
   best_span = 0;
   best_balance = DOUBLE_INF; // max possible value
   best_area = DOUBLE_INF; // max possible value
   best_dim = 0; 

   //int sorted_entry_list[DIR_NODE_SIZE + 1];
   int *sorted_entry_list=new int[count + 1];

   // try every dim in order sorted by span
   for(i = 0; i < DIM; i++){

      cur_dim = sorted_dim_list[i];
      cur_span = span_list[i];
	
	  //if(cur_span*alphabet_sizes[i] ==1)
		 // continue;

      if(best_overlap < DOUBLE_THRESHOLD && (best_span - cur_span) > DOUBLE_THRESHOLD) // best_overlap == 0 & cur_span < best_span
         break;
    
      // sort entries on the dimension 
      sort_entries_by_size(cur_dim, split_entries, num_of_entries, sorted_entry_list);
      // try every possible group
      int best_j;
      for(j = split_start; j <= split_end; j++){
		 
		  if((RESORT_2_OLD_SPLIT==true)&&(nodeSplitType == TO_DEATH_RANDOM_UTIL)&&(TREE_TYPE == DYNAMIC_TREE))
		  {
			  vector<Dir_entry> left, right;

			  for(int num=0;num<j + 1;num++)
				  left.push_back(split_entries[sorted_entry_list[num]]);
			  for(int num=0;num<num_of_entries - j - 1;num++)
				  right.push_back(split_entries[sorted_entry_list[num+j + 1]]);

			  bool couldCompress;
			  couldCompress=(getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE);

			  if(!couldCompress)
				  continue;
		  }


         cal_DMBR(split_entries, sorted_entry_list, j + 1, cur_DMBR1); // create DMBR for grp1
         cal_DMBR(split_entries, sorted_entry_list + j + 1, num_of_entries - j - 1, cur_DMBR2); // create DMBR for grp2
         cur_overlap = cal_normed_overlap(cur_DMBR1, cur_DMBR2, alphabet_sizes);
         if((best_overlap - cur_overlap) > DOUBLE_THRESHOLD){ //cur_overlap < best_overlap
            best_j = j;
            for(k = 0; k < DMBR_SIZE; k++){
               DMBR1[k] = cur_DMBR1[k];
               DMBR2[k] = cur_DMBR2[k];
            }
            best_overlap = cur_overlap;
            best_span = cur_span;
            tmp_sum1 = 0;
            tmp_sum2 = 0;
            for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
               tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
               tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
            }
            // use normalized best_balance
            best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
            if(best_balance < 0)
                best_balance = - best_balance;
            best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
            best_dim = cur_dim;
         }else if(((cur_overlap - best_overlap) > 0 ? (cur_overlap - best_overlap) : (best_overlap - cur_overlap)) <= DOUBLE_THRESHOLD){ // cur_overlap == best_overlap
		 //else if(false){
            if((cur_span - best_span) > DOUBLE_THRESHOLD){ // never happen since dimensions are sorted by span
               best_j = j;
               for(k = 0; k < DMBR_SIZE; k++){
                  DMBR1[k] = cur_DMBR1[k];
                  DMBR2[k] = cur_DMBR2[k];
               }
               best_span = cur_span;
               tmp_sum1 = 0;
               tmp_sum2 = 0;
               for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
                  tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                  tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
               }
               // use normalized best_balance
               best_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
               if(best_balance < 0)
                  best_balance = - best_balance;
               best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
               best_dim = cur_dim;
            }else if(((cur_span - best_span) > 0 ? (cur_span - best_span) : (best_span - cur_span)) <= DOUBLE_THRESHOLD){ // cur_span == best_span
               tmp_sum1 = 0;
               tmp_sum2 = 0;
               for(k = this->DMBR_start_byte_lut[cur_dim]; k <= this->DMBR_end_byte_lut[cur_dim]; k++){
                  tmp_sum1 += bit_1_count_lut[cur_DMBR1[k]];
                  tmp_sum2 += bit_1_count_lut[cur_DMBR2[k]];
               }
               // use normalized balance
               cur_balance = static_cast<double>(tmp_sum1 - tmp_sum2) / alphabet_sizes[cur_dim];
               if(cur_balance < 0)cur_balance = - cur_balance;
               if((best_balance - cur_balance) > DOUBLE_THRESHOLD){ // cur_balance < best_balance
                  best_j = j;
                  for(k = 0; k < DMBR_SIZE; k++){
                     DMBR1[k] = cur_DMBR1[k];
                     DMBR2[k] = cur_DMBR2[k];
                  }
                  best_balance = cur_balance;
                  best_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                  best_dim = cur_dim;
               }else if(((cur_balance - best_balance) > 0 ? (cur_balance - best_balance) : (best_balance - cur_balance)) <= DOUBLE_THRESHOLD){ // cur_balance == best_balance
                  cur_area = cal_normed_area(cur_DMBR1, alphabet_sizes) + cal_normed_area(cur_DMBR2, alphabet_sizes);
                  if(best_area - cur_area > DOUBLE_THRESHOLD){ // cur_area < best_area
                     best_j = j;
                     for(k = 0; k < DMBR_SIZE; k++){
                        DMBR1[k] = cur_DMBR1[k];
                        DMBR2[k] = cur_DMBR2[k];
                     }
                     best_area = cur_area;
                     best_dim = cur_dim;
                  }
               }
            }
         }
      } // for j
      // remember group1 and group 2 now
      if(best_dim == cur_dim){
         count1 = best_j + 1;
         count2 = num_of_entries - count1;
         for(j = 0; j < num_of_entries; j++){
            if(j < count1)group1[j] = sorted_entry_list[j];
            else group2[j - count1] = sorted_entry_list[j];
         }
      }
   } // for dim
#ifdef LOG_MBR_SPLITTING_DIR
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR1);
   logO.log2File("dir  after  ");Node::log_DMBR(DMBR2);
#endif
delete []sorted_entry_list;
}


int Dir_node::sort_entries_by_size_RANDOM_UTIL(
	int sort_dim, 
	Dir_entry* sort_entries, 
	int num_of_sort_entries, 
	int* sorted_entry_list,
	bool forceSplitonThisDim){
   int i, j, k;
   int tmp_sum, tmp_cnt;
   int tmp_index;

   //set_entry set_entry_array[DIR_NODE_SIZE+1];
set_entry *set_entry_array=new set_entry[count+1];

	for(int t=0;t<(count+1);t++) set_entry_array[t].resizeSetEntry(count+1);

   int num_of_set_entry = 0;

   int cur_dim_byte_size = this->DMBR_end_byte_lut[sort_dim] - this->DMBR_start_byte_lut[sort_dim] + 1; // number of bytes for the bitmap used to represent the set
   // cal histogram of the entries
   // initalize tmp_set

   // set_entry_array
   bool found;
   for(i = 0; i < num_of_sort_entries; i++){
      // find the position of the component set of the current entry in the set_array
      found = false;
      for(j = 0; j < num_of_set_entry; j++)
         if(set_equal(sort_entries[i].DMBR+this->DMBR_start_byte_lut[sort_dim], set_entry_array[j].set, cur_dim_byte_size)){
            found = true;
            break;
         }
      if(found){
         set_entry_array[j].entry_set[set_entry_array[j].freq] = i;
         set_entry_array[j].freq++;
      }else{ // new set entry found
         tmp_sum = 0;
         for(j = 0; j < cur_dim_byte_size; j++){
            set_entry_array[num_of_set_entry].set[j] = sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
            tmp_sum += bit_1_count_lut[sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j]];
         }
         set_entry_array[num_of_set_entry].element_size = tmp_sum;
         set_entry_array[num_of_set_entry].entry_set[0] = i;
         set_entry_array[num_of_set_entry].freq = 1;

         num_of_set_entry++;
         
         // collect the letters used into the temp set
      }
   }

   // sort sets by element_size asc, break ties by freq desc.  Selection sort
   //int sorted_set_entry_list[DIR_NODE_SIZE+1];  // an array of indices into the set_entry_array
   int * sorted_set_entry_list = new int[count+1];  // an array of indices into the set_entry_array
   
   // populate
   for(i = 0; i < num_of_set_entry; i++)sorted_set_entry_list[i] = i;
   int largest, current, position;
   for(position = num_of_set_entry - 1; position > 0; position--){
      largest = 0;
      for(current = 1; current <= position; current++)
         if(set_entry_array[sorted_set_entry_list[current]].element_size > set_entry_array[sorted_set_entry_list[largest]].element_size)
			 largest = current;
      // swap
      tmp_index = sorted_set_entry_list[largest];
      sorted_set_entry_list[largest] = sorted_set_entry_list[position];
      sorted_set_entry_list[position] = tmp_index;
   }



vector<int> setOverlapFreeWithRests;
//////////////////debug
//for(int di=0;di<num_of_set_entry;di++)
//	   for(int dj = 0; dj < set_entry_array[sorted_set_entry_list[di]].freq ; dj++)
//		  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[di]].entry_set[dj]],set_entry_array[sorted_set_entry_list[di]].entry_set[dj],	sort_dim, A);

///////////////////

   //vector<int> candidatePositions;

   bool findTheOne=false;
	   int leftIndex;
   //if(set_entry_array[sorted_set_entry_list[0]].element_size!=1)//at least one set of entries with only 1 letter on this dim
	  // findTheOne=false;
   //else if (num_of_set_entry==1)//if this dim has only 1 set of letters
	   if (num_of_set_entry==1)//if this dim has only 1 set of letters
		   findTheOne=false;
   //else if(num_of_set_entry==2)//if this dim has 2 sets of single letter
   //{


	  // int index=0;


	  // for(i = 0; i < set_entry_array[sorted_set_entry_list[0]].freq ; i++)
		 //  sorted_entry_list[index++]=set_entry_array[sorted_set_entry_list[0]].entry_set[i];

	  // for(i = 0; i < set_entry_array[sorted_set_entry_list[1]].freq ; i++)
		 //  sorted_entry_list[index++]=set_entry_array[sorted_set_entry_list[1]].entry_set[i];

	  // findTheOne=true;
	  // //return set_entry_array[sorted_set_entry_list[0]].freq;
	  // leftIndex=0;
   //}
	   else
	   {
		   for (i=0;i<num_of_set_entry;i++)
		   {
			   //unsigned char tmp_set1[MAX_DMBR_DIM_SIZE];
			   unsigned char tmp_set2[MAX_DMBR_DIM_SIZE];

			   //if(set_entry_array[sorted_set_entry_list[i]].element_size>1)//if thas tried all sets with one letters, then just quit, do not handle sets with more than 1 sets
			   // break;

			   //for(int s = 0; s < cur_dim_byte_size; s++)
			   // tmp_set1[s]=set_entry_array[sorted_set_entry_list[i]].set[s];
			   for(int s = 0; s < cur_dim_byte_size; s++)
				   tmp_set2[s]=0x00;

			   for(j=0;j<num_of_set_entry;j++)
			   {

				   if(j==i)
					   continue;
				   else
					   for(int s = 0; s < cur_dim_byte_size; s++)
						   tmp_set2[s]|=set_entry_array[sorted_set_entry_list[j]].set[s];

			   }


			   bool isThisSetOverlapWithRests=true;

			   for(int s = 0; s < cur_dim_byte_size; s++)
				   if(tmp_set2[s]&set_entry_array[sorted_set_entry_list[i]].set[s])
					   isThisSetOverlapWithRests=false;
			   if (isThisSetOverlapWithRests)
				   setOverlapFreeWithRests.push_back(i);
		   }

		   //note setOverlapFreeWithRests still keep qualified sets sorted first by size asc then by frequency dsc
		   if (setOverlapFreeWithRests.size()>0)
		   {

			   if((TREE_TYPE == DYNAMIC_TREE))
			   {

				   int tmp_counter = 0;

				   while(set_entry_array[sorted_set_entry_list[setOverlapFreeWithRests.at(tmp_counter)]].element_size==1)
				   {
					   vector<Dir_entry> left, right;

					   for(int m = 0; m < num_of_set_entry; m++)
					   {
						   if(m ==setOverlapFreeWithRests.at(tmp_counter))
						   {
							   for(int n = 0; n < set_entry_array[sorted_set_entry_list[m]].freq ; n++)
								   left.push_back(sort_entries[set_entry_array[sorted_set_entry_list[m]].entry_set[n]]);
						   }
						   else 
						   {
							   for(int n = 0; n < set_entry_array[sorted_set_entry_list[m]].freq ; n++)
								   right.push_back(sort_entries[set_entry_array[sorted_set_entry_list[m]].entry_set[n]]);
						   }				   
					   }

					   assert((left.size()+right.size())==num_of_sort_entries);

					   if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
					   {//compression ok, so find the right one
						   findTheOne=true;
						   leftIndex=setOverlapFreeWithRests.at(tmp_counter);
						   break;
					   }
						
						tmp_counter++;

						if(setOverlapFreeWithRests.size()==tmp_counter)
							break;
				   }
			   }
			   else//do not need compression, so find the right one
			   {
				   if(set_entry_array[sorted_set_entry_list[setOverlapFreeWithRests.at(0)]].element_size==1)
				   {
					   findTheOne=true;
					   leftIndex=setOverlapFreeWithRests.at(0);
					   //break;
				   }
			   }

		   }

	   }//end of else


	   //findTheOne is false, forceSplitonThisDim is true, and there are sets overlap free with others
	   if ((findTheOne == false)&&(forceSplitonThisDim == true)&&(setOverlapFreeWithRests.size()>0))
	   {//if there exists set overlap free with the rest, just pick the first set could compress
		   if((TREE_TYPE == DYNAMIC_TREE))
		   {//DYNAMIC_TREE
			   int tmp_counter=0;
			   while(1)
			   {
				   vector<Dir_entry> left, right;

				   for(int m = 0; m < num_of_set_entry; m++)
				   {
					   if(m ==setOverlapFreeWithRests.at(tmp_counter))
					   {
						   for(int n = 0; n < set_entry_array[sorted_set_entry_list[m]].freq ; n++)
							   left.push_back(sort_entries[set_entry_array[sorted_set_entry_list[m]].entry_set[n]]);
					   }
					   else 
					   {
						   for(int n = 0; n < set_entry_array[sorted_set_entry_list[m]].freq ; n++)
							   right.push_back(sort_entries[set_entry_array[sorted_set_entry_list[m]].entry_set[n]]);
					   }				   
				   }

				   assert((left.size()+right.size())==num_of_sort_entries);

				   if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
				   {//compression ok, so find the right one
					   findTheOne=true;
					   leftIndex=setOverlapFreeWithRests.at(tmp_counter);
					   break;
				   }

				   tmp_counter++;

				   if(setOverlapFreeWithRests.size()==tmp_counter)
					   break;
			   }
		   }
		   else 
		   {//static tree
			   findTheOne=true;
			   leftIndex=setOverlapFreeWithRests.at(0);

		   }
	   }
	


	if(findTheOne)
	{//find
		   int index=0;

		   for(i = 0; i < set_entry_array[sorted_set_entry_list[leftIndex]].freq ; i++)
			   sorted_entry_list[index++]=set_entry_array[sorted_set_entry_list[leftIndex]].entry_set[i];

		   //candidatePositions.push_back(index);
		   for (j=0;j<num_of_set_entry;j++)
		   {
			   if (j==leftIndex)
				   continue;

			   for(i = 0; i < set_entry_array[sorted_set_entry_list[j]].freq ; i++)
			   {
				   sorted_entry_list[index++]=set_entry_array[sorted_set_entry_list[j]].entry_set[i];
				   //candidatePositions.push_back(index);
			   }

		   }
		   assert(index==num_of_sort_entries);
		   delete []set_entry_array;
		   delete []sorted_set_entry_list;

		   return set_entry_array[sorted_set_entry_list[leftIndex]].freq;
	
	}
	else
	{// cannot find, this means 
		if(forceSplitonThisDim==false)
		{
			delete []set_entry_array;
			delete []sorted_set_entry_list;

			return 0;		

		}
		else
		{//cannot find , but forced to split, at this point:
	    //is static tree but cannot find any single set overlap free with others
	    //is dynamic tree but either cannot find any single set overlap free with others or there exists but cannot compress sucessfully
		
		
		}
	
	}


   //if((!findTheOne)&&(forceSplitonThisDim==false))//if not find and not forced to find, return false
   //{
	  // delete []set_entry_array;
	  // delete []sorted_set_entry_list;

	  // return 0;
   //}
   //else
   //{
	  // if(findTheOne)//in this case doesn't matter if forceSplitonThisDim is true or false
	  // {
		 //  //if TREE_TYPE == STATIC_TREE and findTheOne==false/true
		 //  //or if TREE_TYPE == DYNAMIC_TREE and findTheOne==true
		 //  // for dynamic tree already make sure it could compress before set findTheOne to true
		 //  int index=0;

		 //  for(i = 0; i < set_entry_array[sorted_set_entry_list[leftIndex]].freq ; i++)
			//   sorted_entry_list[index++]=set_entry_array[sorted_set_entry_list[leftIndex]].entry_set[i];

		 //  //candidatePositions.push_back(index);
		 //  for (j=0;j<num_of_set_entry;j++)
		 //  {
			//   if (j==leftIndex)
			//	   continue;

			//   for(i = 0; i < set_entry_array[sorted_set_entry_list[j]].freq ; i++)
			//   {
			//	   sorted_entry_list[index++]=set_entry_array[sorted_set_entry_list[j]].entry_set[i];
			//	   //candidatePositions.push_back(index);
			//   }

		 //  }
		 //  assert(index==num_of_sort_entries);
		 //  delete []set_entry_array;
		 //  delete []sorted_set_entry_list;

		 //  return set_entry_array[sorted_set_entry_list[leftIndex]].freq;
	  // }
	  // else
	  // {//still cannot find the one,
	  //  //is static tree but cannot find any set overlap free with others
	  //  //is dynamic tree but either cannot find any set overlap free with others or there exists but cannot compress sucessfully
		 //  //note setOverlapFreeWithRests still keep qualified sets sorted first by size asc then by frequency dsc
		 //  if (setOverlapFreeWithRests.size()>0)
		 //  {//if there exists set overlap free with the rest, just pick the first set could compress

			//   while(set_entry_array[sorted_set_entry_list[setOverlapFreeWithRests.at(tmp_counter)]].element_size==1)
			//   {
			//	   vector<Dir_entry> left, right;

			//	   for(int m = 0; m < num_of_set_entry; m++)
			//	   {
			//		   if(m ==setOverlapFreeWithRests.at(tmp_counter))
			//		   {
			//			   for(int n = 0; n < set_entry_array[sorted_set_entry_list[m]].freq ; n++)
			//				   left.push_back(sort_entries[set_entry_array[sorted_set_entry_list[m]].entry_set[n]]);
			//		   }
			//		   else 
			//		   {
			//			   for(int n = 0; n < set_entry_array[sorted_set_entry_list[m]].freq ; n++)
			//				   right.push_back(sort_entries[set_entry_array[sorted_set_entry_list[m]].entry_set[n]]);
			//		   }				   
			//	   }

			//	   assert((left.size()+right.size())==num_of_sort_entries);

			//	   if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
			//	   {//compression ok, so find the right one
			//		   findTheOne=true;
			//		   leftIndex=setOverlapFreeWithRests.at(tmp_counter);
			//		   break;
			//	   }

			//	   tmp_counter++;

			//	   if(setOverlapFreeWithRests.size()==tmp_counter)
			//		   break;
			//   }
		 //  }
		 //  else // not set is overlap free with the rest
		 //  {
		 //  
		 //  
		 //  
		 //  
		 //  }

	  // }//end of else



	  // //if(!findTheOne)//if not found a suitable one, just force to split on this dime
		 // // leftIndex=0;

	  // //if((TREE_TYPE == STATIC_TREE)||findTheOne)
	  // //{
		 // // delete []set_entry_array;
		 // // delete []sorted_set_entry_list;
		 // // return set_entry_array[sorted_set_entry_list[leftIndex]].freq;
	  // //}
	  // //else
	  // //{
		 // // leftSideEntryNum=set_entry_array[sorted_set_entry_list[leftIndex]].freq;

		 // // vector<Dir_entry> left, right;				
		 // // for(int m = 0; m < num_of_sort_entries; m++)
		 // // {
			// //  if(m < count1)
			//	//   left.push_back(sort_entries[sorted_entry_list[m]]);
			// //  else 
			//	//   right.push_back(sort_entries[sorted_entry_list[m]]);
		 // // }
		 // // bothNodesCouldCompress==((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))

			// //  if(!bothNodesCouldCompress)
			// //  {


			// //  }

	  // //}


   //}

return -1;//this function is aborted
}



void Dir_node::find_minOverlap_by_aux_tree(int sort_dim, Dir_entry* sort_entries,aux_tree_node* tree, int root, int letters_bitmap_byte_size, set_entry* set_entry_array, int set_entry_array_size, vector<int> & leftSets,  double & rst_overlap)
{
	int i, j, k, l, m;
	int tmp_index;
	int tmp_start_byte, tmp_end_byte;

	if(tree[root].level == 1)
	{ // leaf
		if((tree[root].set_count!=0)&&(tree[root].set_count!=set_entry_array_size))
		{
			leftSets.assign(tree[root].sets,tree[root].sets+tree[root].set_count);
			bool couldCompress=false;
			vector<Dir_entry> left,right;

			for(i = 0; i < set_entry_array_size; i++)
			{
				if(find(leftSets.begin(),leftSets.end(),i)!=leftSets.end())
					for(j = 0; j < set_entry_array[i].freq; j++)
						left.push_back(sort_entries[set_entry_array[i].entry_set[j]]);
				else
					for(j = 0; j < set_entry_array[i].freq; j++)
						right.push_back(sort_entries[set_entry_array[i].entry_set[j]]);
			}

			assert(left.size()+right.size()==count+1);
			if(TREE_TYPE == DYNAMIC_TREE) 
				if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
					couldCompress=true;

			if(couldCompress||(TREE_TYPE == STATIC_TREE))
			{//calculate overlap
				unsigned char right_DMBR[DMBR_SIZE];
				cal_DMBR(right, right.size(),  right_DMBR); // create DMBR for grp1
				rst_overlap = cal_normed_overlap(tree[root].letters, right_DMBR, A);
			}
			else//cannot compress
			{
				leftSets.clear();
				rst_overlap = DOUBLE_INF;
			}
		}
		else //no sets in this node of forest or sets full
		{
			leftSets.clear();
			rst_overlap = DOUBLE_INF;
		}
	}//if is leaf
	else
	{//else not a leaf

		//vector<int> leftSets;
		//double rst_overlap;
		vector<int> tmp_leftSets;
		double tmp_rst_overlap;

		find_minOverlap_by_aux_tree(sort_dim, sort_entries, tree, tree[root].children[0], letters_bitmap_byte_size, set_entry_array, set_entry_array_size, leftSets, rst_overlap);
		
		//assert(leftSets.size()>0);
		assert(leftSets.size()<set_entry_array_size);
		for(int numChildren=1;numChildren<tree[root].cnt;numChildren++)
		{

			find_minOverlap_by_aux_tree(sort_dim, sort_entries, tree, tree[root].children[numChildren], letters_bitmap_byte_size, set_entry_array, set_entry_array_size, tmp_leftSets, tmp_rst_overlap);
			//assert(tmp_leftSets.size()>0);
			assert(tmp_leftSets.size()<set_entry_array_size);
			if(tmp_rst_overlap<rst_overlap)
			{
				leftSets=tmp_leftSets;
				rst_overlap=tmp_rst_overlap;
			}

		}

		if((tree[root].set_count!=0)&&(tree[root].set_count!=set_entry_array_size))
		{
			tmp_leftSets.assign(tree[root].sets,tree[root].sets+tree[root].set_count);
			bool couldCompress=false;
			vector<Dir_entry> left,right;

			for(i = 0; i < set_entry_array_size; i++)
			{
				if(find(tmp_leftSets.begin(),tmp_leftSets.end(),i)!=tmp_leftSets.end())
					for(j = 0; j < set_entry_array[i].freq; j++)
						left.push_back(sort_entries[set_entry_array[i].entry_set[j]]);
				else
					for(j = 0; j < set_entry_array[i].freq; j++)
						right.push_back(sort_entries[set_entry_array[i].entry_set[j]]);

			}

			assert(left.size()+right.size()==count+1);
			if(TREE_TYPE == DYNAMIC_TREE) 
				if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
					couldCompress=true;

			if(couldCompress||(TREE_TYPE == STATIC_TREE))
			{//calculate overlap
				unsigned char right_DMBR[DMBR_SIZE];
				cal_DMBR(right, right.size(),  right_DMBR); // create DMBR for grp1
				tmp_rst_overlap = cal_normed_overlap(tree[root].letters, right_DMBR, A);
				if(tmp_rst_overlap<rst_overlap)
				{
					leftSets=tmp_leftSets;
					rst_overlap=tmp_rst_overlap;
				}
			}
			//else//cannot compress
			//{
			//	tmp_leftSets.clear();
			//	tmp_rst_overlap = DOUBLE_INF;
			//}



		}//end of if((tree[tree[root].children[numChildren]].set_count!=0)&&(tree[tree[root].children[numChildren]].set_count!=set_entry_array_size))
		//else //no sets in this node of forest or sets full
		//{
		//	tmp_leftSets.clear();
		//	tmp_rst_overlap = DOUBLE_INF;
		//}


	}//end of else not a leaf
		return;
}




//derived from sort_entries_by_size_Random_Util_v2
//Knapsack problem
void Dir_node::sort_entries_by_size_Random_Util_v3(
	double dir_min_util,
int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, 
int* sorted_entry_list, 
int & rst_leftNumLetters, 
double & rst_overlap, 
int& rst_cutPoint)
{
   int i, j, k;
   int tmp_sum, tmp_cnt;
   int tmp_index;
   unsigned char tmp_set[MAX_DMBR_DIM_SIZE];
   unsigned char tmp_letters[MAX_DMBR_DIM_SIZE];

   unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
   int size_of_alphabet2 = 0;

   //set_entry set_entry_array[DIR_NODE_SIZE+1];
   set_entry *set_entry_array=new set_entry[count+1];
	for(int t=0;t<count+1;t++) set_entry_array[t].resizeSetEntry(count+1);

   assert(num_of_sort_entries=count+1);

   int num_of_set_entry = 0;

   int cur_dim_byte_size = this->DMBR_end_byte_lut[sort_dim] - this->DMBR_start_byte_lut[sort_dim] + 1; // number of bytes for the bitmap used to represent the set
   // cal histogram of the entries
   // initalize tmp_set
   for(i = 0; i < cur_dim_byte_size; i++)tmp_set[i] = 0;

   // set_entry_array
   bool found;
   for(i = 0; i < num_of_sort_entries; i++)
   {
      // find the position of the component set of the current entry in the set_array
      found = false;
      for(j = 0; j < num_of_set_entry; j++)
         if(set_equal(sort_entries[i].DMBR+this->DMBR_start_byte_lut[sort_dim], set_entry_array[j].set, cur_dim_byte_size)){
            found = true;
            break;
         }
      if(found){
         set_entry_array[j].entry_set[set_entry_array[j].freq] = i;
         set_entry_array[j].freq++;
      }else{ // new set entry found
         tmp_sum = 0;
         for(j = 0; j < cur_dim_byte_size; j++){
            set_entry_array[num_of_set_entry].set[j] = sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
            tmp_sum += bit_1_count_lut[sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j]];
         }
         set_entry_array[num_of_set_entry].element_size = tmp_sum;
         set_entry_array[num_of_set_entry].entry_set[0] = i;
         set_entry_array[num_of_set_entry].freq = 1;

         num_of_set_entry++;
         
         // collect the letters used into the temp set
         for(j = 0; j < cur_dim_byte_size; j++)tmp_set[j] |= sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
      }
   }

   // create alphabet 2 from the tmp set.
   this->bitmap_to_letters(tmp_set, cur_dim_byte_size, alphabet2, size_of_alphabet2);

   // sort sets by element_size asc, break ties by freq desc.  Selection sort
   //int sorted_set_entry_list[DIR_NODE_SIZE+1];  // an array of indices into the set_entry_array
   int *sorted_set_entry_list=new int[count+1];
   
   // populate
   for(i = 0; i < num_of_set_entry; i++)sorted_set_entry_list[i] = i;
   int largest, current, position;
   for(position = num_of_set_entry - 1; position > 0; position--){
      largest = 0;
      for(current = 1; current <= position; current++)
         if((set_entry_array[sorted_set_entry_list[current]].element_size > set_entry_array[sorted_set_entry_list[largest]].element_size) ||
            (set_entry_array[sorted_set_entry_list[current]].element_size == set_entry_array[sorted_set_entry_list[largest]].element_size &&
            set_entry_array[sorted_set_entry_list[current]].freq < set_entry_array[sorted_set_entry_list[largest]].freq))
            largest = current;
      // swap
      tmp_index = sorted_set_entry_list[largest];
      sorted_set_entry_list[largest] = sorted_set_entry_list[position];
      sorted_set_entry_list[position] = tmp_index;
   }


////////////////////debug
//   cout<<"---------------------------\n";
//	  for(int debugp=0;debugp<num_of_set_entry;debugp++)
//	  {
//		  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[debugp]].entry_set[0]],set_entry_array[sorted_set_entry_list[debugp]].entry_set[0],	sort_dim, A);
//
//	  }
//
///////////////////

   if(num_of_set_entry==1)
   {
		for(i=0;i<num_of_sort_entries;i++)
			sorted_entry_list[i]=i;

		rst_leftNumLetters = set_entry_array[0].element_size;
		

		double tmp_rst_overlap;
		rst_overlap= DOUBLE_INF;
		for(int tmp_rst_cutPoint=1;tmp_rst_cutPoint<num_of_sort_entries-1;tmp_rst_cutPoint++)
		{
			unsigned char DMBR1[DMBR_SIZE],DMBR2[DMBR_SIZE];
			cal_DMBR(sort_entries, sorted_entry_list, tmp_rst_cutPoint, DMBR1); // create DMBR for grp1
			cal_DMBR(sort_entries, sorted_entry_list+tmp_rst_cutPoint, num_of_sort_entries- tmp_rst_cutPoint, DMBR2); // create DMBR for grp1
			tmp_rst_overlap = cal_normed_overlap(DMBR1, DMBR2, A);

			if(tmp_rst_overlap<rst_overlap)
			{
				rst_overlap=tmp_rst_overlap;
				rst_cutPoint=tmp_rst_cutPoint;
			}

		}
   }
   else
   {

	   aux_tree_node *forest;


	   forest = new aux_tree_node[MAX_ALPHABET_SIZE + count + 2];
	   for(int f=0;f<(MAX_ALPHABET_SIZE + count + 2);f++) 
		   forest[f].resizeSets(count + 1);



	   int roots[MAX_ALPHABET_SIZE]; // a list of roots in the forest
	   int root_count = size_of_alphabet2;  // initial number of roots in forest
	   int total_num_of_nodes = size_of_alphabet2; // Total number of aux tree nodes in the forest
	   // create forest based on alphabet2
	   for(i = 0; i < size_of_alphabet2; i++){
		   for(j = 0; j < cur_dim_byte_size; j++)forest[i].letters[j] = 0;
		   forest[i].letters[this->bitmap_byte_lut[alphabet2[i]]] |= MASKS[this->bitmap_bit_lut[alphabet2[i]]];
		   forest[i].set_count = 0;
		   forest[i].freq = 0;
		   forest[i].level = 1; // leaf
		   forest[i].cnt = 0;
		   forest[i].parent = -1; // root parent is -1
		   forest[i].numLetters = 1; // root parent is -1

		   roots[i] = i;
	   }

	   int intersect_root_index[MAX_ALPHABET_SIZE]; // The index of an intersect root in roots
	   int intersect_roots[MAX_ALPHABET_SIZE]; // root (index into forest) itself
	   int intersect_root_count;
	   bool is_intersect;

	   int max_level;
	   for(i = 0; i < num_of_set_entry; i++)
	   { // go through all the sets in sorted order of size
		   tmp_index = sorted_set_entry_list[i];
		   intersect_root_count = 0;
		   for(j = 0; j < root_count; j++)
		   { // go through all the tree roots in the forest
			   Node::set_intersect(set_entry_array[tmp_index].set, forest[roots[j]].letters, cur_dim_byte_size, tmp_set, is_intersect);
			   if(is_intersect){
				   intersect_root_index[intersect_root_count] = j;
				   intersect_roots[intersect_root_count] = roots[j];
				   intersect_root_count++;
			   }
		   }
		   if(intersect_root_count == 1){
			   forest[intersect_roots[0]].freq += set_entry_array[tmp_index].freq;
			   forest[intersect_roots[0]].sets[forest[intersect_roots[0]].set_count] = tmp_index;
			   forest[intersect_roots[0]].set_count++;

		   }else if (intersect_root_count > 1)
		   { // tree node merge
			   for(j = 0; j < cur_dim_byte_size; j++)tmp_letters[j] = 0;
			   int sum_of_freq = 0;
			   //bool union_of_sets[DIR_NODE_SIZE+1]; // bit array
			   bool *union_of_sets_ptr=new bool[count+1];
			   for(j = 0; j < num_of_set_entry; j++)
				   union_of_sets_ptr[j] = false;
			   max_level = 1;
			   // merge trees whose root is in intersect_roots
			   for(j = 0; j < intersect_root_count; j++)
			   {
				   forest[intersect_roots[j]].parent = total_num_of_nodes;

				   for(k = 0; k < cur_dim_byte_size; k++)
					   tmp_letters[k] |= forest[intersect_roots[j]].letters[k];
				   sum_of_freq += forest[intersect_roots[j]].freq;
				   for(k = 0; k < forest[intersect_roots[j]].set_count; k++)
					   union_of_sets_ptr[forest[intersect_roots[j]].sets[k]] = true;
				   if(forest[intersect_roots[j]].level > max_level)
					   max_level = forest[intersect_roots[j]].level;
			   }

			   for(j = 0; j < cur_dim_byte_size; j++)
				   forest[total_num_of_nodes].letters[j] = tmp_letters[j] | set_entry_array[tmp_index].set[j];
			   forest[total_num_of_nodes].freq = sum_of_freq + set_entry_array[tmp_index].freq;
			   tmp_cnt = 0;
			   for(j = 0; j < num_of_set_entry; j++)
				   if(union_of_sets_ptr[j]){
					   forest[total_num_of_nodes].sets[tmp_cnt] = j;
					   tmp_cnt++;
				   }
				   forest[total_num_of_nodes].sets[tmp_cnt] = tmp_index;
				   tmp_cnt++;
				   forest[total_num_of_nodes].set_count = tmp_cnt;
				   forest[total_num_of_nodes].level = max_level + 1;
				   for(j = 0; j < intersect_root_count; j++)
					   forest[total_num_of_nodes].children[j] = intersect_roots[j];
				   forest[total_num_of_nodes].cnt = intersect_root_count;
				   forest[total_num_of_nodes].parent = -1;
				   forest[total_num_of_nodes].numLetters=0;
				   for(int s = 0; s < cur_dim_byte_size; s++)
					   forest[total_num_of_nodes].numLetters += bit_1_count_lut[forest[total_num_of_nodes].letters[s]];


				   total_num_of_nodes++;


				   // remove the old roots after merging
				   bool bit_array_is_root[MAX_ALPHABET_SIZE];
				   for(j = 0; j < root_count; j++)bit_array_is_root[j] = true;
				   for(j = 0; j < intersect_root_count; j++)
					   bit_array_is_root[intersect_root_index[j]] = false;
				   tmp_cnt = 0;
				   for(j = 0; j < root_count; j++)
					   if(bit_array_is_root[j]){
						   roots[tmp_cnt] = roots[j];
						   tmp_cnt++;
					   }
					   root_count = tmp_cnt;
					   // add the new merged root
					   roots[root_count] = total_num_of_nodes - 1;
					   root_count++;

					   delete []union_of_sets_ptr;
		   }
	   }
//now forest generation finished

vector<int> result;

	   //logO.log2File("root_count: ");logO.log2File(root_count);logO.log2File("\n");
	   if(root_count > 1)
	   {//exists more than 1 root


		   int * item_size = new int[root_count];
		   int * item_value = new int[root_count];
		   int  all_entries_size=0;
		   for(int i =0;i<root_count;i++)
		   {
			   item_value[i]=forest[roots[i]].numLetters;

			   vector<Dir_entry> tmp_entries;
			   for(int m = 0; m < forest[roots[i]].set_count; m++)
			   {
				   int tmp_set_index = forest[roots[i]].sets[m];
				   for(int n=0;n<set_entry_array[tmp_set_index].freq;n++)
					   tmp_entries.push_back(sort_entries[set_entry_array[tmp_set_index].entry_set[n]]);
			   }

		   if(TREE_TYPE == DYNAMIC_TREE)			   
			   item_size[i]=getCompressedEntriesSize(tmp_entries,A);
		   else
			   item_size[i]=(DMBR_SIZE + sizeof(unsigned int)) * tmp_entries.size();
		   all_entries_size+=item_size[i];
		   }


		   float knapsack_min_capacity;
		   if(TREE_TYPE == STATIC_TREE)
			   knapsack_min_capacity = (DMBR_SIZE + sizeof(unsigned int))*dir_min_util*(count+1);
		   else
			   knapsack_min_capacity = DISK_BLOCK_SIZE*dir_min_util;


		   result =  Knapsack_recursive_forNDTree(root_count, DISK_BLOCK_SIZE - DIR_NODE_OVERHEAD,item_size,item_value,all_entries_size,dir_min_util,DISK_BLOCK_SIZE - DIR_NODE_OVERHEAD);

#ifdef LOG_KANPSACK_VALUE_AND_WEIGHT
		   logO.log2File("splitting_height ");logO.log2File(splitting_height );logO.log2File(" ");

		   if (result.size()!=0)
			   logO.log2File("success! ");
		   else
			   logO.log2File("fail! ");

		   logO.log2File("knapsack weight:value, ");
		   for(int log_i = 0; log_i<root_count;log_i++)
		   {
		   
			   logO.log2File(" ,");logO.log2File(item_size[log_i]);logO.log2File(" : ");logO.log2File(item_value[log_i]);
   
		   }

		   
		   logO.log2File("\n");
	

#endif
			







		   delete []item_size;
		   delete []item_value;
	   }//end of exists more than 1 root

		   if (result.size()==0)
		   {

			   for(i=0;i<num_of_sort_entries;i++)
				   sorted_entry_list[i]=i;
			   rst_overlap= DOUBLE_INF;
			   rst_cutPoint=rst_leftNumLetters = set_entry_array[0].element_size;

			   //for(i=0;i<num_of_sort_entries;i++)
				  // sorted_entry_list[i]=i;

			   //rst_leftNumLetters = set_entry_array[0].element_size;

			   //double tmp_rst_overlap;
			   //rst_overlap= DOUBLE_INF;
			   //for(int tmp_rst_cutPoint=1;tmp_rst_cutPoint<num_of_sort_entries-1;tmp_rst_cutPoint++)
			   //{
				  // unsigned char DMBR1[DMBR_SIZE],DMBR2[DMBR_SIZE];
				  // cal_DMBR(sort_entries, sorted_entry_list, tmp_rst_cutPoint, DMBR1); // create DMBR for grp1
				  // cal_DMBR(sort_entries, sorted_entry_list+tmp_rst_cutPoint, num_of_sort_entries- tmp_rst_cutPoint, DMBR2); // create DMBR for grp1
				  // tmp_rst_overlap = cal_normed_overlap(DMBR1, DMBR2, A);

				  // if(tmp_rst_overlap<rst_overlap)
				  // {
					 //  rst_overlap=tmp_rst_overlap;
					 //  rst_cutPoint=tmp_rst_cutPoint;
				  // }

			   //}


		   }
		   else
		   {//result size not equal 0

			   int leftIndex = 0, right_index = count;
			   rst_leftNumLetters = 0;
			   for(int i =0; i< root_count; i++)
			   {
				   if(find(result.begin(),result.end(),i)!=result.end())
					   for(int m = 0; m < forest[roots[i]].set_count; m++)
					   {
						   int tmp_set_index = forest[roots[i]].sets[m];
						   for(int n=0;n<set_entry_array[tmp_set_index].freq;n++)
						   {
							   sorted_entry_list[right_index]=(set_entry_array[tmp_set_index].entry_set[n]);
							   right_index--;
						   }
					   }
				   else
				   {
					   for(int m = 0; m < forest[roots[i]].set_count; m++)
					   {
						   int tmp_set_index = forest[roots[i]].sets[m];
						   for(int n=0;n<set_entry_array[tmp_set_index].freq;n++)
						   {
							   sorted_entry_list[leftIndex]=(set_entry_array[tmp_set_index].entry_set[n]);
							   leftIndex++;
						   }
					   }

					   rst_leftNumLetters+=forest[roots[i]].numLetters;
				   }



			   }//end of for

////////////////////debug
//cout<<"------------------\n";
//for (int di=0;di<(count+1);di++)
//{
//	cout<<sorted_entry_list[di]<<endl;
//
//}
//////////////////
				assert(leftIndex==(right_index+1));

				unsigned char left_DMBR[DMBR_SIZE],right_DMBR[DMBR_SIZE];
				cal_DMBR(sort_entries, sorted_entry_list, leftIndex, left_DMBR); // create DMBR for grp1
				cal_DMBR(sort_entries, sorted_entry_list+leftIndex, num_of_sort_entries- leftIndex, right_DMBR); // create DMBR for grp1

				rst_overlap = cal_normed_overlap(left_DMBR, right_DMBR, A);

				rst_cutPoint=leftIndex;
		   }


	   ////////////////////debug
	   //   cout<<"---------------------------"<<endl;
	   //	  for(int debugp=0;debugp<leftSets.size();debugp++)
	   //	  {
	   //		  cout<<set_entry_array[sorted_set_entry_list[leftSets.at(debugp)]].freq<<" entries have set :";
	   //		  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[leftSets.at(debugp)]].entry_set[0]],set_entry_array[sorted_set_entry_list[leftSets.at(debugp)]].entry_set[0],	sort_dim, A);
	   //
	   //	  }
	   //
	   //		cout<<"cut here"<<endl;
	   //
	   //	  for(int debugp=0;debugp<num_of_set_entry;debugp++)
	   //	  {
	   //		  if(find(leftSets.begin(),leftSets.end(),debugp)==leftSets.end())
	   //		  {
	   //			  cout<<set_entry_array[sorted_set_entry_list[debugp]].freq<<" entries have set :";
	   //			  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[debugp]].entry_set[0]],set_entry_array[sorted_set_entry_list[debugp]].entry_set[0],	sort_dim, A);
	   //		  }
	   //
	   //	  }
	   //
	   //
	   //
	   ///////////////////

	   delete[] forest;
	   delete[] set_entry_array;		
	   delete[] sorted_set_entry_list;		
   }//end of else

}//end of function


void Dir_node::sort_entries_by_size_Random_Util_v2(int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, 
int* sorted_entry_list, int & rst_leftNumLetters, double & rst_overlap, int& rst_cutPoint)
{
   int i, j, k;
   int tmp_sum, tmp_cnt;
   int tmp_index;
   unsigned char tmp_set[MAX_DMBR_DIM_SIZE];
   unsigned char tmp_letters[MAX_DMBR_DIM_SIZE];

   unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
   int size_of_alphabet2 = 0;

   //set_entry set_entry_array[DIR_NODE_SIZE+1];
   set_entry *set_entry_array=new set_entry[count+1];
	for(int t=0;t<count+1;t++) set_entry_array[t].resizeSetEntry(count+1);

   assert(num_of_sort_entries=count+1);

   int num_of_set_entry = 0;

   int cur_dim_byte_size = this->DMBR_end_byte_lut[sort_dim] - this->DMBR_start_byte_lut[sort_dim] + 1; // number of bytes for the bitmap used to represent the set
   // cal histogram of the entries
   // initalize tmp_set
   for(i = 0; i < cur_dim_byte_size; i++)tmp_set[i] = 0;

   // set_entry_array
   bool found;
   for(i = 0; i < num_of_sort_entries; i++)
   {
      // find the position of the component set of the current entry in the set_array
      found = false;
      for(j = 0; j < num_of_set_entry; j++)
         if(set_equal(sort_entries[i].DMBR+this->DMBR_start_byte_lut[sort_dim], set_entry_array[j].set, cur_dim_byte_size)){
            found = true;
            break;
         }
      if(found){
         set_entry_array[j].entry_set[set_entry_array[j].freq] = i;
         set_entry_array[j].freq++;
      }else{ // new set entry found
         tmp_sum = 0;
         for(j = 0; j < cur_dim_byte_size; j++){
            set_entry_array[num_of_set_entry].set[j] = sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
            tmp_sum += bit_1_count_lut[sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j]];
         }
         set_entry_array[num_of_set_entry].element_size = tmp_sum;
         set_entry_array[num_of_set_entry].entry_set[0] = i;
         set_entry_array[num_of_set_entry].freq = 1;

         num_of_set_entry++;
         
         // collect the letters used into the temp set
         for(j = 0; j < cur_dim_byte_size; j++)tmp_set[j] |= sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
      }
   }

   // create alphabet 2 from the tmp set.
   this->bitmap_to_letters(tmp_set, cur_dim_byte_size, alphabet2, size_of_alphabet2);

   // sort sets by element_size asc, break ties by freq desc.  Selection sort
   //int sorted_set_entry_list[DIR_NODE_SIZE+1];  // an array of indices into the set_entry_array
   int *sorted_set_entry_list=new int[count+1];
   
   // populate
   for(i = 0; i < num_of_set_entry; i++)sorted_set_entry_list[i] = i;
   int largest, current, position;
   for(position = num_of_set_entry - 1; position > 0; position--){
      largest = 0;
      for(current = 1; current <= position; current++)
         if((set_entry_array[sorted_set_entry_list[current]].element_size > set_entry_array[sorted_set_entry_list[largest]].element_size) ||
            (set_entry_array[sorted_set_entry_list[current]].element_size == set_entry_array[sorted_set_entry_list[largest]].element_size &&
            set_entry_array[sorted_set_entry_list[current]].freq < set_entry_array[sorted_set_entry_list[largest]].freq))
            largest = current;
      // swap
      tmp_index = sorted_set_entry_list[largest];
      sorted_set_entry_list[largest] = sorted_set_entry_list[position];
      sorted_set_entry_list[position] = tmp_index;
   }


////////////////////debug
//   cout<<"---------------------------\n";
//	  for(int debugp=0;debugp<num_of_set_entry;debugp++)
//	  {
//		  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[debugp]].entry_set[0]],set_entry_array[sorted_set_entry_list[debugp]].entry_set[0],	sort_dim, A);
//
//	  }
//
///////////////////

   if(num_of_set_entry==1)
   {
		for(i=0;i<num_of_sort_entries;i++)
			sorted_entry_list[i]=i;

		rst_leftNumLetters = set_entry_array[0].element_size;
		
		double tmp_rst_overlap;
		rst_overlap= DOUBLE_INF;
		for(int tmp_rst_cutPoint=1;tmp_rst_cutPoint<num_of_sort_entries-1;tmp_rst_cutPoint++)
		{
			unsigned char DMBR1[DMBR_SIZE],DMBR2[DMBR_SIZE];
			cal_DMBR(sort_entries, sorted_entry_list, tmp_rst_cutPoint, DMBR1); // create DMBR for grp1
			cal_DMBR(sort_entries, sorted_entry_list+tmp_rst_cutPoint, num_of_sort_entries- tmp_rst_cutPoint, DMBR2); // create DMBR for grp1
			tmp_rst_overlap = cal_normed_overlap(DMBR1, DMBR2, A);

			if(tmp_rst_overlap<rst_overlap)
			{
				rst_overlap=tmp_rst_overlap;
				rst_cutPoint=tmp_rst_cutPoint;
			}

		}
   }
   else
   {

	   aux_tree_node *forest;


	   forest = new aux_tree_node[MAX_ALPHABET_SIZE + count + 2];
	   for(int f=0;f<(MAX_ALPHABET_SIZE + count + 2);f++) forest[f].resizeSets(count + 1);



	   int roots[MAX_ALPHABET_SIZE]; // a list of roots in the forest
	   int root_count = size_of_alphabet2;  // initial number of roots in forest
	   int total_num_of_nodes = size_of_alphabet2; // Total number of aux tree nodes in the forest
	   // create forest based on alphabet2
	   for(i = 0; i < size_of_alphabet2; i++){
		   for(j = 0; j < cur_dim_byte_size; j++)forest[i].letters[j] = 0;
		   forest[i].letters[this->bitmap_byte_lut[alphabet2[i]]] |= MASKS[this->bitmap_bit_lut[alphabet2[i]]];
		   forest[i].set_count = 0;
		   forest[i].freq = 0;
		   forest[i].level = 1; // leaf
		   forest[i].cnt = 0;
		   forest[i].parent = -1; // root parent is -1
		   forest[i].numLetters = 1; // root parent is -1

		   roots[i] = i;
	   }

	   int intersect_root_index[MAX_ALPHABET_SIZE]; // The index of an intersect root in roots
	   int intersect_roots[MAX_ALPHABET_SIZE]; // root (index into forest) itself
	   int intersect_root_count;
	   bool is_intersect;

	   int max_level;
	   for(i = 0; i < num_of_set_entry; i++)
	   { // go through all the sets in sorted order of size
		   tmp_index = sorted_set_entry_list[i];
		   intersect_root_count = 0;
		   for(j = 0; j < root_count; j++)
		   { // go through all the tree roots in the forest
			   Node::set_intersect(set_entry_array[tmp_index].set, forest[roots[j]].letters, cur_dim_byte_size, tmp_set, is_intersect);
			   if(is_intersect){
				   intersect_root_index[intersect_root_count] = j;
				   intersect_roots[intersect_root_count] = roots[j];
				   intersect_root_count++;
			   }
		   }
		   if(intersect_root_count == 1){
			   forest[intersect_roots[0]].freq += set_entry_array[tmp_index].freq;
			   forest[intersect_roots[0]].sets[forest[intersect_roots[0]].set_count] = tmp_index;
			   forest[intersect_roots[0]].set_count++;

		   }else if (intersect_root_count > 1)
		   { // tree node merge
			   for(j = 0; j < cur_dim_byte_size; j++)tmp_letters[j] = 0;
			   int sum_of_freq = 0;
			   //bool union_of_sets[DIR_NODE_SIZE+1]; // bit array
			   bool *union_of_sets_ptr=new bool[count+1];
			   for(j = 0; j < num_of_set_entry; j++)
				   union_of_sets_ptr[j] = false;
			   max_level = 1;
			   // merge trees whose root is in intersect_roots
			   for(j = 0; j < intersect_root_count; j++)
			   {
				   forest[intersect_roots[j]].parent = total_num_of_nodes;

				   for(k = 0; k < cur_dim_byte_size; k++)
					   tmp_letters[k] |= forest[intersect_roots[j]].letters[k];
				   sum_of_freq += forest[intersect_roots[j]].freq;
				   for(k = 0; k < forest[intersect_roots[j]].set_count; k++)
					   union_of_sets_ptr[forest[intersect_roots[j]].sets[k]] = true;
				   if(forest[intersect_roots[j]].level > max_level)
					   max_level = forest[intersect_roots[j]].level;
			   }

			   for(j = 0; j < cur_dim_byte_size; j++)
				   forest[total_num_of_nodes].letters[j] = tmp_letters[j] | set_entry_array[tmp_index].set[j];
			   forest[total_num_of_nodes].freq = sum_of_freq + set_entry_array[tmp_index].freq;
			   tmp_cnt = 0;
			   for(j = 0; j < num_of_set_entry; j++)
				   if(union_of_sets_ptr[j]){
					   forest[total_num_of_nodes].sets[tmp_cnt] = j;
					   tmp_cnt++;
				   }
				   forest[total_num_of_nodes].sets[tmp_cnt] = tmp_index;
				   tmp_cnt++;
				   forest[total_num_of_nodes].set_count = tmp_cnt;
				   forest[total_num_of_nodes].level = max_level + 1;
				   for(j = 0; j < intersect_root_count; j++)
					   forest[total_num_of_nodes].children[j] = intersect_roots[j];
				   forest[total_num_of_nodes].cnt = intersect_root_count;
				   forest[total_num_of_nodes].parent = -1;
				   forest[total_num_of_nodes].numLetters=0;
				   for(int s = 0; s < cur_dim_byte_size; s++)
					   forest[total_num_of_nodes].numLetters += bit_1_count_lut[forest[total_num_of_nodes].letters[s]];


				   total_num_of_nodes++;


				   // remove the old roots after merging
				   bool bit_array_is_root[MAX_ALPHABET_SIZE];
				   for(j = 0; j < root_count; j++)bit_array_is_root[j] = true;
				   for(j = 0; j < intersect_root_count; j++)
					   bit_array_is_root[intersect_root_index[j]] = false;
				   tmp_cnt = 0;
				   for(j = 0; j < root_count; j++)
					   if(bit_array_is_root[j]){
						   roots[tmp_cnt] = roots[j];
						   tmp_cnt++;
					   }
					   root_count = tmp_cnt;
					   // add the new merged root
					   roots[root_count] = total_num_of_nodes - 1;
					   root_count++;

					   delete []union_of_sets_ptr;
		   }
	   }


	   bool foundBest=false;
	   vector<Dir_entry> left, right;
	   vector<int> leftSets;

	   double numLettersArray[MAX_ALPHABET_SIZE];
	   int sortedRootsByNumLetters[MAX_ALPHABET_SIZE];
	   if(root_count > 1)
	   {//exists more than 1 root
		   for(i = 0; i < root_count; i++)
			   numLettersArray[i]=forest[roots[i]].numLetters;

		   //sort roots by the number of letters on this dimension
		   NDT_stable_sort(root_count, numLettersArray, sortedRootsByNumLetters, true);

		   if(TREE_TYPE == STATIC_TREE) 
		   {//NO compression needed, so find the right one
			   foundBest=true;
			   for(int t=0;t<forest[roots[sortedRootsByNumLetters[0]]].set_count;t++)
				   leftSets.push_back(forest[roots[sortedRootsByNumLetters[0]]].sets[t]);
			   //leftIndex=sortedRootsByNumLetters[0];
			   rst_leftNumLetters=(int)numLettersArray[0];
			   rst_overlap=0.0;
			   rst_cutPoint=forest[roots[sortedRootsByNumLetters[0]]].freq;
		   }
		   else
		   {// is DYNAMIC_TREE
			   for(i=0;i<root_count;i++)
			   {// for i
				   left.clear();
				   right.clear();
				   int left_root_Index=sortedRootsByNumLetters[i];
				   for(int m = 0; m < forest[roots[left_root_Index]].set_count; m++)
				   {
					   int tmp_set_index = forest[roots[left_root_Index]].sets[m];
					   for(int n=0;n<set_entry_array[tmp_set_index].freq;n++)
						   left.push_back(sort_entries[set_entry_array[tmp_set_index].entry_set[n]]);

				   }
				   for(j=0;j<root_count;j++)
				   { // for j
					   if(j!=left_root_Index)
					   {
						   for(int m = 0; m < forest[roots[j]].set_count; m++)
						   {
							   int tmp_set_index = forest[roots[j]].sets[m];
							   for(int n=0;n<set_entry_array[tmp_set_index].freq;n++)
								   right.push_back(sort_entries[set_entry_array[tmp_set_index].entry_set[n]]);

						   }
					   }

				   }//end of for j

				   assert((left.size()+right.size())==num_of_sort_entries);
				   assert(forest[roots[left_root_Index]].numLetters==numLettersArray[i]);
				   if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
				   {//compression ok, so find the right one
					   foundBest=true;
					   for(int t=0;t<forest[roots[left_root_Index]].set_count;t++)
						   leftSets.push_back(forest[roots[left_root_Index]].sets[t]);
					   rst_leftNumLetters=(int)numLettersArray[i];
					   rst_overlap=0.0;
					   rst_cutPoint=forest[roots[left_root_Index]].freq;
					   break;
				   }

			   }//end of for i
		   }//end of is DYNAMIC_TREE
	   }

int debugForestNum=root_count;

	   if (!foundBest)
	   {// can not find.
		   //either because only 1 root or more root but cannot compress

		   // merge if more than 1 tree left in the forest 
		   if(root_count > 1){
			   // merge trees
			   for(i = 0; i < cur_dim_byte_size; i++)tmp_letters[i] = 0;
			   int sum_of_freq = 0;
			   //bool union_of_sets[DIR_NODE_SIZE+1]; // bit array
			   bool *union_of_sets_ptr=new bool[count+1]; // bit array

			   for(i = 0; i < num_of_set_entry; i++)union_of_sets_ptr[i] = false;
			   max_level = 1;
			   // merge trees whose root is in roots
			   for(i = 0; i < root_count; i++){
				   forest[roots[i]].parent = total_num_of_nodes;
				   for(j = 0; j < cur_dim_byte_size; j++)
					   tmp_letters[j] |= forest[roots[i]].letters[j];
				   sum_of_freq += forest[roots[i]].freq;
				   for(j = 0; j < forest[roots[i]].set_count; j++)
					   union_of_sets_ptr[forest[roots[i]].sets[j]] = true;
				   if(forest[roots[i]].level > max_level)
					   max_level = forest[roots[i]].level;
			   }
			   for(i = 0; i < cur_dim_byte_size; i++)
				   forest[total_num_of_nodes].letters[i] = tmp_letters[i];
			   forest[total_num_of_nodes].freq = sum_of_freq;
			   tmp_cnt = 0;
			   for(i = 0; i < num_of_set_entry; i++)
				   if(union_of_sets_ptr[i]){
					   forest[total_num_of_nodes].sets[tmp_cnt] = i;
					   tmp_cnt++;
				   }
				   forest[total_num_of_nodes].set_count = tmp_cnt;
				   forest[total_num_of_nodes].level = max_level + 1;
				   for(i = 0; i < root_count; i++)
					   forest[total_num_of_nodes].children[i] = roots[i];
				   forest[total_num_of_nodes].cnt = root_count;
				   forest[total_num_of_nodes].parent = -1;
				   total_num_of_nodes++;

				   // remove the old merged roots
				   root_count = 1;
				   roots[0] = total_num_of_nodes - 1;
				   delete []union_of_sets_ptr;	   
		   }
		   assert(forest[roots[0]].freq == count + 1);
		   assert(forest[roots[0]].set_count == num_of_set_entry);

		   find_minOverlap_by_aux_tree(sort_dim, sort_entries,forest, roots[0], cur_dim_byte_size, set_entry_array, num_of_set_entry, leftSets, rst_overlap);

		   //when is static tree and forest num >1, overlap should be 0
			assert((TREE_TYPE == DYNAMIC_TREE)||((debugForestNum==1)||(rst_overlap==0)));
		   //assert((TREE_TYPE == DYNAMIC_TREE)||(((debugForestNum==1)&&(rst_overlap>0))||((debugForestNum>1)&&(rst_overlap==0))));//because no overlap free partition could be found
		   //assert((debugForestNum==1)||(rst_overlap>0));//because no overlap free partition could be found
		   //assert((TREE_TYPE == STATIC_TREE)||((debugForestNum==1)&&(rst_overlap>0)));//because no overlap free partition could be found


		   //assert(leftSets.size()>0);
		   //assert(leftSets.size()<num_of_set_entry);

		   if((leftSets.size()==0)||(leftSets.size()==num_of_set_entry))
		   {//resort to decide_order_by_aux_tree
			   int list_size = 0;
			   decide_order_by_aux_tree(forest, roots[0], cur_dim_byte_size, set_entry_array, num_of_set_entry, sorted_set_entry_list, list_size);
			   assert(list_size == num_of_set_entry);

			   int tmp_cnt = 0;
			   for(i = 0; i < num_of_set_entry; i++)
			   {
				   for(j = 0; j < set_entry_array[sorted_set_entry_list[i]].freq; j++)
				   {
					   sorted_entry_list[tmp_cnt] = set_entry_array[sorted_set_entry_list[i]].entry_set[j];
					   tmp_cnt++;
				   }
			   }
			   assert(tmp_cnt == num_of_sort_entries); 


			   rst_overlap=DOUBLE_INF;
			   for(int tmp_cutPoint=0;tmp_cutPoint<num_of_sort_entries-1;tmp_cutPoint++)
			   { 
				   vector<Dir_entry> left,right;

				   for(int x=0;x<num_of_sort_entries;x++)
				   {
					   if(x<=tmp_cutPoint)
						   left.push_back(sort_entries[sorted_entry_list[x]]);
					   else
						   right.push_back(sort_entries[sorted_entry_list[x]]);
				   }
				   assert((left.size()+right.size())==(count+1));


				   bool couldCompress=false;
				   if(TREE_TYPE == DYNAMIC_TREE) 
					   if((getCompressedDiskSize(left,A)<= DISK_BLOCK_SIZE)&&(getCompressedDiskSize(right,A)<= DISK_BLOCK_SIZE))
						   couldCompress=true;


				   if(couldCompress||(TREE_TYPE == STATIC_TREE))
				   {//calculate overlap
					   unsigned char left_DMBR[DMBR_SIZE];
					   unsigned char right_DMBR[DMBR_SIZE];
					   cal_DMBR(left, left.size(),  left_DMBR); // create DMBR for grp1
					   cal_DMBR(right, right.size(),  right_DMBR); // create DMBR for grp1

					   double tmp_overlap;		 
					   tmp_overlap = cal_normed_overlap(left_DMBR, right_DMBR, A);
					   if(tmp_overlap<rst_overlap)
					   {
						   assert((tmp_cutPoint+1)==left.size());
						   rst_overlap=tmp_overlap;
						   rst_cutPoint=left.size();
						   rst_leftNumLetters=A[0]*10;//meaningless since when overlap is not 0, only check overlap
					   }
				   }

			   }//end of for(int tmp_cutPoint=0;tmp_cutPoint<num_of_sort_entries-1;tmp_cutPoint++)

			   assert(rst_cutPoint>0);
			   assert(rst_cutPoint<num_of_sort_entries);

			   delete[] forest;
			   delete[] set_entry_array;		
			   delete[] sorted_set_entry_list;		

			   return;
		   }//end of //resort to decide_order_by_aux_tree
		   else
		   {
			   for(int t=0;t<leftSets.size();t++)
			   {
				   rst_cutPoint+=set_entry_array[leftSets.at(t)].freq;
				   rst_leftNumLetters=A[0]*10;//meaningless since when overlap is not 0, only check overlap

			   }
		   }


	   }

	   //int list_size = 0;
	   //decide_order_by_aux_tree(forest, roots[0], cur_dim_byte_size, set_entry_array, num_of_set_entry, sorted_set_entry_list, list_size);
	   //assert(list_size == num_of_set_entry);

	   int tmp_cnt1 = 0,tmp_cnt2=num_of_sort_entries-1;
	   for(i = 0; i < num_of_set_entry; i++)
	   {
		   if(find(leftSets.begin(),leftSets.end(),i)!=leftSets.end())
			   for(j = 0; j < set_entry_array[i].freq; j++)
			   {
				   sorted_entry_list[tmp_cnt1] = set_entry_array[i].entry_set[j];
				   tmp_cnt1++;
			   }
		   else
			   for(j = 0; j < set_entry_array[i].freq; j++)
			   {
				   sorted_entry_list[tmp_cnt2] = set_entry_array[i].entry_set[j];
				   tmp_cnt2--;
			   }

	   }

	   rst_cutPoint=tmp_cnt1;
	   assert(tmp_cnt1 == (tmp_cnt2+1));

	   assert(rst_cutPoint>0);
	   assert(rst_cutPoint<num_of_sort_entries);


	   ////////////////////debug
	   //   cout<<"---------------------------"<<endl;
	   //	  for(int debugp=0;debugp<leftSets.size();debugp++)
	   //	  {
	   //		  cout<<set_entry_array[sorted_set_entry_list[leftSets.at(debugp)]].freq<<" entries have set :";
	   //		  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[leftSets.at(debugp)]].entry_set[0]],set_entry_array[sorted_set_entry_list[leftSets.at(debugp)]].entry_set[0],	sort_dim, A);
	   //
	   //	  }
	   //
	   //		cout<<"cut here"<<endl;
	   //
	   //	  for(int debugp=0;debugp<num_of_set_entry;debugp++)
	   //	  {
	   //		  if(find(leftSets.begin(),leftSets.end(),debugp)==leftSets.end())
	   //		  {
	   //			  cout<<set_entry_array[sorted_set_entry_list[debugp]].freq<<" entries have set :";
	   //			  print_OneEntry_OnOneDscDim(&sort_entries[set_entry_array[sorted_set_entry_list[debugp]].entry_set[0]],set_entry_array[sorted_set_entry_list[debugp]].entry_set[0],	sort_dim, A);
	   //		  }
	   //
	   //	  }
	   //
	   //
	   //
	   ///////////////////


	   delete[] forest;
	   delete[] set_entry_array;		
	   delete[] sorted_set_entry_list;		
   }//end of else

}//end of function







void Dir_node::sort_entries_by_size(int sort_dim, Dir_entry* sort_entries, int num_of_sort_entries, int* sorted_entry_list){
   int i, j, k;
   int tmp_sum, tmp_cnt;
   int tmp_index;
   unsigned char tmp_set[MAX_DMBR_DIM_SIZE];
   unsigned char tmp_letters[MAX_DMBR_DIM_SIZE];

   unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
   int size_of_alphabet2 = 0;

   //set_entry set_entry_array[DIR_NODE_SIZE+1];
   set_entry *set_entry_array=new set_entry[count+1];
	for(int t=0;t<count+1;t++) set_entry_array[t].resizeSetEntry(count+1);

   
   int num_of_set_entry = 0;

   int cur_dim_byte_size = this->DMBR_end_byte_lut[sort_dim] - this->DMBR_start_byte_lut[sort_dim] + 1; // number of bytes for the bitmap used to represent the set
   // cal histogram of the entries
   // initalize tmp_set
   for(i = 0; i < cur_dim_byte_size; i++)tmp_set[i] = 0;

   // set_entry_array
   bool found;
   for(i = 0; i < num_of_sort_entries; i++){
      // find the position of the component set of the current entry in the set_array
      found = false;
      for(j = 0; j < num_of_set_entry; j++)
         if(set_equal(sort_entries[i].DMBR+this->DMBR_start_byte_lut[sort_dim], set_entry_array[j].set, cur_dim_byte_size)){
            found = true;
            break;
         }
      if(found){
         set_entry_array[j].entry_set[set_entry_array[j].freq] = i;
         set_entry_array[j].freq++;
      }else{ // new set entry found
         tmp_sum = 0;
         for(j = 0; j < cur_dim_byte_size; j++){
            set_entry_array[num_of_set_entry].set[j] = sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
            tmp_sum += bit_1_count_lut[sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j]];
         }
         set_entry_array[num_of_set_entry].element_size = tmp_sum;
         set_entry_array[num_of_set_entry].entry_set[0] = i;
         set_entry_array[num_of_set_entry].freq = 1;

         num_of_set_entry++;
         
         // collect the letters used into the temp set
         for(j = 0; j < cur_dim_byte_size; j++)tmp_set[j] |= sort_entries[i].DMBR[this->DMBR_start_byte_lut[sort_dim]+j];
      }
   }
   // create alphabet 2 from the tmp set.
   this->bitmap_to_letters(tmp_set, cur_dim_byte_size, alphabet2, size_of_alphabet2);

   // sort sets by element_size asc, break ties by freq desc.  Selection sort
   //int sorted_set_entry_list[DIR_NODE_SIZE+1];  // an array of indices into the set_entry_array
   int *sorted_set_entry_list=new int[count+1];
   
   // populate
   for(i = 0; i < num_of_set_entry; i++)sorted_set_entry_list[i] = i;
   int largest, current, position;
   for(position = num_of_set_entry - 1; position > 0; position--){
      largest = 0;
      for(current = 1; current <= position; current++)
         if((set_entry_array[sorted_set_entry_list[current]].element_size > set_entry_array[sorted_set_entry_list[largest]].element_size) ||
            (set_entry_array[sorted_set_entry_list[current]].element_size == set_entry_array[sorted_set_entry_list[largest]].element_size &&
            set_entry_array[sorted_set_entry_list[current]].freq < set_entry_array[sorted_set_entry_list[largest]].freq))
            largest = current;
      // swap
      tmp_index = sorted_set_entry_list[largest];
      sorted_set_entry_list[largest] = sorted_set_entry_list[position];
      sorted_set_entry_list[position] = tmp_index;
   }

aux_tree_node *forest;


forest = new aux_tree_node[MAX_ALPHABET_SIZE + count + 2];
for(int f=0;f<(MAX_ALPHABET_SIZE + count + 2);f++) forest[f].resizeSets(count + 1);



   int roots[MAX_ALPHABET_SIZE]; // a list of roots in the forest
   int root_count = size_of_alphabet2;  // initial number of roots in forest
   int total_num_of_nodes = size_of_alphabet2; // Total number of aux tree nodes in the forest
   // create forest based on alphabet2
   for(i = 0; i < size_of_alphabet2; i++){
      for(j = 0; j < cur_dim_byte_size; j++)forest[i].letters[j] = 0;
      forest[i].letters[this->bitmap_byte_lut[alphabet2[i]]] |= MASKS[this->bitmap_bit_lut[alphabet2[i]]];
      forest[i].set_count = 0;
      forest[i].freq = 0;
      forest[i].level = 1; // leaf
      forest[i].cnt = 0;
      forest[i].parent = -1; // root parent is -1

      roots[i] = i;
   }

   int intersect_root_index[MAX_ALPHABET_SIZE]; // The index of an intersect root in roots
   int intersect_roots[MAX_ALPHABET_SIZE]; // root (index into forest) itself
   int intersect_root_count;
   bool is_intersect;

   int max_level;
   for(i = 0; i < num_of_set_entry; i++)
   { // go through all the sets in sorted order of size
      tmp_index = sorted_set_entry_list[i];
      intersect_root_count = 0;
      for(j = 0; j < root_count; j++)
	  { // go through all the tree roots in the forest
         Node::set_intersect(set_entry_array[tmp_index].set, forest[roots[j]].letters, cur_dim_byte_size, tmp_set, is_intersect);
         if(is_intersect){
            intersect_root_index[intersect_root_count] = j;
            intersect_roots[intersect_root_count] = roots[j];
            intersect_root_count++;
         }
      }
      if(intersect_root_count == 1){
         forest[intersect_roots[0]].freq += set_entry_array[tmp_index].freq;
         forest[intersect_roots[0]].sets[forest[intersect_roots[0]].set_count] = tmp_index;
         forest[intersect_roots[0]].set_count++;
      }else if (intersect_root_count > 1)
	  { // tree node merge
         for(j = 0; j < cur_dim_byte_size; j++)tmp_letters[j] = 0;
         int sum_of_freq = 0;
         //bool union_of_sets[DIR_NODE_SIZE+1]; // bit array
		 bool *union_of_sets_ptr=new bool[count+1];
         for(j = 0; j < num_of_set_entry; j++)
			union_of_sets_ptr[j] = false;
         max_level = 1;
         // merge trees whose root is in intersect_roots
         for(j = 0; j < intersect_root_count; j++)
		 {
            forest[intersect_roots[j]].parent = total_num_of_nodes;
            
            for(k = 0; k < cur_dim_byte_size; k++)
               tmp_letters[k] |= forest[intersect_roots[j]].letters[k];
            sum_of_freq += forest[intersect_roots[j]].freq;
            for(k = 0; k < forest[intersect_roots[j]].set_count; k++)
               union_of_sets_ptr[forest[intersect_roots[j]].sets[k]] = true;
            if(forest[intersect_roots[j]].level > max_level)
               max_level = forest[intersect_roots[j]].level;
         }

         for(j = 0; j < cur_dim_byte_size; j++)
            forest[total_num_of_nodes].letters[j] = tmp_letters[j] | set_entry_array[tmp_index].set[j];
         forest[total_num_of_nodes].freq = sum_of_freq + set_entry_array[tmp_index].freq;
         tmp_cnt = 0;
         for(j = 0; j < num_of_set_entry; j++)
            if(union_of_sets_ptr[j]){
               forest[total_num_of_nodes].sets[tmp_cnt] = j;
               tmp_cnt++;
            }
         forest[total_num_of_nodes].sets[tmp_cnt] = tmp_index;
         tmp_cnt++;
         forest[total_num_of_nodes].set_count = tmp_cnt;
         forest[total_num_of_nodes].level = max_level + 1;
         for(j = 0; j < intersect_root_count; j++)
            forest[total_num_of_nodes].children[j] = intersect_roots[j];
         forest[total_num_of_nodes].cnt = intersect_root_count;
         forest[total_num_of_nodes].parent = -1;
         total_num_of_nodes++;
        
         // remove the old roots after merging
         bool bit_array_is_root[MAX_ALPHABET_SIZE];
         for(j = 0; j < root_count; j++)bit_array_is_root[j] = true;
         for(j = 0; j < intersect_root_count; j++)
            bit_array_is_root[intersect_root_index[j]] = false;
         tmp_cnt = 0;
         for(j = 0; j < root_count; j++)
            if(bit_array_is_root[j]){
               roots[tmp_cnt] = roots[j];
               tmp_cnt++;
            }
         root_count = tmp_cnt;
         // add the new merged root
         roots[root_count] = total_num_of_nodes - 1;
         root_count++;
		 
		 delete []union_of_sets_ptr;
      }
   }

   // merge if more than 1 tree left in the forest 
   if(root_count > 1){
      // merge trees
      for(i = 0; i < cur_dim_byte_size; i++)tmp_letters[i] = 0;
	   int sum_of_freq = 0;
      //bool union_of_sets[DIR_NODE_SIZE+1]; // bit array
	  bool *union_of_sets_ptr=new bool[count+1]; // bit array
  
      for(i = 0; i < num_of_set_entry; i++)union_of_sets_ptr[i] = false;
	   max_level = 1;
      // merge trees whose root is in roots
      for(i = 0; i < root_count; i++){
         forest[roots[i]].parent = total_num_of_nodes;
         for(j = 0; j < cur_dim_byte_size; j++)
            tmp_letters[j] |= forest[roots[i]].letters[j];
         sum_of_freq += forest[roots[i]].freq;
         for(j = 0; j < forest[roots[i]].set_count; j++)
            union_of_sets_ptr[forest[roots[i]].sets[j]] = true;
         if(forest[roots[i]].level > max_level)
            max_level = forest[roots[i]].level;
      }
      for(i = 0; i < cur_dim_byte_size; i++)
         forest[total_num_of_nodes].letters[i] = tmp_letters[i];
      forest[total_num_of_nodes].freq = sum_of_freq;
      tmp_cnt = 0;
      for(i = 0; i < num_of_set_entry; i++)
         if(union_of_sets_ptr[i]){
            forest[total_num_of_nodes].sets[tmp_cnt] = i;
            tmp_cnt++;
         }
      forest[total_num_of_nodes].set_count = tmp_cnt;
	   forest[total_num_of_nodes].level = max_level + 1;
      for(i = 0; i < root_count; i++)
         forest[total_num_of_nodes].children[i] = roots[i];
      forest[total_num_of_nodes].cnt = root_count;
      forest[total_num_of_nodes].parent = -1;
      total_num_of_nodes++;
	
   	// remove the old merged roots
      root_count = 1;
	   roots[0] = total_num_of_nodes - 1;
	  delete []union_of_sets_ptr;	   
   }
   assert(forest[roots[0]].freq == count + 1);
   assert(forest[roots[0]].set_count == num_of_set_entry);

   int list_size = 0;
   decide_order_by_aux_tree(forest, roots[0], cur_dim_byte_size, set_entry_array, num_of_set_entry, sorted_set_entry_list, list_size);
   assert(list_size == num_of_set_entry);

   tmp_cnt = 0;
   for(i = 0; i < num_of_set_entry; i++)
      for(j = 0; j < set_entry_array[sorted_set_entry_list[i]].freq; j++){
         sorted_entry_list[tmp_cnt] = set_entry_array[sorted_set_entry_list[i]].entry_set[j];
         tmp_cnt++;
      }
   assert(tmp_cnt == num_of_sort_entries);


		delete[] forest;
		delete[] set_entry_array;		
		delete[] sorted_set_entry_list;		
}



void Dir_node::decide_order_by_aux_tree(aux_tree_node* tree, int root, int letters_bitmap_byte_size, set_entry* set_entry_array, int set_entry_array_size, int* sorted_set_list, int& list_size){
   int i, j, k, l, m;
   int tmp_index;
   int tmp_start_byte, tmp_end_byte;

   if(tree[root].level == 1){ // leaf
      sorted_set_list[list_size] = tree[root].sets[0]; // guranteed to be one set only at leaf
      list_size++;
      return;
   }else{
      //mixed_list_entry mixed_list[DIR_NODE_SIZE+1];
      mixed_list_entry *mixed_list=new mixed_list_entry[count+1];	  
      int mixed_list_size = 0;

      if(tree[root].cnt == 1){ // only one sub tree
         if(tree[tree[root].children[0]].freq > 0){ // only count non-empty trees
            mixed_list[mixed_list_size].index = tree[root].children[0];
            mixed_list[mixed_list_size].is_set = false; // it is a root
            mixed_list_size++;
         }
      }else{ // several subtree
         // create the list to be sorted
         for(i = 0; i < tree[root].cnt; i++){
            tmp_index = tree[root].children[i];
            if(tree[tmp_index].freq > 0){ // only subtree with corresponding entries count
               mixed_list[mixed_list_size].index = tmp_index;
               mixed_list[mixed_list_size].is_set = false; // root
               mixed_list[mixed_list_size].weight = tree[tmp_index].freq;
               mixed_list_size++;
            }
         }
         // decide order of subtree based on histogram freq
         int unsorted_entry_count = mixed_list_size;
         int weight1 = 0, weight2 = 0;
         int left = 0, right = mixed_list_size - 1; // indices
         int largest, current;
         mixed_list_entry tmp_list_entry;
         while(unsorted_entry_count > 0){
            if(weight1 <= weight2){ // next child should be put to the left side
               largest = left;
               for(current = left + 1; current <= right; current++) // pick the letter with the largest freq
                  if(mixed_list[current].weight > mixed_list[largest].weight)
                     largest = current;
               // swap
               tmp_list_entry = mixed_list[largest];
               mixed_list[largest] = mixed_list[left];
               mixed_list[left] = tmp_list_entry;

               weight1 += mixed_list[left].weight;
               left++;
            }else{
               largest = right;
               for(current = left; current <= right - 1; current++) // pick the letter with the largest freq
                  if(mixed_list[current].weight > mixed_list[largest].weight)
                     largest = current;
               // swap
               tmp_list_entry = mixed_list[largest];
               mixed_list[largest] = mixed_list[right];
               mixed_list[right] = tmp_list_entry;

               weight2 += mixed_list[right].weight;
               right--;
            }
            unsorted_entry_count--;
         }
      }

      // insert crossing sets
      // size of letter_union_arrays is mixed_list_size + 1
      //unsigned char letter_union_array_left[MAX_DMBR_DIM_SIZE * DIR_NODE_SIZE];
      //unsigned char letter_union_array_right[MAX_DMBR_DIM_SIZE * DIR_NODE_SIZE];
      unsigned char *letter_union_array_left=new unsigned char[MAX_DMBR_DIM_SIZE * count];
      unsigned char *letter_union_array_right=new unsigned char[MAX_DMBR_DIM_SIZE * count];	  
	  
      // calculate union of letters at each interval
      // left
      // initialize leftmost item
      for(i = 0; i < letters_bitmap_byte_size; i++)
         letter_union_array_left[i] = 0;
      for(i = 0; i < mixed_list_size; i++){
         tmp_start_byte = (i + 1) * letters_bitmap_byte_size;
         tmp_end_byte = (i + 2) * letters_bitmap_byte_size - 1;
         for(j = tmp_start_byte, k = 0; j <= tmp_end_byte; j++, k++)
            letter_union_array_left[j] = letter_union_array_left[j - letters_bitmap_byte_size] | tree[mixed_list[i].index].letters[k];
      }
      // right
      tmp_start_byte = mixed_list_size * letters_bitmap_byte_size;
      tmp_end_byte = (mixed_list_size + 1) * letters_bitmap_byte_size - 1;
      // initialize rightmost item
      for(i = tmp_start_byte; i <= tmp_end_byte; i++)
         letter_union_array_right[i] = 0;
      for(i = mixed_list_size - 1; i >= 0; i--){
         tmp_start_byte = i * letters_bitmap_byte_size;
         tmp_end_byte = (i + 1) * letters_bitmap_byte_size - 1;
         for(j = tmp_start_byte, k = 0; j <= tmp_end_byte; j++, k++)
            letter_union_array_right[j] = letter_union_array_right[j + letters_bitmap_byte_size] | tree[mixed_list[i].index].letters[k];
      }
    
      // get crossing sets.  All sets in parent minus sets in children
      //int crossing_sets[DIR_NODE_SIZE+1]; // set indices
      int *crossing_sets=new int[count+1]; // set indices	  
	  
      int crossing_set_count;

      //bool bit_array_sets[DIR_NODE_SIZE+1];
      bool *bit_array_sets=new bool[count+1];	  
      for(i = 0; i < set_entry_array_size; i++)bit_array_sets[i] = false;
      for(i = 0; i < tree[root].set_count; i++)bit_array_sets[tree[root].sets[i]] = true;
      for(i = 0; i < tree[root].cnt; i++){
         tmp_index = tree[root].children[i];
         for(j = 0; j < tree[tmp_index].set_count; j++)
            bit_array_sets[tree[tmp_index].sets[j]] = false;
      }
      crossing_set_count = 0;
      for(i = 0; i < set_entry_array_size; i++)
         if(bit_array_sets[i]){
            crossing_sets[crossing_set_count] = i;
            crossing_set_count++;
         }

      // insert every crossing set
//      unsigned char tmp_letter_union_array_left[MAX_DMBR_DIM_SIZE * DIR_NODE_SIZE];
//      unsigned char tmp_letter_union_array_right[MAX_DMBR_DIM_SIZE * DIR_NODE_SIZE];
//      unsigned char best_letter_union_array_left[MAX_DMBR_DIM_SIZE * DIR_NODE_SIZE];
//      unsigned char best_letter_union_array_right[MAX_DMBR_DIM_SIZE * DIR_NODE_SIZE];
	  
      unsigned char *tmp_letter_union_array_left=new unsigned char[MAX_DMBR_DIM_SIZE * count];
      unsigned char *tmp_letter_union_array_right=new unsigned char[MAX_DMBR_DIM_SIZE * count];
      unsigned char *best_letter_union_array_left=new unsigned char[MAX_DMBR_DIM_SIZE * count];
      unsigned char *best_letter_union_array_right=new unsigned char[MAX_DMBR_DIM_SIZE * count];	  
	  
	  
      int best_sum_overlap, sum_overlap;
      int best_pos;
      unsigned char tmp_set[MAX_DMBR_DIM_SIZE];
      bool is_intersect;
      mixed_list_entry tmp_list_entry;
      for(i = 0; i < crossing_set_count; i++){ // insert every crossing set into the mixed list
         // decide the best interval to put the crossing set in
         if(mixed_list_size == 0){
            for(j = 0; j < letters_bitmap_byte_size; j++)
               letter_union_array_left[j] = 0;
            tmp_start_byte = letters_bitmap_byte_size;
            tmp_end_byte = 2 * letters_bitmap_byte_size - 1;
            for(j = tmp_start_byte, k = 0; j <= tmp_end_byte; j++, k++)
               letter_union_array_left[j] = set_entry_array[crossing_sets[i]].set[k];

            for(j = 0; j < letters_bitmap_byte_size; j++)
               letter_union_array_right[j] = set_entry_array[crossing_sets[i]].set[j];
            tmp_start_byte = letters_bitmap_byte_size;
            tmp_end_byte = 2 * letters_bitmap_byte_size - 1;
            for(j = tmp_start_byte; j <= tmp_end_byte; j++)
               letter_union_array_right[j] = 0;
            
            // put current set into mixed_list
            mixed_list[0].index = crossing_sets[i];
            mixed_list[0].is_set = true;
            mixed_list_size = 1;
         }else{
            // same for all locations
            for(j = 0; j < letters_bitmap_byte_size; j++)
               tmp_letter_union_array_left[j] = 0;
            tmp_start_byte = (mixed_list_size + 1) * letters_bitmap_byte_size;
            tmp_end_byte = (mixed_list_size + 2) * letters_bitmap_byte_size - 1;
            for(j = tmp_start_byte, k = 0; j <= tmp_end_byte; j++, k++){
               tmp_letter_union_array_left[j] = letter_union_array_right[k] | set_entry_array[crossing_sets[i]].set[k];
               tmp_letter_union_array_right[j] = 0;
            }
            for(j = 0, k = tmp_start_byte; j < letters_bitmap_byte_size; j++, k++)
               tmp_letter_union_array_right[j] = tmp_letter_union_array_left[k];

            best_sum_overlap = INT_INF; // max possible value

            for(j = 0; j <= mixed_list_size; j++){ // try all position
               sum_overlap = 0;
	
               for(k = 1; k <= mixed_list_size; k++){
                  tmp_start_byte = k * letters_bitmap_byte_size;
                  tmp_end_byte = (k + 1) * letters_bitmap_byte_size - 1;
                  if(k <= j){
                     for(l = tmp_start_byte, m = 0; l <= tmp_end_byte; l++, m++){
                        tmp_letter_union_array_left[l] = letter_union_array_left[l];
                        tmp_letter_union_array_right[l] = letter_union_array_right[l] | set_entry_array[crossing_sets[i]].set[m];
                     }
                  }else{
                     for(l = tmp_start_byte, m = 0; l <= tmp_end_byte; l++, m++){
                        tmp_letter_union_array_left[l] = letter_union_array_left[l - letters_bitmap_byte_size] | set_entry_array[crossing_sets[i]].set[m];
                        tmp_letter_union_array_right[l] = letter_union_array_right[l - letters_bitmap_byte_size];
                     }
                  }
                  Node::set_intersect(tmp_letter_union_array_left + tmp_start_byte, tmp_letter_union_array_right + tmp_start_byte, letters_bitmap_byte_size, tmp_set, is_intersect);
                  sum_overlap += this->set_size(tmp_set, letters_bitmap_byte_size);
               }
               if(sum_overlap < best_sum_overlap){
                  best_pos = j;
                  tmp_start_byte = 0;
                  tmp_end_byte = (mixed_list_size + 2) * letters_bitmap_byte_size - 1;
                  for(k = tmp_start_byte; k <= tmp_end_byte; k++){
                     best_letter_union_array_left[k] = tmp_letter_union_array_left[k];
                     best_letter_union_array_right[k] = tmp_letter_union_array_right[k];
                  }
                  best_sum_overlap = sum_overlap;
               }
            }
            // put current set into mixed_list
            tmp_list_entry.index = crossing_sets[i];
            tmp_list_entry.is_set = true;
            for(j = mixed_list_size; j > best_pos; j--)
               mixed_list[j] = mixed_list[j - 1];
            mixed_list[best_pos] = tmp_list_entry;

            tmp_start_byte = 0;
            tmp_end_byte = (mixed_list_size + 2) * letters_bitmap_byte_size - 1;
            for(j = tmp_start_byte; j <= tmp_end_byte; j++){
               letter_union_array_left[j] = best_letter_union_array_left[j];
               letter_union_array_right[j] = best_letter_union_array_right[k];
            }

            mixed_list_size++;
         }
      }
    
      // process subtrees
      for(i = 0; i < mixed_list_size; i++){
         if(mixed_list[i].is_set){
            sorted_set_list[list_size] = mixed_list[i].index;
            list_size++;
         }else{
            decide_order_by_aux_tree(tree, mixed_list[i].index, letters_bitmap_byte_size, set_entry_array, set_entry_array_size, sorted_set_list, list_size); 
         }
      }
	  
	  
	delete []mixed_list;
	delete []letter_union_array_left;
	delete []letter_union_array_right;
	delete []crossing_sets;
	delete []bit_array_sets;
	delete []tmp_letter_union_array_left;
	delete []tmp_letter_union_array_right;
	delete []best_letter_union_array_left;
	delete []best_letter_union_array_right;  
	  	  
	  
	  
   }
}            
            
         



void Dir_node::logDMBR( )
{

	unsigned char		tmp_DMBR[DMBR_SIZE];

	cal_DMBR(entriesVec, count, tmp_DMBR);

	logO.log2File("dir edge ");Node::log_DMBR(tmp_DMBR);

} ////end of function

void Dir_node::logSingleChildDMBR(const unsigned char* const DMBR)
{

	//unsigned char		tmp_DMBR[DMBR_SIZE];



	//cal_DMBR(entriesVec, count, tmp_DMBR);

	logO.log2File("dir edge ");Node::log_DMBR(DMBR);

} ////end of function





void Dir_node::log2file( )
{

	unsigned char		tmp_DMBR[DMBR_SIZE];



	cal_DMBR(entriesVec, count, tmp_DMBR);


	unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
	int size_of_alphabet2 = 0;


	for(int i=0;i<DIM;i++)
	{
	
		
		bitmap_to_letters(&tmp_DMBR[i*BYTES_PER_DIM_IN_DMBR], BYTES_PER_DIM_IN_DMBR, alphabet2, size_of_alphabet2);
		
		for(int j=0;j<size_of_alphabet2;j++)
		{
			logO.log2File(alphabet2[j]); logO.log2File(" ");
			
		}

		for(int j=size_of_alphabet2;j<10;j++)

		{
			logO.log2File("  ");
			
		}


	
		logO.log2File(";\t");	
	}



		logO.log2File("\n");

} ////end of function


void Dir_node::log_OneEntry_OnOneDscDim( Dir_entry*  oneEntry,int entryIndex, int dimNum,int* alphabet_sizes )
{


	unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
	int size_of_alphabet2 = 0;


	assert(dimNum<DIM);

	Node::bitmap_to_letters(&((*oneEntry).DMBR[dimNum*BYTES_PER_DIM_IN_DMBR]), BYTES_PER_DIM_IN_DMBR, alphabet2, size_of_alphabet2);


	for(int j=0;j<size_of_alphabet2;j++)
	{
		logO.log2File((int)alphabet2[j]);logO.log2File(" ");
	}

} ////end of function


void Dir_node::print_OneEntry_OnOneDscDim( Dir_entry*  oneEntry,int entryIndex, int dimNum,int* alphabet_sizes )
{


	unsigned char alphabet2[MAX_ALPHABET_SIZE]; // letters that appear in all components.  May be a subset of the current alphabet
	int size_of_alphabet2 = 0;


	assert(dimNum<DIM);

	Node::bitmap_to_letters(&((*oneEntry).DMBR[dimNum*BYTES_PER_DIM_IN_DMBR]), BYTES_PER_DIM_IN_DMBR, alphabet2, size_of_alphabet2);



	cout<<"dir entry "<<fixed<<setw(3)<<entryIndex<<" on dim "<<dimNum<<":";

	for(int j=0;j<size_of_alphabet2;j++)
		cout<<(int)alphabet2[j]<<" ";

	cout<<"\n";

} ////end of function


//derived from getCompressedDiskSize
//do not consider DIR_NODE_OVERHEAD
int Dir_node::getCompressedEntriesSize(vector<Dir_entry> & tmp_entries, int *alphabet_sizes)
{

	int totalBytes=tmp_entries.size()*(DMBR_SIZE + sizeof(unsigned int)+usedBytesForEntryHeader);

	for (int e =0;e<tmp_entries.size();e++)
	{
		for(int i =0;i<DIM;i++)
		{
			int tmp_sum = 0;
			for(int j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
				tmp_sum += bit_1_count_lut[tmp_entries.at(e).DMBR[j]];

			if(tmp_sum==alphabet_sizes[i])
				totalBytes -= BYTES_PER_DIM_IN_DMBR;
		}
	}

	return totalBytes;

}





// check the size of bytes after compression
int Dir_node::getCompressedDiskSize(vector<Dir_entry> & tmp_entries, int *alphabet_sizes)
{

	int totalBytes=DIR_NODE_OVERHEAD+tmp_entries.size()*(DMBR_SIZE + sizeof(unsigned int)+usedBytesForEntryHeader);

	for (int e =0;e<tmp_entries.size();e++)
	{
		for(int i =0;i<DIM;i++)
		{
			int tmp_sum = 0;
			for(int j = this->DMBR_start_byte_lut[i]; j <= this->DMBR_end_byte_lut[i]; j++)
				tmp_sum += bit_1_count_lut[tmp_entries.at(e).DMBR[j]];

			if(tmp_sum==alphabet_sizes[i])
				totalBytes -= BYTES_PER_DIM_IN_DMBR;
		}
	}

	return totalBytes;

}



#endif
