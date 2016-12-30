#ifndef ND_TREE_H
#define ND_TREE_H

#include "logClass.h"
#include "Dir_node.h"
#include <stack>
#include <queue>
#include <map>

const int ND_FILE_NAME_LENGTH = 40; // Define the maximum length of the ND-tree file name
extern logClass logO;
extern int debug_height;
/*
DIM, DISK_BLOCK_SIZE and DMBR_SIZE need to be specified at the compilation time
   DIM: Number of dimensions
   DMBR_SIZE: The size (in bytes) of one DMBR

   Note that signature size is not needed.  It equals DIM, since we use 
      one byte for each character in the key.  Hence we do not allow the 
      alphabet size to be greater than 256
*/

class ND_tree{
public:
   ND_tree();
   ~ND_tree();

   void create_empty_tree(int* alphabet_s, double dir_min_u, double leaf_min_u, string ND_file_n);
   void read_existing_tree(string ND_file_n);
   //Error_code exact_query(Leaf_entry &query_data, int &number_of_io);
   Error_code exact_query_use_link(Leaf_entry &query_data, int &number_of_io);


   Error_code range_query_by_Hamming_dist(const Leaf_entry& query_data, int range, Leaf_entry* query_results, int &number_of_query_results, int &number_of_io,int &number_of_dist_computation);
   Error_code range_query_by_Hamming_dist_knn(const Leaf_entry& query_data, int range, Leaf_entry* query_results, int &number_of_query_results, int &number_of_io,int &number_of_dist_computation,priority_queue<double>& knns,vector<map<char,double> >& freq);
   Error_code insert_use_link(Leaf_entry &new_data, int &number_of_io);

   void print_information();

   unsigned int get_tree_size();

	void box_query( const Dir_entry box_query_data,	Leaf_entry* query_results, int &number_of_query_results,int & number_of_io);


   // public member
   static const int DIR_NODE_SIZE =  (DISK_BLOCK_SIZE - DIR_NODE_OVERHEAD) / (DMBR_SIZE + sizeof(unsigned int));
 


   //static const int LEAF_NODE_SIZE = (DISK_BLOCK_SIZE - LEAF_NODE_OVERHEAD) / (DIM + sizeof(ND_tree_record));
	//static const int LEAF_NODE_SIZE = (DISK_BLOCK_SIZE - LEAF_NODE_OVERHEAD) / (DIM + sizeof(float)*(TOTAL_DIM-DIM)+ sizeof(ND_tree_record));/** devidor should be the same as sizeof(Leaf_entry) **/

	//use link
	static const int LEAF_NODE_SIZE =  (DISK_BLOCK_SIZE - LEAF_NODE_OVERHEAD) / (DIM + sizeof(float)*(TOTAL_DIM-DIM)+ sizeof(ND_tree_record) );/** devidor should be the same as sizeof(Leaf_entry) **/
	//static const int LEAF_NODE_SIZE = (DISK_BLOCK_SIZE - LEAF_NODE_OVERHEAD-sizeof(int)/*pointer to link nodes*/) / (DIM +  sizeof(ND_tree_record)*2);/** devidor should be the same as sizeof(Leaf_entry) **/


   // debug methods
   void verify();
private:
   unsigned int get_block_number(); // Get the next available block number
   unsigned int choose_subtree(Leaf_entry &new_data, int &number_of_io, stack<unsigned int>& node_stack, stack<int>& index_stack);

   void read_tree_info(fstream& ND_file); // Read the tree control data into disk
   void write_tree_info(fstream& ND_file); // Write the tree control data into disk

   // Data members
   int alphabet_sizes[DIM]; // Dynamic array of alphabet sizes

   double dir_min_util;
   double leaf_min_util;

   unsigned int max_block_number; // The maximal block number now used in the tree, always increasing, start from 0
   unsigned int num_of_available_block_numbers; // Number of available block numbers stored in stack available_block_numbers
   stack<unsigned int> available_block_numbers; // The block numbers (<= max_block_number) of the nodes being deleted and become available again

   unsigned int tree_size; // Total number of entries inserted into the tree
   int height; // The height of the tree

   char ND_file_name[ND_FILE_NAME_LENGTH]; // The name of the ND-tree file
   fstream ND_file;

   unsigned int root; // Root block number
};



ND_tree::ND_tree(){
}



ND_tree::~ND_tree(){
   if(ND_file.is_open()){
      write_tree_info(ND_file);
      ND_file.close();
   }
}


void ND_tree::create_empty_tree(
     int* alphabet_s, double dir_min_u, double leaf_min_u, string ND_file_n)
{
   if(ND_file.is_open()){ // close previous tree
      write_tree_info(ND_file);
      ND_file.close();
   }

   int i;
   for(i = 0; i < DIM; i++) alphabet_sizes[i] = alphabet_s[i];

   dir_min_util = dir_min_u;
   leaf_min_util = leaf_min_u;

   max_block_number = 0;
   num_of_available_block_numbers = 0;
   tree_size = 0;
   height = 1; // The empty root

   int file_name_length = (int)ND_file_n.copy(ND_file_name,(int) ND_file_n.length());
   ND_file_name[file_name_length] = '\0'; // Append a null character
   ND_file.open(ND_file_name, ios_base::binary | ios_base::out);

   root = get_block_number();
   Leaf_node* root_node = new Leaf_node(alphabet_sizes);
   root_node->write_node(ND_file, root);
   delete root_node;

   write_tree_info(ND_file);

   ND_file.close();
   ND_file.open(ND_file_name, ios_base::binary | ios_base::in | ios_base::out); 
}



void ND_tree::read_existing_tree(string ND_file_n)
{
   if(ND_file.is_open()){ // close previous tree
      write_tree_info(ND_file);
      ND_file.close();
   }

   int file_name_length = (int)ND_file_n.copy(ND_file_name, (int)ND_file_n.length());
   ND_file_name[file_name_length] = '\0'; // Append a null character
   ND_file.open(ND_file_name, ios_base::binary | ios_base::in | ios_base::out); 
   read_tree_info(ND_file);
}


Error_code ND_tree::exact_query_use_link(Leaf_entry &query_data, int &number_of_io){
   int i;
   Error_code result;

   number_of_io = 0;
   if(height == 1){ // root is a leaf
      Leaf_node* root_node = new Leaf_node(alphabet_sizes);
      root_node->read_node(ND_file, root);
      number_of_io++;
      Error_code result = root_node->retrieve(query_data);
	  if(result == success)
		root_node->write_node(ND_file, root);
      delete root_node;
      return result;
   }else{
      stack<unsigned int> dir_check_list; // dir_check_list contains a list of block number of dir nodes that may contain the query point
      stack<int> height_list; // The height of the dir node corresponds to dir_check_list
      unsigned int check_block_no;
      int tmp_height;
      dir_check_list.push(root);
      height_list.push(height);

      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);


	  vector<int> covered_block_numbers;
      bool not_found = true;
      while(!dir_check_list.empty() && not_found){
         check_block_no = dir_check_list.top(); // fetch the fisrt item from the check list
         dir_check_list.pop(); // remove the item from the list
         tmp_height = height_list.top();
         height_list.pop();
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, check_block_no);
		 else
			d_node->read_node_static(ND_file, check_block_no);
         number_of_io++;

		covered_block_numbers.clear();
         d_node->find_covering_entries(query_data, covered_block_numbers);

         // check if the dir node is one level higher than the leaves
         if(tmp_height == 2) // One level higher than leaf
            for(i = 0; i < covered_block_numbers.size(); i++){ // go through all the covering entries
               l_node->read_node(ND_file, covered_block_numbers.at(i));
               number_of_io++;
               result = l_node->retrieve(query_data);
               if(result == success){
				  //l_node->write_node(ND_file, covering_block_numbers[i]);
                  not_found = false;
                  break;
               }
            }
         else{
            assert(tmp_height > 2);
            for(i = 0; i < covered_block_numbers.size(); i++){ // go through all the covering entries
               dir_check_list.push(covered_block_numbers.at(i));
               height_list.push(tmp_height - 1);
            }
         }
      }
      delete l_node;
      delete d_node;

      if(not_found)return not_present;
      else return success;
   }
}

void ND_tree::box_query(
	const Dir_entry box_query_data,
	Leaf_entry* query_results, 
	int &number_of_query_results,
	int &	number_of_io)
{
	int i;   
	number_of_query_results = 0;
	number_of_io = 0;


	if(height == 1)
	{ // root is a leaf

      Leaf_node* root_node = new Leaf_node(alphabet_sizes);
      root_node->read_node(ND_file, root);


		number_of_io++; 

		//if(USE_LINK_FOR_BOX_QUERY==1)
		//	//root_node->retrieve_by_box_query_use_link(box_query_data.DMBR,  query_results, number_of_query_results,number_of_io);

			root_node->retrieve_by_box_query_use_link_v2(box_query_data.DMBR,  query_results, number_of_query_results,number_of_io);
		//else
		//	root_node->retrieve_by_box_query(box_query_data.DMBR,  query_results, number_of_query_results);

		delete root_node;
		return;
	}
	else
	{

		stack<unsigned int> dir_check_list; // dir_check_list contains a list of block number of dir nodes that may contain the query point
		stack<int> height_list; // The height of the dir node corresponds to dir_check_list
		unsigned int check_block_no;
		int tmp_height;
		dir_check_list.push(root);
		height_list.push(height);

		Dir_node* d_node = new Dir_node(alphabet_sizes);
		Leaf_node* l_node = new Leaf_node(alphabet_sizes);

		vector<int> block_numbers_in_box;
		while(!dir_check_list.empty())
		{
            check_block_no = dir_check_list.top(); // fetch the fisrt item from the check list
			dir_check_list.pop(); // remove the item from the list
            tmp_height = height_list.top();
			height_list.pop();

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, check_block_no);
		 else
			d_node->read_node_static(ND_file, check_block_no);

			number_of_io++;


			block_numbers_in_box.clear();

			d_node->find_entries_by_box_query( box_query_data, block_numbers_in_box);

			// check if the dir node is one level higher than the leaves
			if(tmp_height == 2) // One level higher than leaf
				for(i = 0; i < block_numbers_in_box.size(); i++)
				{ // go through all the relevant entries

					l_node->read_node(ND_file, block_numbers_in_box.at(i));
					number_of_io++;

					l_node->retrieve_by_box_query_use_link_v2(box_query_data.DMBR,  query_results, number_of_query_results,number_of_io);
	

				}
			else{
				assert(tmp_height > 2);
				for(i = 0; i < block_numbers_in_box.size(); i++)
				{ // go through all the relevant entries
					dir_check_list.push(block_numbers_in_box.at(i));
					height_list.push(tmp_height - 1);
				}
			}
		}
		//delete []covering_block_numbers;
		delete l_node;
		delete d_node;

	}



}


Error_code ND_tree::range_query_by_Hamming_dist_knn(const Leaf_entry& query_data, int range, Leaf_entry* query_results, int &number_of_query_results, int &number_of_io,int &number_of_dist_computation,priority_queue<double>& knns,vector<map<char,double> > & freq){
   int i;
   Error_code result;

   double r = 999999.0;
   number_of_query_results = 0;
   number_of_io = 0;
   number_of_dist_computation = 0;

   if(height == 1){ // root is a leaf
      Leaf_node* root_node = new Leaf_node(alphabet_sizes);
      root_node->read_node(ND_file, root);
      number_of_io++;
      Error_code result = root_node->retrieve_by_hamming_dist_knn(query_data, range, query_results, number_of_query_results, number_of_io, number_of_dist_computation,r,knns,freq);
      delete root_node;
      return result;
   }else{
      stack<unsigned int> dir_check_list; // dir_check_list contains a list of block number of dir nodes that may contain the query point
      stack<int> height_list; // The height of the dir node corresponds to dir_check_list
      unsigned int check_block_no;
      int tmp_height;
      dir_check_list.push(root);
      height_list.push(height);

      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);

	  vector<int> block_numbers_in_range;

      while(!dir_check_list.empty()){
         check_block_no = dir_check_list.top(); // fetch the fisrt item from the check list
         dir_check_list.pop(); // remove the item from the list
         tmp_height = height_list.top();
         height_list.pop();
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, check_block_no);
		 else
			d_node->read_node_static(ND_file, check_block_no);

         number_of_io++;

		 block_numbers_in_range.clear();
         d_node->find_entries_by_hamming_dist(query_data, range, block_numbers_in_range,number_of_dist_computation);

         // check if the dir node is one level higher than the leaves
         if(tmp_height == 2) // One level higher than leaf
            for(i = 0; i < block_numbers_in_range.size(); i++){ // go through all the relevant entries
               l_node->read_node(ND_file, block_numbers_in_range.at(i));
               number_of_io++;
               result = l_node->retrieve_by_hamming_dist_knn(query_data, range, query_results, number_of_query_results, number_of_io,number_of_dist_computation,r,knns,freq);
            }
         else{
            assert(tmp_height > 2);
            for(i = 0; i < block_numbers_in_range.size(); i++){ // go through all the relevant entries
               dir_check_list.push(block_numbers_in_range.at(i));
               height_list.push(tmp_height - 1);
            }
         }
      }
      //delete []covering_block_numbers;
      delete l_node;
      delete d_node;

      if(number_of_query_results == 0)return not_present;
      else return success;
   }
}

Error_code ND_tree::range_query_by_Hamming_dist(const Leaf_entry& query_data, int range, Leaf_entry* query_results, int &number_of_query_results, int &number_of_io,int &number_of_dist_computation){
   int i;
   Error_code result;

   number_of_query_results = 0;
   number_of_io = 0;
   number_of_dist_computation = 0;

   if(height == 1){ // root is a leaf
      Leaf_node* root_node = new Leaf_node(alphabet_sizes);
      root_node->read_node(ND_file, root);
      number_of_io++;
      Error_code result = root_node->retrieve_by_hamming_dist(query_data, range, query_results, number_of_query_results, number_of_io, number_of_dist_computation);
      delete root_node;
      return result;
   }else{
      stack<unsigned int> dir_check_list; // dir_check_list contains a list of block number of dir nodes that may contain the query point
      stack<int> height_list; // The height of the dir node corresponds to dir_check_list
      unsigned int check_block_no;
      int tmp_height;
      dir_check_list.push(root);
      height_list.push(height);

      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);

	  vector<int> block_numbers_in_range;

      while(!dir_check_list.empty()){
         check_block_no = dir_check_list.top(); // fetch the fisrt item from the check list
         dir_check_list.pop(); // remove the item from the list
         tmp_height = height_list.top();
         height_list.pop();
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, check_block_no);
		 else
			d_node->read_node_static(ND_file, check_block_no);

         number_of_io++;

		 block_numbers_in_range.clear();
         d_node->find_entries_by_hamming_dist(query_data, range, block_numbers_in_range,number_of_dist_computation);

         // check if the dir node is one level higher than the leaves
         if(tmp_height == 2) // One level higher than leaf
            for(i = 0; i < block_numbers_in_range.size(); i++){ // go through all the relevant entries
               l_node->read_node(ND_file, block_numbers_in_range.at(i));
               number_of_io++;
               result = l_node->retrieve_by_hamming_dist(query_data, range, query_results, number_of_query_results, number_of_io,number_of_dist_computation);
            }
         else{
            assert(tmp_height > 2);
            for(i = 0; i < block_numbers_in_range.size(); i++){ // go through all the relevant entries
               dir_check_list.push(block_numbers_in_range.at(i));
               height_list.push(tmp_height - 1);
            }
         }
      }
      //delete []covering_block_numbers;
      delete l_node;
      delete d_node;

      if(number_of_query_results == 0)return not_present;
      else return success;
   }
}

//if(USE_LINK_FOR_BOX_QUERY==1) then call this function

Error_code ND_tree::insert_use_link(Leaf_entry &new_data, int &number_of_io){
/////////////////////////
//	/*debug*/
//	if(height>1)
//	{
//		Dir_node* debug_node= new Dir_node(alphabet_sizes);
//		debug_node->read_node(ND_file, 3);
//		cout<<debug_node->get_node_count()<<endl;
//		delete debug_node;
//	}
/////////////////////////
   if(exact_query_use_link(new_data, number_of_io) == success)
   {
      return duplicate_error;
   
   }

   Error_code result;
   stack<unsigned int> ins_path_nodes; // sequence of dir node block nos from root
   stack<int> ins_path_child_indices; // sequence of child indices from root
   unsigned int insert_block_no = choose_subtree(new_data, number_of_io, ins_path_nodes, ins_path_child_indices);


   Leaf_node* l_node = new Leaf_node(alphabet_sizes);
   l_node->read_node(ND_file, insert_block_no);
   number_of_io++;
   unsigned char cur_DMBR[DMBR_SIZE], new_DMBR[DMBR_SIZE]; // Digital Minimum Bounding Rectangle.  
   Leaf_node* new_l_node = NULL; // used if split
#ifdef LOG_SPLITTING_HEIGHT
		 splitting_height = 1;
#endif
   result = l_node->insert_new_data(new_data, leaf_min_util, alphabet_sizes, cur_DMBR, new_l_node, new_DMBR); // new_l_node and new_DMBR used only if split
   l_node->write_node(ND_file, insert_block_no);
   number_of_io++;



   // overflow and adjust DMBR
   Dir_node* d_node = new Dir_node(alphabet_sizes);
   unsigned int parent_block_no;
   int child_index;
   int cur_level;
   bool need_to_change;
   if(result == success){ // no overflow.  Adjust DMBR only
      cur_level = 1; // leaf
      need_to_change = true;
      while(cur_level < height && need_to_change){ // leaf is not root
         // get the parent
         parent_block_no = ins_path_nodes.top();
         ins_path_nodes.pop();
         child_index = ins_path_child_indices.top();
         ins_path_child_indices.pop();

         // Modify the DMBR in parent
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, parent_block_no);
		 else
			d_node->read_node_static(ND_file, parent_block_no);

         number_of_io++;
         d_node->set_DMBR(child_index, cur_DMBR, need_to_change);
		 
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->write_node_dynamic(ND_file, parent_block_no);
		 else
			d_node->write_node_static(ND_file, parent_block_no);
         number_of_io++;

         cur_level++;
      }
   }else{ // overflows
      Dir_node* new_d_node = NULL;

      unsigned int new_block_no;
      // get a block no for the new leaf 
      new_block_no = get_block_number();
      new_l_node->write_node(ND_file, new_block_no);

#ifdef LOG_OLD_AND_NEW_BLOCK_LEAF
	  logO.log2File("new leaf node block:");logO.log2File((int)new_block_no);logO.log2File("\n");

	  logO.log2File("old leaf node block:");logO.log2File((int)insert_block_no);logO.log2File("\n");
	  //l_node->log2file();
	  //logO.log2File("new leaf node block:");logO.log2File((int)new_block_no);logO.log2File("\n");
	  //new_l_node->log2file();
#endif


	//l_node->log2file();
	//new_l_node->log2file();

#ifdef LOG_ENTRY_ON_SPLITTING_DIM

	logO.log2File("level: ");logO.log2File((int)1);logO.log2File("\n");

#endif   


      number_of_io++;
      delete new_l_node;

      bool overflow_now = true;
      cur_level = 1; // leaf level
      while(overflow_now && cur_level < height){ // not root yet
		  // get the parent and child index
		  parent_block_no = ins_path_nodes.top();
		  ins_path_nodes.pop();
		  child_index = ins_path_child_indices.top();
		  ins_path_child_indices.pop();


		  if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			  d_node->read_node_dynamic(ND_file, parent_block_no);
		  else
			  d_node->read_node_static(ND_file, parent_block_no);



#ifdef ShorterThanSplitDim
logO.log2File("height: ");logO.log2File(cur_level);logO.log2File("\n");
		  unsigned char twoNewNode_DMBR[DMBR_SIZE],otherNodes_DMBR[DMBR_SIZE],d_tmp_DMBR[DMBR_SIZE];

		  Node aNode(A);

		  for(int d = 0; d < DMBR_SIZE; d++)
			  twoNewNode_DMBR[d]=otherNodes_DMBR[d]=0;

		  for(int d = 0; d < DMBR_SIZE; d++)
		  {
			  twoNewNode_DMBR[d]|=cur_DMBR[d];
			  twoNewNode_DMBR[d]|=new_DMBR[d];
		  }

		  logO.log2File("length: \t");
		  aNode.log_DMBR(twoNewNode_DMBR);


		  for(int t=0;t<d_node->get_node_count();t++)
		  {
			  if(t!=child_index)
			  {
				  d_node->get_DMBR(d_tmp_DMBR,t);
				  for(int d = 0; d < DMBR_SIZE; d++)
					  otherNodes_DMBR[d]|=d_tmp_DMBR[d];
			  }
		  }

		  for(int d = 0; d < DMBR_SIZE; d++)
			  d_tmp_DMBR[d]=otherNodes_DMBR[d]&twoNewNode_DMBR[d];
		  
		 //aNode.log_OneDirEntry_OnOneDscDim(twoNewNode_DMBR,0,A);logO.log2File("\n");
		 //aNode.log_OneDirEntry_OnOneDscDim(otherNodes_DMBR,0,A);logO.log2File("\n");
logO.log2File("overlap: \t");

		  aNode.log_DMBR(d_tmp_DMBR);

#endif




         number_of_io++;
         d_node->set_DMBR(child_index, cur_DMBR, need_to_change);

//////////////////debug
//		 cout<<d_node->get_node_count();
//////////////////////
#ifdef LOG_SPLITTING_HEIGHT
		 splitting_height = cur_level+1;
#endif
         result = d_node->insert_new_entry(new_block_no, new_DMBR, dir_min_util, alphabet_sizes, cur_DMBR, new_d_node, new_DMBR);

//////////////////debug
//		 cout<<d_node->get_node_count();
//////////////////////


		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->write_node_dynamic(ND_file, parent_block_no);
		 else
			d_node->write_node_static(ND_file, parent_block_no);
         number_of_io++;
         if(result == success){ // no overflow
            overflow_now = false;
         }else{ // overflow
            // get a block no for the new dir node
            new_block_no = get_block_number();
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			new_d_node->write_node_dynamic(ND_file, new_block_no);
		 else
			new_d_node->write_node_static(ND_file, new_block_no);
            number_of_io++;
            delete new_d_node;

#ifdef LOG_ENTRY_ON_SPLITTING_DIM

	logO.log2File("level: ");logO.log2File((int)(cur_level+1));logO.log2File("\n");

#endif   

         }
         cur_level++;
      }


      if(overflow_now){ // cur_level == height
         // need to create a new root
         Dir_node* new_root_node = new Dir_node(alphabet_sizes);
         new_root_node->set_node_count(2);
         int new_root_block_no = get_block_number();
         if(height > 1) 
            new_root_node->set_node_entry(0, parent_block_no, cur_DMBR);
            // new_root_node->insert_new_entry(parent_block_no, cur_DMBR, dir_min_util, alphabet_sizes, NULL, NULL, NULL); // insert the old root as the 1 st entry.  Knowing it will not overflow
         else // leaf is split
            new_root_node->set_node_entry(0, insert_block_no, cur_DMBR);
            //new_root_node->insert_new_entry(insert_block_no, cur_DMBR, dir_min_util, alphabet_sizes, NULL, NULL, NULL); // insert the old root as the 1 st entry.  Knowing it will not overflow
         new_root_node->set_node_entry(1, new_block_no, new_DMBR);
         //new_root_node->insert_new_entry(new_block_no, new_DMBR, dir_min_util, alphabet_sizes, NULL, NULL, NULL); // insert the new node from split.  Knowing it will not overflow
                  
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			new_root_node->write_node_dynamic(ND_file, new_root_block_no);
		 else
			new_root_node->write_node_static(ND_file, new_root_block_no);
         number_of_io++;
         delete new_root_node;

//////////////////////////
///*debug*/
//Dir_node* debug_node= new Dir_node(alphabet_sizes);
//debug_node->read_node(ND_file, 3);
//cout<<debug_node->get_node_count()<<endl;
//delete debug_node;
/////////////////////////

         root = new_root_block_no;
         height++;
      }else{ // no overflow, continue adjusting DMBR
         need_to_change = true;
         while(cur_level < height && need_to_change){ // leaf is not root
            // get the parent
            parent_block_no = ins_path_nodes.top();
            ins_path_nodes.pop();
            child_index = ins_path_child_indices.top();
            ins_path_child_indices.pop();

            // Modify the DMBR in parent

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, parent_block_no);
		 else
			d_node->read_node_static(ND_file, parent_block_no);

            number_of_io++;
            d_node->set_DMBR(child_index, cur_DMBR, need_to_change);
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->write_node_dynamic(ND_file, parent_block_no);
		 else
			d_node->write_node_static(ND_file, parent_block_no);
            number_of_io++;

            cur_level++;
         }
      }
   }
   tree_size++;

   delete l_node;
   delete d_node;
   return success;
}





void ND_tree::print_information(){
   cout << endl;
   cout << "Dsc dim: " << DIM << endl;
   cout << "Cnt dim: " << CNTDIM << endl;
   if(TREE_TYPE != DYNAMIC_TREE)
	cout << "Dir node size: " << DIR_NODE_SIZE << endl;
   cout << "Leaf node size: " << LEAF_NODE_SIZE << endl;
   cout << endl;
   cout << "Dir node minimum utilization: " << dir_min_util << endl;
   cout << "Leaf node minimum utilization: " << leaf_min_util << endl;
   cout << endl;
   cout << "Number of disk blocks: " << max_block_number - num_of_available_block_numbers << endl;
   cout << "Number of vectors indexed: " << tree_size << endl;
   cout << "Height: " << height << endl;
   cout << endl;

   cout<<"nodeSplitType "<<nodeSplitType<<endl;
   cout<<"DISK_BLOCK_SIZE "<<DISK_BLOCK_SIZE<<endl;

   if(nodeSplitType  == TO_DEATH_RANDOM_UTIL)
   {
	   cout<<  "heuristics_overlap_used_dir after knapsack: " << heuristics_overlap_used_dir << endl;
	   cout<<  "heuristics_area_used_dir after knapsack: " << heuristics_area_used_dir << endl;
	   cout<<  "heuristics_overlap_used_leaf after knapsack: " << heuristics_overlap_used_leaf << endl;
	   cout<<  "heuristics_area_used_leaf after knapsack: " << heuristics_area_used_leaf << endl;
   }
	cout<<  "BestChild heuristics used for covered_area: " << BestChild_covered_area << endl;
	cout<<  "BestChild heuristics used for notcovered_overlap_enlarge: " << BestChild_notcovered_overlap_enlarge << endl;
	cout<<  "BestChild heuristics used for notcovered_area_enlarge:    " << BestChild_notcovered_area_enlarge << endl;
	cout<<  "BestChild heuristics used for notcovered_area:    " << BestChild_notcovered_area << endl;

   cout << endl;

   int i;

   int dir_node_count, leaf_node_count;
   int dir_node_entry_count, leaf_node_entry_count;
   double dir_util, leaf_util, util;

   if(height == 1){ // root is a leaf
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);
      l_node->read_node(ND_file, root);

      leaf_node_entry_count = l_node->get_node_count(); 
      dir_node_count = 0; 
      leaf_node_count = 1; // the root

      dir_util = -1; // n/a
      leaf_util = static_cast<double>(leaf_node_entry_count) / (leaf_node_count * LEAF_NODE_SIZE);
      util = leaf_util;

      cout << "Number of dir nodes: " << dir_node_count << endl;
      cout << "Number of leaf nodes: " << leaf_node_count << endl;
      cout << endl;
      cout << "Dir node utilization: " << dir_util << endl;
      cout << "Leaf node utilization: " << leaf_util << endl;
      cout << "Utilization: " << util << endl;
      cout << endl;

      delete l_node;
      return;
   }
   else if(height == 2){
      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);

      int tmp_count;
      Dir_entry tmp_d_entry;
      int child_count;

      dir_node_count = 0; // root is not counted
      leaf_node_entry_count = 0;

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, root);
		 else
			d_node->read_node_static(ND_file, root);
#ifdef LOG_EDGE_LENGTH_DIR
		 logO.log2File("height: ");logO.log2File(2); logO.log2File(" ");
		 d_node->logDMBR();
#endif
      child_count = d_node->get_node_count();
      leaf_node_count = child_count;

      for(i = 0; i < child_count; i++){ // go through all the children
         d_node->get_node_entry(i, tmp_d_entry);
         l_node->read_node(ND_file, tmp_d_entry.child);
#ifdef LOG_EDGE_LENGTH_LEAF
		 logO.log2File("height: ");logO.log2File(1); logO.log2File(" ");

	  l_node->logDMBR();
#endif
         tmp_count = l_node->get_node_count();
         leaf_node_entry_count += tmp_count;
      }

      dir_util = -1; 
      leaf_util = static_cast<double>(leaf_node_entry_count) / (leaf_node_count * LEAF_NODE_SIZE);
      util = leaf_util;

      cout << "Number of dir nodes: " << dir_node_count << " (root is excluded)" << endl;
      cout << "Number of leaf nodes: " << leaf_node_count << endl;
      cout << endl;
      cout << "Dir node utilization: " << dir_util << endl;
      cout << "Leaf node utilization: " << leaf_util << endl;
      cout << "Utilization: " << util << endl;
      cout << endl;

      delete l_node;
      delete d_node;
   }
   else{
      stack< Dir_entry > dir_check_list; // dir_check_list contains a list of block number of dir nodes that may contain the query point
      stack<int> height_list; // The height of the dir node corresponds to dir_check_list

      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);

      Dir_entry check_d_entry;
      int tmp_height;
      int tmp_count;
      Dir_entry tmp_d_entry;
      int child_count;

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, root);
		 else
			d_node->read_node_static(ND_file, root);
      tmp_count = d_node->get_node_count();
//#ifdef LOG_UTILIZATION
//	  logO.log2File("dir Util ");logO.log2File((float)tmp_count/(float)DIR_NODE_SIZE);logO.log2File("\n");
//#endif
#ifdef LOG_EDGE_LENGTH_DIR
		 logO.log2File("height: ");logO.log2File(height); logO.log2File(" ");
	  d_node->logDMBR();
#endif
      for(i = 0; i < tmp_count; i++){
         d_node->get_node_entry(i, tmp_d_entry);
         dir_check_list.push(tmp_d_entry);
         height_list.push(height - 1);
      }

      dir_node_count = 0; // root is not counted
      dir_node_entry_count = 0;
      leaf_node_count = 0;
      leaf_node_entry_count = 0;

      while(!dir_check_list.empty()){
         check_d_entry = dir_check_list.top(); // fetch the fisrt item from the check list
         dir_check_list.pop(); // remove the item from the list
         tmp_height = height_list.top();
         height_list.pop();

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, check_d_entry.child);
		 else
			d_node->read_node_static(ND_file, check_d_entry.child);
         tmp_count = d_node->get_node_count();
//#ifdef LOG_UTILIZATION
//	  logO.log2File("dir Util ");logO.log2File((float)tmp_count/(float)DIR_NODE_SIZE);logO.log2File("\n");
//#endif
#ifdef LOG_EDGE_LENGTH_DIR
		 logO.log2File("height: ");logO.log2File(tmp_height); logO.log2File(" ");
	  d_node->logDMBR();
#endif
         dir_node_count++;
         dir_node_entry_count += tmp_count;

         // check if the dir node is one level higher than the leaves
         if(tmp_height == 2){ // One level higher than leaf
            child_count = d_node->get_node_count();
            leaf_node_count += child_count;

            for(i = 0; i < child_count; i++){ // go through all the children
               d_node->get_node_entry(i, tmp_d_entry);
               l_node->read_node(ND_file, tmp_d_entry.child);

               tmp_count = l_node->get_node_count();
#ifdef LOG_UTILIZATION
	  logO.log2File("leaf Util ");logO.log2File((float)tmp_count/(float)LEAF_NODE_SIZE);logO.log2File("\n");
#endif

#ifdef LOG_EDGE_LENGTH_LEAF
		 logO.log2File("height: ");logO.log2File(1); logO.log2File(" ");
	  l_node->logDMBR();
#endif
               leaf_node_entry_count += tmp_count;
            }
         }
         else{
            assert(tmp_height > 2);
            child_count = d_node->get_node_count();
            for(i = 0; i < child_count; i++){ // go through all the children
               d_node->get_node_entry(i, tmp_d_entry);
               dir_check_list.push(tmp_d_entry);
               height_list.push(tmp_height - 1);
            }
         }
      }

      dir_util = static_cast<double>(dir_node_entry_count) / (dir_node_count * DIR_NODE_SIZE); 
      leaf_util = static_cast<double>(leaf_node_entry_count) / (leaf_node_count * LEAF_NODE_SIZE);
      util = static_cast<double>(leaf_node_entry_count + dir_node_entry_count) / (leaf_node_count * LEAF_NODE_SIZE + dir_node_count * DIR_NODE_SIZE);

      cout << "Number of dir nodes: " << dir_node_count << " (root is excluded)" << endl;
      cout << "Number of leaf nodes: " << leaf_node_count << endl;
      cout << "AVG number of entries in dir nodes: " << static_cast<double>(dir_node_entry_count)/dir_node_count << endl;
      cout << "AVG number of entries in leaf nodes: " << static_cast<double>(leaf_node_entry_count)/leaf_node_count << endl;

	  cout << endl;
	  if(TREE_TYPE != DYNAMIC_TREE)
		  cout << "Dir node utilization: " << dir_util << endl;
	  cout << "Leaf node utilization: " << leaf_util << endl;
	  if(TREE_TYPE != DYNAMIC_TREE)
		  cout << "Utilization: " << util << endl;
	  cout << endl;

	  delete l_node;
      delete d_node;
   }
}



unsigned int ND_tree::get_tree_size(){
   return tree_size;
}



void ND_tree::verify(){
   int i, j;

   if(height == 1){ // root is a leaf
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);
      l_node->read_node(ND_file, root);
      assert(tree_size == l_node->get_node_count()); 
      delete l_node;
      return;
   }
   else if(height == 2){
      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);

      int tmp_count;
      Dir_entry tmp_d_entry;
      unsigned char tmp_DMBR[DMBR_SIZE];
      int child_count;
      int object_count;

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, root);
		 else
			d_node->read_node_static(ND_file, root);

      object_count = 0;
      child_count = d_node->get_node_count();
      for(i = 0; i < child_count; i++)
	  { // go through all the children
         d_node->get_node_entry(i, tmp_d_entry);
         l_node->read_node(ND_file, tmp_d_entry.child);

         tmp_count = l_node->get_node_count();
         object_count += tmp_count;

         // check DMBR
         l_node->get_DMBR(tmp_DMBR);
         for(j = 0; j < DMBR_SIZE; j++)
            assert(tmp_DMBR[j] == tmp_d_entry.DMBR[j]);
      }
      assert(object_count == tree_size);

      delete l_node;
      delete d_node;
   }
   else{
      stack< Dir_entry > dir_check_list; // dir_check_list contains a list of block number of dir nodes that may contain the query point
      stack<int> height_list; // The height of the dir node corresponds to dir_check_list
      stack<unsigned int> parent_list; // The block number of the parent of the dir node corresponds to dir_check_list
      stack<int> index_list; // The child index in parent of the dir node corresponds to dir_check_list

      Dir_node* d_node = new Dir_node(alphabet_sizes);
      Leaf_node* l_node = new Leaf_node(alphabet_sizes);

      Dir_entry check_d_entry;
      int tmp_height;
      unsigned int tmp_parent;
      int tmp_index;
      int tmp_count;
      Dir_entry tmp_d_entry;
      unsigned char tmp_DMBR[DMBR_SIZE];
      int child_count;
      int object_count;

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, root);
		 else
			d_node->read_node_static(ND_file, root);
      tmp_count = d_node->get_node_count();
      for(i = 0; i < tmp_count; i++){
         d_node->get_node_entry(i, tmp_d_entry);
         dir_check_list.push(tmp_d_entry);
         height_list.push(height - 1);
         parent_list.push(root);
         index_list.push(i);
      }

      object_count = 0;
      while(!dir_check_list.empty()){
         check_d_entry = dir_check_list.top(); // fetch the fisrt item from the check list
         dir_check_list.pop(); // remove the item from the list
         tmp_height = height_list.top();
         height_list.pop();
         tmp_parent = parent_list.top();
         parent_list.pop();
         tmp_index = index_list.top();
         index_list.pop();

		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, check_d_entry.child);
		 else
			d_node->read_node_static(ND_file, check_d_entry.child);
         // check DMBR
         tmp_count = d_node->get_node_count();
         d_node->get_DMBR(tmp_DMBR); 
         for(i = 0; i < DMBR_SIZE; i++)
            assert(tmp_DMBR[i] == check_d_entry.DMBR[i]);

         // check if the dir node is one level higher than the leaves
         if(tmp_height == 2){ // One level higher than leaf
            child_count = d_node->get_node_count();
            for(i = 0; i < child_count; i++){ // go through all the children
               d_node->get_node_entry(i, tmp_d_entry);
               l_node->read_node(ND_file, tmp_d_entry.child);

               tmp_count = l_node->get_node_count();
               object_count += tmp_count;

               // check DMBR
               l_node->get_DMBR(tmp_DMBR);
               for(j = 0; j < DMBR_SIZE; j++)
                  assert(tmp_DMBR[j] == tmp_d_entry.DMBR[j]);
            }
         }
         else{
            assert(tmp_height > 2);
            child_count = d_node->get_node_count();
            for(i = 0; i < child_count; i++){ // go through all the children
               d_node->get_node_entry(i, tmp_d_entry);
               dir_check_list.push(tmp_d_entry);
               height_list.push(tmp_height - 1);
               parent_list.push(check_d_entry.child); 
               index_list.push(i); // n/a for root
            }
         }
      }
      assert(object_count == tree_size);

      delete l_node;
      delete d_node;
   }
}



unsigned int ND_tree::get_block_number()
/*
Get a block number for the new node
if there is available block number from stack, use it first
otherwise use a new max block number
return the block number, and modify available_block_numbers and max_block_number accordingly
note that we have no concurrency control here.
*/
{
   if(num_of_available_block_numbers == 0){
      max_block_number++;
      return max_block_number;
   }else{
      unsigned int block_number;
      block_number = available_block_numbers.top();
      available_block_numbers.pop();
      num_of_available_block_numbers--;
      return block_number;
   }
}



unsigned int ND_tree::choose_subtree(Leaf_entry &new_data, int &number_of_io, stack<unsigned int>& node_stack, stack<int>& index_stack){
   unsigned int cur_block_no = root; // point to the current node
   int cur_child_index;

   Dir_node* d_node = new Dir_node(alphabet_sizes);

#ifdef LOG_CHOOSE_SUB_TREE
	   logO.log2File(height);logO.log2File(":");logO.log2File((int)cur_block_no);logO.log2File(" ");
#endif

   int cur_height = height;
   while(cur_height > 1){

      node_stack.push(cur_block_no);
		 if((nodeSplitType!=ORIGINAL)&&(TREE_TYPE == DYNAMIC_TREE))
			d_node->read_node_dynamic(ND_file, cur_block_no);
		 else
			d_node->read_node_static(ND_file, cur_block_no);
      number_of_io++;
debug_height = cur_height;
      d_node->find_best_child(new_data, alphabet_sizes, cur_block_no, cur_child_index);

      index_stack.push(cur_child_index);

      cur_height--;
#ifdef LOG_CHOOSE_SUB_TREE
	   logO.log2File(cur_height);logO.log2File(":");logO.log2File((int)cur_block_no);logO.log2File(" ");
#endif


   }
#ifdef LOG_CHOOSE_SUB_TREE
	   logO.log2File("\n");
#endif
   delete d_node;
   return cur_block_no;
}


/* Read tree info from block 0 of the file.  stack available_block_number, if not empty,
   is read from the last block */

void ND_tree::read_tree_info(fstream& ND_file)
{
   ND_file.seekg(0);
   ND_file.read((char*)alphabet_sizes, sizeof(alphabet_sizes));
   ND_file.read((char*)(&dir_min_util), sizeof(double));
   ND_file.read((char*)(&leaf_min_util), sizeof(double));
   ND_file.read((char*)(&max_block_number), sizeof(unsigned int));
   ND_file.read((char*)(&num_of_available_block_numbers), sizeof(unsigned int));
   ND_file.read((char*)(&tree_size), sizeof(unsigned int));
   ND_file.read((char*)(&height), sizeof(int));
   ND_file.read((char*)(&root), sizeof(unsigned int));

   unsigned int dummy;
  ND_file.read((char*)(&dummy), sizeof(unsigned int));//read in DIR_NODE_SIZE and Leaf_NODE_SIZE
  ND_file.read((char*)(&dummy), sizeof(unsigned int));

   if(num_of_available_block_numbers > 0){
      unsigned int block_number;
      unsigned int count = 0;
      ND_file.seekg((max_block_number + 1) * DISK_BLOCK_SIZE);
      while(count < num_of_available_block_numbers){
         ND_file.read((char*)(&block_number), sizeof(unsigned int));
         available_block_numbers.push(block_number);
         count++;
      }
   }
}


/* write tree info to block 0 of the file.  stack available_block_number, if not empty,
   is written to the last block */

void ND_tree::write_tree_info(fstream& ND_file)
{
   ND_file.seekp(0);
   ND_file.write((const char*)alphabet_sizes, sizeof(alphabet_sizes));
   ND_file.write((const char*)(&dir_min_util), sizeof(double));
   ND_file.write((const char*)(&leaf_min_util), sizeof(double));
   ND_file.write((const char*)(&max_block_number), sizeof(unsigned int));
   ND_file.write((const char*)(&num_of_available_block_numbers), sizeof(unsigned int));
   ND_file.write((const char*)(&tree_size), sizeof(unsigned int));
   ND_file.write((const char*)(&height), sizeof(int));
   ND_file.write((const char*)(&root), sizeof(unsigned int));
	 int temp1 = 0;
	 int temp2 = 0;
  ND_file.write((const char*)(&temp1 ), sizeof(unsigned int));
  ND_file.write((const char*)(&temp2 ), sizeof(unsigned int));



   if(num_of_available_block_numbers > 0){
      unsigned int count = num_of_available_block_numbers;
      stack<unsigned int> tmp_stack = available_block_numbers;
      unsigned int block_number;
      ND_file.seekp((max_block_number + 1) * DISK_BLOCK_SIZE);
      while(count > 0){
         block_number = tmp_stack.top();
         tmp_stack.pop();
         ND_file.write((const char*)(&block_number), sizeof(unsigned int));
         count--;
      }
   }
}

#endif
