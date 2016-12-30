#include "utility.h"

int string_to_int(string s){
   istringstream instr(s);
   int n;
   instr >> n;
   return n;
}

string int_to_string(int i){
   ostringstream outstr;
   outstr << i;
   return outstr.str();
}


string lpad(string s, unsigned int str_len, string pad_s){
   string result = s;
   if(s.size() < str_len){
      int difference = str_len - s.size();
      for(int i = 0; i < difference; i++)
         result = pad_s + result;
   }
   return result;
}


//supports only ascending order
void NDT_stable_sort(int count, double* weight_list, int* sorted_index_list, bool asc){
	// populate
	//for(int i = 0; i < count; i++)
	//   sorted_index_list[i] = i;
		std::vector<std::pair<double, int> > myVector;

		for(int i = 0; i < count; i++)
			myVector.push_back(std::pair<double, int>(weight_list[i], i));


	if(asc==false)//descending order
	{
		for(int i =count-1;i>=0;i--)
		{

			int tempIndex=i;
			for(int j =i-1;j>=0;j--)
			{
				if(myVector.at(j).first <myVector.at(tempIndex).first)
					tempIndex=j;
			}
			double tempWeight;
			std::pair<double, int> tempItem;
			tempItem=myVector.at(tempIndex);
			myVector.at(tempIndex)=myVector.at(i);
			myVector.at(i)=tempItem;
			sorted_index_list[i]=tempItem.second;
		}

		//int x = 0, y=count-1;

		//while(x<y)
		//{
		// double d;
		// int t;
		// d=weight_list[x];
		// weight_list[x]=weight_list[y];
		// weight_list[y]=d;
		// t=sorted_index_list[x];
		// sorted_index_list[x]=sorted_index_list[y];
		// sorted_index_list[y]=t;

		// x++;
		// y--;
		//}
	}
	else
	{

		std::stable_sort(myVector.begin(), myVector.end());

		for(int i = 0; i < count; i++)
		{
			//if(asc)
			//{
			weight_list[i]=myVector.at(i).first;
			sorted_index_list[i]=myVector.at(i).second;
			//}
			//else
			//{
			// weight_list[i]=myVector.at(count-1-i).first;
			// sorted_index_list[i]=myVector.at(count-1-i).second;
			//}
		}

	}
}


//not stable 
void NDT_sort(int count, double* weight_list, int* sorted_index_list, bool asc){
   // populate
   for(int i = 0; i < count; i++)
      sorted_index_list[i] = i;
   
	int position, largest, current;
   int tmp_index;
   double tmp_weight;
   if(asc){ // ascending order
      // Selection sort
      for(position = count - 1; position > 0; position--){
         // pick the index to the largest element
	      largest = 0;
	      for(current = 1; current <= position; current++)
		      if (weight_list[largest] < weight_list[current])largest = current;

         // swap
         tmp_weight = weight_list[largest];
         weight_list[largest] = weight_list[position];
         weight_list[position] = tmp_weight;
         tmp_index = sorted_index_list[largest];
         sorted_index_list[largest] = sorted_index_list[position];
         sorted_index_list[position] = tmp_index;
	   }
   }else{ // descending order
      // Selection sort
      for(position = 0; position < count - 1; position++){
         // pick the index to the largest element
	      largest = count - 1;
	      for(current = position; current <= count - 2; current++)
		      if(weight_list[largest] < weight_list[current])largest = current;

         // swap
         tmp_weight = weight_list[largest];
         weight_list[largest] = weight_list[position];
         weight_list[position] = tmp_weight;
         tmp_index = sorted_index_list[largest];
         sorted_index_list[largest] = sorted_index_list[position];
         sorted_index_list[position] = tmp_index;
	   }
   }
}

void NDT_sort(int count, unsigned int* weight_list){
	int position, largest, current;
   unsigned int tmp_weight;

   // Selection sort
   for(position = count - 1; position > 0; position--){
      // pick the index to the largest element
	   largest = 0;
	   for(current = 1; current <= position; current++)
		   if (weight_list[largest] < weight_list[current])largest = current;

      // swap
      tmp_weight = weight_list[largest];
      weight_list[largest] = weight_list[position];
      weight_list[position] = tmp_weight;
	}
}

void NDT_sort(unsigned int* a, int lo, int hi)
/* sort a[lo..hi] */
 { int left, right;
   unsigned int median, temp;

   if( hi > lo ) /* i.e. at least 2 elements, then */
    { left=lo; right=hi;
      median=a[lo];  /* NB. just an estimate! */

      while(right >= left) /* partition a[lo..hi] */
      /* a[lo..left-1]<=median and a[right+1..hi]>=median */
       { while(a[left] < median) left++;
         /* a[left] >= median */

         while(a[right] > median) right--;
         /* a[left] >= median >= a[right] */

         if(left > right) break;

         //swap:
         temp=a[left]; a[left]=a[right]; a[right]=temp;
         left++; right--;
       }
      /* a[lo..left-1]<=median and a[right+1..hi]>=median
         and left > right */

      NDT_sort(a, lo, right);// divide and conquer
      NDT_sort(a, left,  hi);
    }

}

void string2binary(const string &origString, unsigned char* const ptr)
{
	int strSizeInBytes = ceil((float)origString.length()/8);


	for(int i =0;i<strSizeInBytes;i++)
	{
		ptr[i]=0x00;

		for (int j =0;j<8;j++)
		{
			if(origString.at(i*8+j)=='1')
				ptr[i]^=bit_mask[7-j];
		}	
	
	}

}

void binary2string(const unsigned char* const ptr, int lengthInBytes, string & resultString)
{
	resultString="";
	
	for(int i=0;i<lengthInBytes;i++)
	{
		for (int j =7;j>=0;j--)
		{
			if(ptr[i]&bit_mask[j])
				resultString+='1';
			else
				resultString+='0';
		}
	}

}



vector<int> Knapsack_recursive_forNDTree(
	int n, 
	int available_size, //availabe knapsack capacity, initialy equals knapsack_capacity = block size - overhead
const int * const item_size,
const int * const item_value,	
const int 	all_entries_size,
const float knapsack_min_capacity,
const int knapsack_capacity)
{
	//static int debug_level = 0;
	//debug_level++;

	vector<int>result;
	assert(n>=0);
	assert(available_size>=0);
	if((n==0)||(available_size==0))
	{
		//cout<<"debug_level: "<<debug_level<<endl;
		//cout<<" Knapsack_recursive "<<n<<","<<available_size<<" "<<0<<endl;

		//debug_level--;
		return result;
	}
	else
	{
		vector<int> ans1_vec,ans2_vec;
		ans1_vec=Knapsack_recursive_forNDTree(n-1,available_size,item_size,item_value,all_entries_size,knapsack_min_capacity,knapsack_capacity);

		if(available_size>=item_size[n-1])
		{
			ans2_vec = Knapsack_recursive_forNDTree(n-1,available_size-item_size[n-1],item_size,item_value,all_entries_size,knapsack_min_capacity,knapsack_capacity);
			ans2_vec.push_back(n-1);
		}
		else
			ans2_vec.clear();

		int ans1_value=0,ans2_value=0;
		for(int i =0;i<ans1_vec.size();i++)
			ans1_value+=item_value[ans1_vec.at(i)];
		for(int i =0;i<ans2_vec.size();i++)
			ans2_value+=item_value[ans2_vec.at(i)];


		float ans1_diskSize=0,ans2_diskSize=0;

		for(int i =0;i<ans1_vec.size();i++)
			ans1_diskSize+=item_size[ans1_vec.at(i)];

		ans1_diskSize+=knapsack_capacity - available_size;//plus disk sizes taken by parents

		for(int i =0;i<ans2_vec.size();i++)
			ans2_diskSize+=item_size[ans2_vec.at(i)];

		ans2_diskSize+=knapsack_capacity - available_size;//plus disk sizes taken by parents


		float ans1_diskSize_counterPart = all_entries_size -ans1_diskSize;
		float ans2_diskSize_counterPart = all_entries_size -ans2_diskSize;

		bool ans1_ok=(ans1_diskSize>=knapsack_min_capacity)&&(ans1_diskSize_counterPart>=knapsack_min_capacity)&&(ans1_diskSize_counterPart<=knapsack_capacity);
		bool ans2_ok=(ans2_diskSize>=knapsack_min_capacity)&&(ans2_diskSize_counterPart>=knapsack_min_capacity)&&(ans2_diskSize_counterPart<=knapsack_capacity);

		if(!ans1_ok)
		{
			ans1_value=0;
			ans1_vec.clear();
		}
		if(!ans2_ok)
		{
			ans2_value=0;
			ans2_vec.clear();
		}
		//cout<<"debug_level: "<<debug_level<<endl;
		//cout<<" Knapsack_recursive "<<n<<","<<available_size<<" "<<(ans1_value>ans2_value?ans1_value:ans2_value)<<endl;
		//debug_level--;
		return ans1_value>ans2_value?ans1_vec:ans2_vec;
	}
}








