#include "ND_tree.h"
#include <vector>
#include "logClass.h"

ofstream OutStream;

//#define MAX_COPIES_OF_DATA_FILE 50

//double dir_min_util = 0.0000000000003;
const double dir_min_util = ((nodeSplitType==ORIGINAL)||(enforce_minUtil_for_exhaustSplit==1))?0.3:0.0000000000003;
const double leaf_min_util = 0.3;

int tmp_DMBR_byte_lut[DIM][MAX_ALPHABET_SIZE]; // Given a dim and a letter, store the coresponding byte in DMBR
int tmp_DMBR_bit_lut[DIM][MAX_ALPHABET_SIZE]; // Given a dim and a letter, store the coresponding bit in DMBR

int debug_boxQ_leaf_accessed=0;
int debug_boxQ_leaf_hit=0;
int debug_boxQ_leaf_hit_peak=0;
vector<int> debug_boxQ_leaf_hit_for_all;

logClass logO;

int debug_height;

int duplicateDataPoints; //used in batchBuild(...), batchGrow(...)

void LocalDMBRInfoCalculation()
{
	int tmp_byte = 0;
	int i, j;
	for(i = 0; i < DIM; i++)
	{
		//tmp_DMBR_start_byte_lut[i] = tmp_byte;

		for(j = 0; j < A[i]; j++)
		{
			tmp_DMBR_byte_lut[i][j] = tmp_byte + j / BITS_PER_BYTE;
			tmp_DMBR_bit_lut[i][j] = j % BITS_PER_BYTE; ////letters must starts from 0
		}


		/** how many bytes occupied by this DBMR[i][j] **/
		tmp_byte += (A[i] - 1) / BITS_PER_BYTE + 1;

		//tmp_DMBR_end_byte_lut[i] = tmp_byte - 1;
	}



}

Dir_entry makeBoxQueryData(ifstream & query_file)
{
	
	//int maxDimAndContDim = TOTAL_DIM;//01/09/2007 commented out

	Dir_entry dirEntry;

	for(int j = 0; j < DMBR_SIZE; j++)
		dirEntry.DMBR[j]=0;


	string line;



	for(int j=0;j<DIM;j++)
	{

		//Node<DIM> tempNode;

		int byte_no,bit_no;

		getline(query_file, line);
		istringstream instr(line);

		int v;
		for(int t=0;t<boxSize;t++)
		{
			instr>>v;

			byte_no = tmp_DMBR_byte_lut[j][v];

			bit_no = tmp_DMBR_bit_lut[j][v];/** todo: here DMBR_bit_lut[j] should support enough letters **/
			dirEntry.DMBR[byte_no] |= MASKS[bit_no];
		}
	}


	for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
		getline(query_file, line);


	//consumes the rest lines holding float values
	for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
		getline(query_file, line);


	return dirEntry;
}

void batchGrow_with_duplicate(int skipSize,	int sizeGrowTo)
{


	ND_tree ndt;
	Leaf_entry new_data/*, query_data, db_data*/;

	Error_code result;

	string input="ndTree.dat";
	//string input="tree"+int_to_string(size)+".dat";
	//////////////////////////////
	//ndt.create_empty_tree(A, dir_min_util, leaf_min_util, input);


	ndt.read_existing_tree(input);


	string input_fn;
	if(RUNNING_ENVIRONMENT == UNIX)
		input_fn="sourceData"+int_to_string(TOTAL_DSC_VALUE)+ "+" +int_to_string(TOTAL_CNT_VALUE)+".txt";
	else
		input_fn="C:\\temp\\sourceData"+int_to_string(TOTAL_DSC_VALUE)+ "+" +int_to_string(TOTAL_CNT_VALUE)+".txt";


	string str_num_of_points;



	ifstream data_file;
	data_file.open(input_fn.c_str());
	if (data_file.fail())
	{
		cout<<"cant open file "<<input_fn<<endl;
		exit(1);

	}

	int number_of_io = 0;

	string line;


	int n;

	for(int i = 0; i < sizeGrowTo; i++)
	{

#ifdef LOG_VECTOR_INDEX
logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif

		if(RUNNING_ENVIRONMENT == WINDOWS)
			cout << i << endl;

#ifdef LOG_SPLITTING_VECTOR_INDEX
		vector_index = i;
#endif

		getline(data_file, line);

		if(i<skipSize)
			continue;

		debug_dataIndex = i;

		istringstream instr(line);

		for(int j = 0; j < DIM; j++)
		{
			instr >> n;
			new_data.key[j] = n;
		}

		//new_data.record = i;
		new_data.record = 1; 


		//if(USE_LINK_FOR_BOX_QUERY==1)
			result = ndt.insert_use_link(new_data, number_of_io);
		//else
		//	result = ndt.insert(new_data, number_of_io);


		if(result == duplicate_error)
		{
			//cout << "  Duplicate record" << endl;
			//distinctDataPoints--;
			duplicateDataPoints++;
		}


	}
	data_file.close();

	//finish = clock();
	//elapsed_time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "Time (in seconds): " << elapsed_time << endl;
	cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;

	cout<<"DistinctDataPoints data points indexed:"<<(sizeGrowTo-duplicateDataPoints )<<endl;
	cout<<"Total data points read:"<<sizeGrowTo <<endl;

	ndt.print_information( );
	//cout<<"------------------------------------"<<endl;
}





//this one only reads data file up to the give number of data points 
void batchBuild_with_duplicate(int size)
{
		logO.clearLogs();
	duplicateDataPoints=0;

	ND_tree ndt;
	Leaf_entry new_data/*, query_data, db_data*/;

	Error_code result;

	string input="ndTree.dat";
	//string input="tree"+int_to_string(size)+".dat";
	//////////////////////////////
	ndt.create_empty_tree(A, dir_min_util, leaf_min_util, input);


	ndt.read_existing_tree(input);


	//string input_fn="C:\\temp\\sourceData8+0.txt";
	//string input_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+0.txt";
	//string input_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+" +int_to_string(CNTDIM)+".txt";
	string input_fn;
	if(RUNNING_ENVIRONMENT == UNIX)
		input_fn="sourceData"+int_to_string(TOTAL_DSC_VALUE)+ "+" +int_to_string(TOTAL_CNT_VALUE)+".txt";
	else
		input_fn="C:\\temp\\sourceData"+int_to_string(TOTAL_DSC_VALUE)+ "+" +int_to_string(TOTAL_CNT_VALUE)+".txt";


	string str_num_of_points;
	int num_of_points=size;


	ifstream data_file;
	data_file.open(input_fn.c_str());
	if (data_file.fail())
	{
		cout<<"cant open file "<<input_fn<<endl;
		exit(1);

	}

	int number_of_io = 0;

	string line;


	int n=0;

	for(int i = 0; i < num_of_points; i++)
	{
#ifdef LOG_VECTOR_INDEX
logO.log2File("----------");logO.log2File(i);logO.log2File("\n");
#endif

		if(RUNNING_ENVIRONMENT == WINDOWS)
			cout << i << endl;

#ifdef LOG_SPLITTING_VECTOR_INDEX
		vector_index = i;
#endif

		if(i==82427)
			i=i;
debug_dataIndex = i;

		getline(data_file, line);
		istringstream instr(line);

		for(int j = 0; j < DIM; j++)
		{
			instr >> n;
			new_data.key[j] = n;
		}
		//new_data.record = i;
		new_data.record = 1; 


		//if(USE_LINK_FOR_BOX_QUERY==1)
			result = ndt.insert_use_link(new_data, number_of_io);
		//else
		//	result = ndt.insert(new_data, number_of_io);


		if(result == duplicate_error)
		{
			//cout << "  Duplicate record" << endl;
			//distinctDataPoints--;
			duplicateDataPoints++;
		}


	}
	data_file.close();

	//finish = clock();
	//elapsed_time = (double)(finish - start) / CLOCKS_PER_SEC;
	//cout << "Time (in seconds): " << elapsed_time << endl;
	cout<<"Duplicate data points encountered:"<<duplicateDataPoints<<endl;

	cout<<"DistinctDataPoints data points indexed:"<<(size-duplicateDataPoints )<<endl;
	cout<<"Total data points read:"<<size <<endl;

	ndt.print_information( );
	//cout<<"------------------------------------"<<endl;
}



//this one tries its best to get the give number of identical data points out of the data file
//void batchBuild_best_effort(	int size)
//{
//
//
//	ND_tree ndt;
//	Leaf_entry new_data/*, query_data, db_data*/;
//
//	Error_code result;
//
//
//	string input="tree"+int_to_string(size)+".dat";
//	//////////////////////////////
//	ndt.create_empty_tree(A, dir_min_util, leaf_min_util, input);
//
//
//	ndt.read_existing_tree(input);
//
//	int distinctDataPoints=0;
//	//string input_fn="C:\\temp\\sourceData8+0.txt";
//	//string db_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+0.txt";
//	string input_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+" +int_to_string(CNTDIM)+".txt";
//	//string input_fn="/egr/scratch/chencha3/data2/sourceData"+int_to_string(DIM) + "+" +int_to_string(CNTDIM)+".txt";
//
//
//	string str_num_of_points;
//	int num_of_points=size;
//
//
//	ifstream data_file;
//	data_file.open(input_fn.c_str());
//
//	if (data_file.fail())
//	{
//		cout<<"cant open file "<<input_fn<<endl;
//		exit(1);
//
//	}
//
//
//
//	int number_of_io = 0;
//
//	string line;
//
//
//
//	int n,i;
//
//	for( i = 0, distinctDataPoints=0; distinctDataPoints < num_of_points && i<MAX_LINE_IN_SOURCE_FILE; i++)
//	{
//		cout << i << endl;
//		getline(data_file, line);
//		istringstream instr(line);
//		for(int j = 0; j < DIM; j++){
//			instr >> n;
//			new_data.key[j] = n;
//		}
//		new_data.record = i;
//		result = ndt.insert(new_data, number_of_io);
//		if(result == duplicate_error)
//		{
//			;
//			//cout << "  Duplicate record" << endl;
//			//distinctDataPoints--;
//			//num_of_points--;
//		}
//		else
//			distinctDataPoints++;
//	}
//
//
//
//	data_file.close();
//
//	//finish = clock();
//	//elapsed_time = (double)(finish - start) / CLOCKS_PER_SEC;
//	//cout << "Time (in seconds): " << elapsed_time << endl;
//	cout<<"Duplicate data points encountered:"<<i-distinctDataPoints<<endl;
//
//	cout<<"DistinctDataPoints data points indexed:"<<distinctDataPoints <<endl;
//	ndt.print_information( );
//	//cout<<"------------------------------------"<<endl;
//}
//
//void batchBoxQuery_V2(	int size ) //size is tree size, not query size
//{
//	string	input="tree"+int_to_string(size)+".dat";
//
//	ND_tree ndt;
//	ndt.read_existing_tree(input);
//
//
//
//
//	LocalDMBRInfoCalculation();
//
//
//
//
//	for(int discrete_DimToBeConstrained=DIM;discrete_DimToBeConstrained>=1;discrete_DimToBeConstrained--)
//	{
//
//		Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
//		int query_results_size;
//
//
//		Dir_entry boxQueryData;
//
//		//string query_fn="c:\\temp\\boxqueryAll.txt";
//		string query_fn="/egr/scratch/chencha3/data2/boxqueryAll.txt";
//
//		int num_of_points=200;
//
//		ifstream query_file;
//		query_file.open(query_fn.c_str());
//	if (query_file.fail())
//	{
//		cout<<"cant open file "<<query_fn<<endl;
//		exit(1);
//
//	}
//
//
//		debug_boxQ_leaf_hit_for_all.clear();
//
//		debug_boxQ_leaf_accessed=0;
//
//		int number_of_io ;
//		int total_number_of_io = 0;
//
////		int n;
//		int total_results_size=0;
//		for(int i = 0; i < num_of_points; i++)
//		{
//			debug_boxQ_leaf_hit_peak=0;
//			boxQueryData = makeBoxQueryData_V2(query_file,discrete_DimToBeConstrained);
//			ndt.box_query(boxQueryData, query_results, query_results_size, number_of_io);
//			total_number_of_io += number_of_io;
//			total_results_size += query_results_size;
//			//cout<<"for point "<<i<<" # of leaf node accessed: "<<debug_boxQ_leaf_accessed<<"and i/o: "<<number_of_io<<endl;
//		debug_boxQ_leaf_hit_for_all.push_back(debug_boxQ_leaf_hit_peak);
//
//		}
//
//
//		query_file.close();
//		cout<<"Discrete dim:"<<discrete_DimToBeConstrained<<" of "<<DIM<<" AvG matched data point found="<< static_cast<double>(total_results_size)/num_of_points;
//
//		cout << " AvG Disk I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
//		
//
//		//system("pause");
//	}
//
//}
//

void batchRangeQuery(	)
{

	Error_code result;

	string input="ndTree.dat";
	ND_tree ndt;
	Leaf_entry query_data;

	ndt.read_existing_tree(input);


				Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
				int query_results_size;

				string query_fn;
				////cout << endl << "Enter Query File Name:" << endl;
				////getline(cin, query_fn);

				//query_fn="c:\\temp\\queryData16+0.dat";
	if(RUNNING_ENVIRONMENT == UNIX)
		query_fn="rangequeryAll.txt";
	else
		query_fn="C:\\temp\\rangequeryAll.txt";

	//string str_range;

	//cout << endl << "Distance Range:" << endl;
	//getline(cin, str_range);
	//range = string_to_int(str_range);

	string str_num_of_points;
	int num_of_points;
	////cout << endl << "How many query records:" << endl;
	////getline(cin, str_num_of_points);
	////num_of_points = string_to_int(str_num_of_points);

	num_of_points=TOTAL_RANGE_QUERY_NUM;

	for(int range=1;range<RANGE_STOP_BEFORE;range+=RANGE_SIZE_STEP)
	{
		//start = clock();

		ifstream query_file;
		query_file.open(query_fn.c_str(),ios::in);
		if (query_file.fail())
		{
			cout<<"cant open file "<<query_fn<<endl;
			exit(1);

		}				

		int number_of_io = 0;
        int number_of_dist_computation = 0;
		int total_number_of_io = 0;
		int total_results_size = 0;
        int total_number_of_dist_computation = 0;

		string line;
		int n;
        clock_t start = clock();
        /*
         * sisobus modi
         */
        //num_of_points = 1;
		for(int i = 0; i < num_of_points; i++)
		{
			query_results_size = 0;

			//cout <<endl<<"range query for query point "<< i << endl;
			getline(query_file, line);

			if(query_file.bad()||query_file.fail() )

			{
				cout<<"reading query file error\n";
				exit(0);
			}

			istringstream instr(line);



			for(int j = 0; j < DIM; j++)
			{
				instr >> n;
				query_data.key[j] = n;
			}


			result = ndt.range_query_by_Hamming_dist(query_data, range, query_results, query_results_size, number_of_io, number_of_dist_computation);

			//cout << "  " << query_results_size << " " << number_of_io << endl;
			total_number_of_io += number_of_io;
			total_results_size += query_results_size;
            total_number_of_dist_computation += number_of_dist_computation;
		}
		query_file.close();

		clock_t finish = clock();
		double elapsed_time = (double)(finish - start) / CLOCKS_PER_SEC;
		cout << "Query Time (in seconds): " << elapsed_time / num_of_points << endl;
		cout << "For range: "<<range<<" range query I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
		cout << "Result size: " << static_cast<double>(total_results_size) / num_of_points << endl;
        cout << "distance computation: " << static_cast<double>(total_number_of_dist_computation) / num_of_points << endl;



        /*
         *sisobus modi
         */
        //break;

	}//end of for(int range=1;range<RANGE_STOP_BEFORE;++range)


	if(RUNNING_ENVIRONMENT == WINDOWS)
		system("pause");
}



void batchBoxQuery(/*	int size */) //size is tree size, not query size
{
	//string	input="tree"+int_to_string(size)+".dat";
	string input="ndTree.dat";
	ND_tree ndt;
	ndt.read_existing_tree(input);




	LocalDMBRInfoCalculation();




	for(boxSize=BOX_SIZE_START_AT;boxSize<BOX_SIZE_STOP_BEFORE;boxSize+=BOX_SIZE_STEP)
	{
		//logO.clearLogs();
		Leaf_entry query_results[QUERY_RESULTS_BUFFER_SIZE];
		int query_results_size;


		Dir_entry boxQueryData;

		string query_fn="c:\\temp\\boxqueryAll.txt";
		if(RUNNING_ENVIRONMENT == UNIX)
			query_fn="boxqueryAll.txt";
		else
			query_fn="C:\\temp\\boxqueryAll.txt";
		int num_of_points=TOTAL_BOX_QUERY_NUM;

		ifstream query_file;
		query_file.open(query_fn.c_str());
	if (query_file.fail())
	{
		cout<<"cant open file "<<query_fn<<endl;
		exit(1);

	}




		debug_boxQ_leaf_hit_for_all.clear();

		debug_boxQ_leaf_accessed=0;

		int number_of_io =0;
		int total_number_of_io = 0;

//		int n;
		int total_results_size=0;
		for(int i = 0; i < num_of_points; i++)
		{
			debug_boxQ_leaf_hit_peak=0;
			boxQueryData = makeBoxQueryData(query_file);
			ndt.box_query(boxQueryData, query_results, query_results_size, number_of_io);
			total_number_of_io += number_of_io;
			total_results_size += query_results_size;
			//cout<<"for point "<<i<<" # of leaf node accessed: "<<debug_boxQ_leaf_accessed<<"and i/o: "<<number_of_io<<endl;
			debug_boxQ_leaf_hit_for_all.push_back(debug_boxQ_leaf_hit_peak);

		}


		query_file.close();
		cout<<"boxSize= "<<boxSize<< " AvG boxquery I/O: " << static_cast<double>(total_number_of_io) / num_of_points << endl;
		cout<<" AvG matched data point found="<< static_cast<double>(total_results_size)/num_of_points<< endl; 


		cout << " AvG leaf node accessed: " << static_cast<double>(debug_boxQ_leaf_accessed) / num_of_points << endl;

		assert(debug_boxQ_leaf_hit_for_all.size()==num_of_points);
		int debug_tatol=0;
		for(unsigned int i=0;i<debug_boxQ_leaf_hit_for_all.size();i++)
		{
			
			debug_tatol+=debug_boxQ_leaf_hit_for_all.at(i);
		
		}

		cout << " AvG leaf node hit peak: " << static_cast<double>(debug_tatol) / num_of_points << endl;




	}





	
		if(RUNNING_ENVIRONMENT == WINDOWS)
			system("pause");
}


vector< vector<int> >  makeBoxQueryData_for_linearScan(ifstream & query_file)
{

	//int maxDimAndContDim = 16;//01/09/2007


	string line;

	vector< vector<int> >  boxQueryData;
	boxQueryData.resize(DIM);

	for(int j=0;j<DIM;j++)
	{

		getline(query_file, line);
		istringstream instr(line);
		boxQueryData[j].clear();

		int v;
		for(int t=0;t<(boxSize*A[j])/10;t++)
		{
			instr>>v;
			boxQueryData[j].push_back(v);

		}
	}


	for(int i=0;i<MAX_DIM_AND_CONTDIM-DIM;i++)
		getline(query_file, line);


	//consumes the rest lines holding float values
	for(int i=0;i<MAX_DIM_AND_CONTDIM;i++)
		getline(query_file, line);

	assert(boxQueryData.size()==DIM);
	assert(boxQueryData[0].size()==boxSize);

	return boxQueryData;
}








//result of this LinearScanBoxQuery should be 
void LinearScanBoxQuery(int dataNUM)
{
//	Leaf_entry query_results [QUERY_RESULTS_BUFFER_SIZE];


//	int query_results_size;


	int db_size=dataNUM;
	////string db_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+0.txt";
	//string db_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+" +int_to_string(CNTDIM)+".txt";
	//ifstream db_file;

	//db_file.open(db_fn.c_str());

	string query_fn="c:\\temp\\boxqueryAll.txt";


	ifstream query_file;
	query_file.open(query_fn.c_str());
	if (query_file.fail())
	{
		cout<<"cant open file "<<query_fn<<endl;
		exit(1);

	}




	for(boxSize=1;boxSize<=9;boxSize+=1)
	{


		query_file.seekg(0, ios::beg);


		int num_of_points=200;

		int totalNumOfMatches=0;

		for(int i = 0; i < num_of_points; i++) // for every query point
		{

			vector< vector<int> > matchedPoints;
			matchedPoints.clear();

			vector< vector<int> >  boxQueryData;
			boxQueryData = makeBoxQueryData_for_linearScan(query_file);



			vector<int> db_data;

			//string db_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+0.txt";
			string db_fn="C:\\temp\\sourceData"+int_to_string(DIM) + "+" +int_to_string(CNTDIM)+".txt";
			ifstream db_file;
			db_file.open(db_fn.c_str());

	if (db_file.fail())
	{
		cout<<"cant open file "<<db_fn<<endl;
		exit(1);

	}




			for(int j = 0; j < db_size; j++)
			{

				//cout<<"query data"<<i<<"db data"<<j<<endl;
				db_data.clear();
				string line;
				getline(db_file, line);
				istringstream dbstr(line);
				for(int k = 0; k < DIM; k++)
				{
					int n;
					dbstr >> n;
					db_data.push_back(n);
				}


				bool matchOnAllDim = true;		
				for(unsigned int s=0;(s<boxQueryData.size())&&matchOnAllDim ;s++)
				{
					bool matchOnOneDim=false;
					for(unsigned int x=0;x<boxQueryData[s].size();x++)	
						if(boxQueryData[s][x]==db_data[s])
							matchOnOneDim=true;

					if(!matchOnOneDim)
						matchOnAllDim=false;

				}



				if(matchOnAllDim)
				{
					bool findDuplicate=false;
					for(unsigned int i=0;i<matchedPoints.size();i++)
						if(matchedPoints.at(i)==db_data)
							findDuplicate=true;

					if(!findDuplicate)
					{
						matchedPoints.push_back(db_data);
						//for(int t1=0;t1<boxQueryData.size();t1++)
						//{
						//	for(int t2=0;t2<boxQueryData.at(t1).size();t2++)
						//		cout<<boxQueryData.at(t1).at(t2)<<" ";
						//	cout<<endl;
 
						//}	

						//for(int t2=0;t2<db_data.size();t2++)
						//	cout<<db_data.at(t2)<<": ";
						//cout<<endl;


					}
				}


			} 
			db_file.close();
			totalNumOfMatches+=(int)matchedPoints.size();

		}// end of 		for(int i = 0; i < num_of_points; i++) // for every query point




		cout<<"linear scan boxSize="<<boxSize<<endl;
		cout <<  "matched data # from linear scan:" << static_cast<double>(totalNumOfMatches) / num_of_points << endl;

		//system("pause");
	}//end of 	for(boxSize=1;boxSize<=9;boxSize+=1)
	query_file.close();



}


int main(int argc, char *argv[])
{
//assert(1==0);
	if(nodeSplitType==ORIGINAL)
	{
		assert(dir_min_util>0.1);
		TREE_TYPE = STATIC_TREE;
	}
	//else
	//{
	//	assert(dir_min_util<0.001);
	//
	//}
	//batchBuild_with_duplicate(22155);return 0;
	//batchGrow_with_duplicate(22155,2000000);return 0;

	batchBuild_with_duplicate(MAX_LINE_IN_SOURCE_FILE);
	//batchBoxQuery();
	batchRangeQuery();
	logO.log2File("1M done\n");
    /*
	batchGrow_with_duplicate(1000000,2000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("2M done\n");
	batchGrow_with_duplicate(2000000,3000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("3M done\n");
	batchGrow_with_duplicate(3000000,4000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("4M done\n");
	batchGrow_with_duplicate(4000000,5000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("5M done\n"); return 0;	
	batchGrow_with_duplicate(5000000,6000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("6M done\n");
	batchGrow_with_duplicate(6000000,7000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("7M done\n");
	batchGrow_with_duplicate(7000000,8000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("8M done\n");
	batchGrow_with_duplicate(8000000,9000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("9M done\n");
	batchGrow_with_duplicate(9000000,10000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("10M done\n");

	batchGrow_with_duplicate(10000000,11000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("11M done\n");
	batchGrow_with_duplicate(11000000,12000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("12M done\n");
	batchGrow_with_duplicate(12000000,13000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("13M done\n");
	batchGrow_with_duplicate(13000000,14000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("14M done\n");
	batchGrow_with_duplicate(14000000,15000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("15M done\n");	
	batchGrow_with_duplicate(15000000,16000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("16M done\n");
	batchGrow_with_duplicate(16000000,17000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("17M done\n");
	batchGrow_with_duplicate(17000000,18000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("18M done\n");
	batchGrow_with_duplicate(18000000,19000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("19M done\n");
	batchGrow_with_duplicate(19000000,20000000);
	batchBoxQuery();
	batchRangeQuery();
	logO.log2File("20M done\n");
    */
	return 0;

}
