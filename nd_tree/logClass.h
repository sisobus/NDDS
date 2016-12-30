#ifndef LOGCLASS
#define LOGCLASS

#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <sstream>
using namespace std;

class logClass
{


public:
	logClass(void);
	
public:
	~logClass(void);
public:
	bool logFlag;
	bool log2File(string infoStr  , string fileName="defaultLogFile.txt");
	bool log2File( int infoInt  , string fileName="defaultLogFile.txt");
	bool log2File(unsigned int infoInt  , string fileName="defaultLogFile.txt");
	bool log2File(float infoFloat  , string fileName="defaultLogFile.txt");



	int string_to_int(string s);
	string int_to_string(int i);

	void setlogFlag(bool value);

	void clearLogs(string fileName="defaultLogFile.txt");
};


#endif
