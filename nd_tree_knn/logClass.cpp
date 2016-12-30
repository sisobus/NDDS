#include "logClass.h"
#include <iomanip>
logClass::logClass(void)
{
	logFlag=true;
}

logClass::~logClass(void)
{
}

bool logClass::log2File(string infoStr , string fileName)
{
	//return 1;
	if(!logFlag)
		return true;
	ofstream fout;

	fout.open(fileName.c_str(),ios::out|ios::app);
	if(fout.fail())
	{
		cout<<"can't open "<<fileName<<endl;
	}

	fout<<infoStr;
	//cout<<"log: "<<infoStr;
	//cout<<"log: "<<infoStr;
	fout.close();

	return true;
}


bool logClass::log2File(float infoFloat , string fileName)
{
		//return 1;
		if(!logFlag)
		return true;

	ofstream fout;

	fout.open(fileName.c_str(),ios::out|ios::app);
	if(fout.fail())
	{
		cout<<"can't open "<<fileName<<endl;
	}

	//if(infoInt==3)
	//	cout<<"";

	fout<<fixed<<showpoint<<setprecision(10)<<infoFloat;
	//cout<<"log: "<<infoInt;
	//cout<<"log: "<<infoInt;

	fout.close();

	return true;
}

bool logClass::log2File(unsigned int infoInt , string fileName)
{
		//return 1;
		if(!logFlag)
		return true;

	ofstream fout;

	fout.open(fileName.c_str(),ios::out|ios::app);
	if(fout.fail())
	{
		cout<<"can't open "<<fileName<<endl;
	}

	//if(infoInt==3)
	//	cout<<"";

	fout<<setw(4)<<infoInt;
	//cout<<"log: "<<infoInt;
	//cout<<"log: "<<infoInt;

	fout.close();

	return true;
}

bool logClass::log2File( int infoInt , string fileName)
{
		//return 1;
		if(!logFlag)
		return true;

	ofstream fout;

	fout.open(fileName.c_str(),ios::out|ios::app);
	if(fout.fail())
	{
		cout<<"can't open "<<fileName<<endl;
	}

	//if(infoInt==3)
	//	cout<<"";

	fout<<setw(4)<<infoInt;
	//cout<<"log: "<<infoInt;
	//cout<<"log: "<<infoInt;

	fout.close();

	return true;
}

int logClass::string_to_int(string s){

   istringstream instr(s);
   int n;
   instr >> n;
   return n;
}


string logClass::int_to_string(int i){
   
	ostringstream outstr;
	outstr << i;
	return outstr.str();
}

void logClass::setlogFlag(bool value)
{
	logFlag=value;

}

void logClass::clearLogs(string fileName)
{
	ofstream fout;

	fout.open(fileName.c_str(),ios::out|ios::trunc );
	if(fout.fail())
	{
		cout<<"can't open "<<fileName<<endl;
	}

	fout.close();
}
