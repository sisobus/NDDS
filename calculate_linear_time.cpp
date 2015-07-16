#include <cstdio>
#include <sys/stat.h>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <dirent.h>
#include <vector>
#include <cstring>
#include <map>
#include <string>
#include <algorithm>
#include <glob.h>
using namespace std;
bool isDataFile(const string& s) {
    char t[1024]={};
    strcpy(t,s.c_str());
    string extension = "";
    for ( char *p = strtok(t,".") ; p ; p=strtok(NULL,".") )
        extension = string(p);
    return extension == "txt";
}   
vector<string> getFileNamesFromPath(const string& pattern) {
    glob_t glob_result;
    glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
    vector<string> files;
    for ( unsigned int i = 0 ; i < glob_result.gl_pathc ; i++ ) 
        files.push_back(string(glob_result.gl_pathv[i]));
    globfree(&glob_result);
    return files;
}
void printVector(vector<string>& v){
    for ( int i = 0 ; i < (int)v.size() ; i++ ) 
        printf("%s\n",v[i].c_str());
}
vector<string> getStringSplit(const string& ts,const string& tpattern) {
    char s[1024]={},pattern[1024]={};
    strcpy(s,ts.c_str());
    strcpy(pattern,tpattern.c_str());
    vector<string> ret;
    for ( char *p=strtok(s,pattern) ; p ; p=strtok(NULL,pattern) )
        ret.push_back(string(p));
    return ret;
}
map<string,string> getOptionsFromFileName(const string& filename) {
    string onlyFileName = getStringSplit(filename,".")[0];
    onlyFileName        = getStringSplit(onlyFileName,"/")[1];
    vector<string> splitedFileName = getStringSplit(onlyFileName,"_");
    map<string,string> ret;
    ret["numberOfData"] = splitedFileName[1];
    ret["numberOfDimension"] = splitedFileName[2];
    ret["distribution"] = splitedFileName[3];
    ret["numberOfAlphabet"] = splitedFileName[4];
    return ret;
}
void printMap(map<string,string>& mp) {
    for ( map<string,string>::iterator it=mp.begin();it!=mp.end();it++ )
        printf("%s %s\n",it->first.c_str(),it->second.c_str());
}
vector<string> getDataFromDataFile(const string& filename,map<string,string>& options) {
    vector<string> ret;
    char s[1024]={};
    FILE *fp = fopen(filename.c_str(),"r");
    fscanf(fp,"%[^\n]\n",s);
    int size,dimension;
    sscanf(options["numberOfData"].c_str(),"%d",&size);
    sscanf(options["numberOfDimension"].c_str(),"%d",&dimension);
    for ( int i = 0 ; i < size ; i++ ) {
        fscanf(fp,"%[^\n]\n",s);
        string onlyData = getStringSplit(string(s),":")[1];
        vector<string> t = getStringSplit(onlyData,",");
        string tt = "";
        for ( int j = 0 ; j < (int)t.size() ; j++ ) 
            tt += t[j];
        ret.push_back(tt);
    }
    fclose(fp);
    return ret;
}
vector<string> getQueryFromQueryFile(const string& filename,map<string,string>& options) {
    vector<string> ret;
    FILE *fp = fopen(filename.c_str(),"r");
    int dimension;
    sscanf(options["numberOfDimension"].c_str(),"%d",&dimension);
    for ( int i = 0 ; i < 100 ; i++ ) {
        string s = "";
        for ( int j = 0 ; j < dimension ; j++ ) {
            char t[11]={};
            fscanf(fp,"%s",t);
            s += t;
        }
        ret.push_back(s);
    }
    fclose(fp);
    return ret;
}
int getHammingDistance(const string& s1,const string& s2) {
    assert(s1.length()==s2.length());
    int ret = 0;
    for ( int i = 0 ; i < (int)s1.length() ; i++ ) 
        ret += (s1[i] != s2[i]);
    return ret;
}
inline bool isFileExists(const string& fileName) {
    struct stat buffer;
    return (stat(fileName.c_str(),&buffer) == 0);
}
void rangeQuery(const vector<string>& data,const vector<string>& query,map<string,string> options) {
    string resultFileName = "./result/linear_"+options["numberOfData"]+"_"+options["numberOfDimension"]+"_"+options["distribution"]+"_"+options["numberOfAlphabet"]+".txt";
    if ( isFileExists(resultFileName) ) {
        printf("%s is exists!\n",resultFileName.c_str());
        return ;
    }
    float rqTime = 0.0;
    for ( int i = 0 ; i < (int)query.size() ; i++ ) {
        printf("query#%d\n",i);
        clock_t start = clock();
        int ans = 0;
        int radius = 1;
        for ( int j = 0 ; j < (int)data.size() ; j++ ) {
            int dist = getHammingDistance(query[i],data[j]);
            if ( dist <= radius ) ans++;
        }
        clock_t end = clock();
        float curTime = (end-start)/float(CLOCKS_PER_SEC);
        rqTime += curTime;
    }
    rqTime /= (int)query.size();
    FILE *fp = fopen(resultFileName.c_str(),"w");
    printf("%lf\n",rqTime);
    fprintf(fp,"QueryTime:%lf\n",rqTime);
    fclose(fp);
}
int main() {
    vector<string> filenames = getFileNamesFromPath("./data/*.txt");
    for ( int i = 0 ; i < (int)filenames.size() ; i++ ) {
        printf("%s\n",filenames[i].c_str());
        map<string,string> options = getOptionsFromFileName(filenames[i]);
        string queryFileName = "./query/query_"+options["numberOfData"]+"_"+options["numberOfDimension"]+"_"+options["distribution"]+"_"+options["numberOfAlphabet"]+".txt";
        vector<string> data = getDataFromDataFile(filenames[i],options);
        vector<string> query = getQueryFromQueryFile(queryFileName,options);
        rangeQuery(data,query,options);
    }
    return 0;
}
