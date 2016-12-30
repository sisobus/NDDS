#include <cstdio>
#include <set>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
using namespace std;
#define numberOfDimension 10
#define numberOfVantagePoint 5
#define numberOfAlphabet 4 
vector<string> list;
string alphabet;
void createAllList(string& now,int pos) {
    if ( pos == numberOfDimension ) {
        list.push_back(now);
        return;
    }
    for ( int i = 0 ; i < (int)alphabet.size() ; i++ ) {
        now.push_back(alphabet[i]);
        createAllList(now,pos+1);
        now.pop_back();
    }
}
void initializing() {
    list.clear();
    alphabet = "";
    for ( int i = 0 ; i < numberOfAlphabet ; i++ ) 
        alphabet += ('A'+i);
    string dummy = "";
    createAllList(dummy,0);
}
const double dabs(double a){return a>=0?a:-a;}
double calculateDistanceVariance(vector<string>& vps) {
    double squareSum = 0;
    double average = 0;
    int n = 0;
    for ( int i = 0 ; i < (int)vps.size() ; i++ ) 
        for ( int j = 0 ; j < (int)vps.size() ; j++ ) 
            if ( i != j ) {
                int now = 0;
                for ( int k = 0 ; k < (int)vps[i].length() ; k++ ) 
                    now += (vps[i][k] != vps[j][k]);
                squareSum += now*now;
                average += now;
                n++;
            }
    average /= n; 
    return squareSum / n - average*average;
}
struct Tuple {
    double data;
    int id;
    Tuple(){}
    Tuple(double _data,int _id):
        data (_data), id (_id){}
};
double calculateMany(vector<string>& vps) {
    set<int> s;
    for ( int i = 0 ; i < (int)vps.size() ; i++ ) 
        for ( int j = 0 ; j < (int)vps.size() ; j++ ) 
            if ( i != j ) {
                int now = 0;
                for ( int k = 0 ; k < (int)vps[i].length() ; k++ ) 
                    now += (vps[i][k] != vps[j][k]);
                s.insert(now);
            }
    return (double)s.size();
}
vector<string> createVantagePointsWithManyAlgorithm() {
    vector<string> vps;
    vps.push_back(list[0]);
    while ( (int)vps.size() < numberOfVantagePoint ) {
        Tuple t(-1,-1);
        for ( int i = 1 ; i < (int)list.size() ; i++ ) {
            bool pass = false;
            for ( int j = 0 ; j < (int)vps.size() ; j++ ) 
                if ( vps[j] == list[i] ) pass = true;
            if ( pass ) continue;
            vps.push_back(list[i]);
            double cur = calculateMany(vps);
            if ( cur > t.data ) 
                t = Tuple(cur,i);
            vps.pop_back();
        }
        vps.push_back(list[t.id]);
    }
    return vps;
}

vector<string> createVantagePointsWithGreedyAlgorithm() {
    vector<string> vps;
    vps.push_back(list[0]);
    while ( (int)vps.size() < numberOfVantagePoint ) {
        Tuple t(-1,-1);
        for ( int i = 1 ; i < (int)list.size() ; i++ ) {
            bool pass = false;
            for ( int j = 0 ; j < (int)vps.size() ; j++ ) 
                if ( vps[j] == list[i] ) pass = true;
            if ( pass ) continue;
            vps.push_back(list[i]);
            double cur = calculateDistanceVariance(vps);
            if ( cur > t.data ) 
                t = Tuple(cur,i);
            vps.pop_back();
        }
        vps.push_back(list[t.id]);
    }
    return vps;
}
int calcuateHammingDistance(string& s1,string& s2) {
    int ret = 0;
    for ( int i = 0 ; i < (int)s1.length() ; i++ ) 
        ret += (s1[i] != s2[i]);
    return ret;
}
void printAllPairDistance(vector<string>& vps) {
    for ( int i = 0 ; i < (int)vps.size() ; i++ ) {
        for ( int j = 0 ; j < (int)vps.size() ; j++ ) {
            printf("(%d,%d,%d)\n",i+1,j+1,calcuateHammingDistance(vps[i],vps[j]));
        }
    }
}
int main() {
    initializing();
    //vector<string> vps = createVantagePointsWithGreedyAlgorithm();
    vector<string> vps = createVantagePointsWithManyAlgorithm();
    for ( int i = 0 ; i < (int)vps.size() ; i++ ) 
        printf("%s\n",vps[i].c_str());
    printf("%lf\n",calculateDistanceVariance(vps));
    printAllPairDistance(vps);
    return 0;
}
