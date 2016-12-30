#include <cstdio>
#include <algorithm>
#include <string>
#include <vector>
#include <cmath>
using namespace std;
double f(vector<double>& v1,vector<double>& v2) {
    double ret = 0;
    for ( int i = 0 ; i < (int)v1.size() ; i++ ) 
        ret += (v1[i]-v2[i])*(v1[i]-v2[i]);
    return sqrt(ret);
}
void printAverage(vector<vector<double> >& v) {
    double average = 0;
    int cnt = 0;
    double sqav = 0;
    for ( int i = 0 ; i < (int)v.size()/100 ; i++ ){ 
        double mn = 987654321;
        for ( int j = 0 ; j < (int)v.size()/100 ; j++ ) 
            if ( i != j ) {
                mn = min(mn,f(v[i],v[j]));
            }
        cnt ++;
        average += mn;
        sqav += mn*mn;
    }
    average /= cnt;
    double sd = sqav/cnt - average*average;
    printf("%lf %lf\n",average,sd);
}
int main() {
    int d;
    scanf("%d",&d);
    char s[128];
    sprintf(s,"data_100000_%d_u_10.txt",d);
    FILE* fp = fopen(s,"r");
    vector<vector<double> > v;
    for ( int i = 0 ; i < 100000 ; i++ ) {
        vector<double> cur;
        double t;
        for ( int j = 0 ; j < d ; j++ ) {
            fscanf(fp,"%lf",&t);
            cur.push_back(t);
        }
        v.push_back(cur);
    }
    fclose(fp);
    printAverage(v);
    
    return 0;
}
