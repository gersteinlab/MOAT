#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <string>
#include <errno.h>
#include <sys/stat.h>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

using namespace std;

bool cmpIntervals (vector<string> a, vector<string> b);
vector<vector<string> > permute_variants (int varcount, vector<string> region);
vector<vector<string> > intersecting_intervals (vector<vector<string> > &intervals, vector<string> &region);
pair<vector<vector<string> >, unsigned int> intersecting_intervals_w (vector<vector<string> > &intervals, vector<string> &region);
pair<unsigned int,unsigned int> intersecting_variants (vector<vector<string> > const &intervals, vector<string> region, unsigned int left_pointer);
pair<unsigned int,unsigned int> intersecting_variantsV (vector<vector<string> > const &intervals, vector<string> region, unsigned int left_pointer);
vector<vector<string> > subtract_intervals (vector<string> region, vector<vector<string> > intervals);
string exec (const char* cmd);
