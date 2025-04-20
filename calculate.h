#pragma once
#include<string>
#include<algorithm>
#include<vector>
using namespace std;
double GCpercent(string dna);
int matchScore(char,char);
pair<int, string>needlemanWunsch(const string&, const string& ,int=-1);