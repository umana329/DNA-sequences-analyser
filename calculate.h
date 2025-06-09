#pragma once
#include<string>
#include<algorithm>
#include<vector>
using namespace std;
double GCpercentMin(string dna);
double GCpercentMax(string dna);
int matchScore(const char&,const char&);
pair<int, pair<string,string>>needlemanWunsch(const string&, const string& ,int=-1);
double turnIntoDouble(const string);
void normalization(vector<double>&);
void checkT(const double&,const int&);
void csvResult0point2no();
void csvResult0point2yes();
void csvResult0point1();
void csvResult0point05();
void csvResult0point02();
void csvResult0point01();
void csvResult0point002();
void csvResult0point001();