#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<iomanip>
#include<algorithm>
#include"calculate.h"
using namespace std;
void readfasta(const string&, vector<string>&, vector<string>&);

int main() {
vector<string> seq;
vector<string> seqInf;
string filename="C:/fasta/",filename1;
cout << "输入要读取的文件名(记得加上\".fa\")：";
cin >> filename1;
filename += filename1;
readfasta(filename, seq, seqInf);
    for (size_t i = 0; i < seqInf.size(); i++) {
	    cout << "Sequence Info:" << seqInf[i] << endl;
	    cout << "Sequence:" << seq[i] << endl;
		cout << "Length:" << seq[i].size() << endl;
		cout << "Proportion of GC:" << setprecision(3) << GCpercent(seq[i]) * 100 << "%\n" << endl;
	}
	string s1 = "GATTACA", s2 = "GCATGCC";//用于测试的序列
	auto result = needlemanWunsch(s1, s2);
	cout << "score " << result.first << endl << "result " << endl << result.second << endl;
return 0;
}

	void readfasta(const string& filename, vector<string>& seq, vector<string>& inf) {
		ifstream file(filename);
		string line, currentSeq = "", currentInf = "";
		if (file.is_open()) {
			while (std::getline(file, line)) {
				if (line.empty()) {//skip empty lines
					continue;
				}
				if (line[0] == '>') {//Info lines start with >
					if (!currentSeq.empty()) {
						seq.push_back(currentSeq);
						inf.push_back(currentInf);
						currentInf = currentSeq = "";
					}
					currentInf = line.substr(1);
				}
				else {//Sequqnce lines
					currentSeq += line;
				}
			}
			if (!currentSeq.empty()) {
				seq.push_back(currentSeq);
				inf.push_back(currentInf);
			}
			file.close();
		}
		else {
			cerr << "Failed to open the file!" << endl;
		}
	}