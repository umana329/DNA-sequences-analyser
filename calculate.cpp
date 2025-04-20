#include"calculate.h"
#include<string>
#include<iostream>
#include<algorithm>
#include<vector>
using namespace std;
double GCpercent(string dna) {//calculate the proprotion of GC
	double gc = 0;
	string::iterator it = dna.begin();
	while (it != dna.end()) {
		if (*it == 'C' || *it == 'G') {
			gc++;
		}
		it++;
	}
	return gc / dna.size();
}
int matchScore(char a, char b) {
	return (a == b) ? 1 : -1;
}
pair<int, string>needlemanWunsch(const string& seq1, const string& seq2, int gapPenalty) {//gapPenalty默认为-1
	int m = seq1.size(), n = seq2.size();
	vector<vector<int>>matrix(m + 1, vector<int>(n + 1, 0));//填充得分矩阵的首行、首列
	for (int i = 1; i <= m; i++) {
		matrix[i][0] = i * gapPenalty;
	}
	for (int i = 1; i <= n; i++) {
		matrix[0][i] = i * gapPenalty;
	}
	for (int i = 1; i <= m; i++) {//填充矩阵中间的部分
		for (int j = 1; j <= n; j++) {
			int match = matrix[i - 1][j - 1] + matchScore(seq1[i - 1], seq2[j - 1]);
			int insert = matrix[i][j - 1] + gapPenalty;
			int Delete = matrix[i - 1][j] + gapPenalty;
			matrix[i][j] = max({ match,insert,Delete });
		}
	}
	string align1, align2;
	int i = m, j = n;
	while (i > 0 || j > 0) {//回溯寻找最佳路径
		if (i > 0 && j > 0 && matrix[i][j] == matrix[i - 1][j - 1] + matchScore(seq1[i - 1], seq2[j - 1])) {
			align1 = seq1[i - 1] + align1;
			align2 = seq2[j - 1] + align2;
			i--, j--;
		}
		else if (i > 0 && matrix[i][j] == matrix[i - 1][j] + gapPenalty) {
			align1 = seq1[i - 1] + align1;
			align2 = '-' + align2;
			i--;
		}
		else {
			align1 = '-' + align1;
			align2 = seq2[j - 1] + align2;
			j--;
		}
	}
	return make_pair(matrix[m][n], align1 + "\n" + align2);
}