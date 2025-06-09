#include"calculate.h"
#include<string>
#include<iostream>
#include<algorithm>
#include<vector>
#include<sstream>
#include<iomanip>
#include<locale>
#include<codecvt>
#include<cmath>
#include<graphics.h>
#undef max
using namespace std;

double GCpercentMin(string dna) {//calculate the minimun proprotion of GC
	double gc = 0;
	string::iterator it = dna.begin();

	while (it != dna.end()) {//only counts chars that are certainly G or C

		if (*it == 'C' || *it == 'G'|| *it=='S') {
			gc++;
		}

		it++;
	}

	return gc / dna.size();

}


double GCpercentMax(string dna) {//calculate the maximun proprotion of GC

	double gc = 0;
	string::iterator it = dna.begin();

	while (it != dna.end()) {//counts every chars that are probably G or C

		if (*it == 'C' || *it == 'G' || *it == 'S' || *it == 'R' || *it == 'Y' || *it == 'K' || *it == 'M' || *it == 'B' || *it == 'D' || *it == 'H' || *it == 'V' || *it == 'N') {
			gc++;
		}

		it++;
	}

	return gc / dna.size();

}


int matchScore(const char& a,const char& b) {
	return (a == b) ? 1 : -1;
}


pair<int, pair<string,string>>needlemanWunsch(const string& seq1, const string& seq2, int gapPenalty) {//gapPenalty默认为-1
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

	return make_pair(matrix[m][n], make_pair(align1 , align2));
}


double turnIntoDouble(const string s) {//turn string into double

	if (s.find('.') != string::npos) {

		string inte, deci;
		istringstream ss(s);
		getline(ss, inte, '.');
		getline(ss, deci, '.');
		double result = 0;

		for (int i = 0; i < inte.length(); i++) {

			switch (inte[i]) {

			case '0':
				break;

			case '1':
				result += 1 * pow(10, s.size() - i - 1);
				break;

			case '2':
				result += 2 * pow(10, s.size() - i - 1);
				break;

			case '3':
				result += 3 * pow(10, s.size() - i - 1);
				break;

			case '4':
				result += 4 * pow(10, s.size() - i - 1);
				break;

			case '5':
				result += 5 * pow(10, s.size() - i - 1);
				break;

			case '6':
				result += 6 * pow(10, s.size() - i - 1);
				break;

			case '7':
				result += 7 * pow(10, s.size() - i - 1);
				break;

			case '8':
				result += 8 * pow(10, s.size() - i - 1);
				break;

			case '9':
				result += 9 * pow(10, s.size() - i - 1);
				break;

			default:
				break;

			}

		}

		for (int i = 0; i < deci.length(); i++) {

			switch (deci[i]) {

			case '0':
				break;

			case '1':
				result += 1 * pow(10, -i - 1);
				break;

			case '2':
				result += 2 * pow(10, -i - 1);
				break;

			case '3':
				result += 3 * pow(10, -i - 1);
				break;

			case '4':
				result += 4 * pow(10, -i - 1);
				break;

			case '5':
				result += 5 * pow(10, -i - 1);
				break;

			case '6':
				result += 6 * pow(10, -i - 1);
				break;

			case '7':
				result += 7 * pow(10, -i - 1);
				break;

			case '8':
				result += 8 * pow(10, -i - 1);
				break;

			case '9':
				result += 9 * pow(10, -i - 1);
				break;

			default:
				break;

			}

		}

		return result;
	}

	else {

		double result = 0;

		for (int i = 0; i < s.length(); i++) {

			switch (s[i]) {

			case '0':
				break;

			case '1':
				result += 1 * pow(10, s.size() - i - 1);
				break;

			case '2':
				result += 2 * pow(10, s.size() - i - 1);
				break;

			case '3':
				result += 3 * pow(10, s.size() - i - 1);
				break;

			case '4':
				result += 4 * pow(10, s.size() - i - 1);
				break;

			case '5':
				result += 5 * pow(10, s.size() - i - 1);
				break;

			case '6':
				result += 6 * pow(10, s.size() - i - 1);
				break;

			case '7':
				result += 7 * pow(10, s.size() - i - 1);
				break;

			case '8':
				result += 8 * pow(10, s.size() - i - 1);
				break;

			case '9':
				result += 9 * pow(10, s.size() - i - 1);
				break;

			default:
				break;

			}

		}
		return result;

	}

}


void normalization(vector<double>& record) {//归一化

	if (record.empty()) {
		return;
	}

	double Max = record[0], Min = record[0];

	for (int i = 1; i < record.size(); i++) {
		Max = max(Max, record[i]);
		Min = min(Min, record[i]);
	}

	for (int i = 0; i < record.size(); i++) {
		record[i] = (record[i] - Min) / (Max - Min);
	}

	return;

}


void checkT(const double& t,const int& n) {//t为计算得到的t值，n为自由度(样本量-1)
	const double a = abs(t);
	switch (n) {
	case 1:
		if (a<3.078) {
			csvResult0point2no();
		}
		else if (a < 6.314) {
			csvResult0point2yes();
		}
		else if (a < 12.706) {
			csvResult0point1();
		}
		else if (a < 31.821) {
			csvResult0point05();
		}
		else if (a < 63.656) {
			csvResult0point02();
		}
		else if (a < 318.289) {
			csvResult0point01();
		}
		else if (a < 636.578) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 2:
		if (a < 1.886) {
			csvResult0point2no();
		}
		else if (a < 2.92) {
			csvResult0point2yes();
		}
		else if (a < 4.303) {
			csvResult0point1();
		}
		else if (a < 6.965) {
			csvResult0point05();
		}
		else if (a < 9.925) {
			csvResult0point02();
		}
		else if (a < 22.328) {
			csvResult0point01();
		}
		else if (a < 31.6) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 3:
		if (a < 1.638) {
			csvResult0point2no();
		}
		else if (a < 2.353) {
			csvResult0point2yes();
		}
		else if (a < 3.182) {
			csvResult0point1();
		}
		else if (a < 4.541) {
			csvResult0point05();
		}
		else if (a < 5.841) {
			csvResult0point02();
		}
		else if (a < 10.214) {
			csvResult0point01();
		}
		else if (a < 12.924) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 4:
		if (a < 1.533) {
			csvResult0point2no();
		}
		else if (a < 2.132) {
			csvResult0point2yes();
		}
		else if (a < 2.776) {
			csvResult0point1();
		}
		else if (a < 3.747) {
			csvResult0point05();
		}
		else if (a < 4.604) {
			csvResult0point02();
		}
		else if (a < 7.173) {
			csvResult0point01();
		}
		else if (a < 8.61) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 5:
		if (a < 1.476) {
			csvResult0point2no();
		}
		else if (a < 2.015) {
			csvResult0point2yes();
		}
		else if (a < 2.571) {
			csvResult0point1();
		}
		else if (a < 3.365) {
			csvResult0point05();
		}
		else if (a < 4.032) {
			csvResult0point02();
		}
		else if (a < 5.894) {
			csvResult0point01();
		}
		else if (a < 6.869) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 6:
		if (a < 1.44) {
			csvResult0point2no();
		}
		else if (a < 1.943) {
			csvResult0point2yes();
		}
		else if (a < 2.447) {
			csvResult0point1();
		}
		else if (a < 3.143) {
			csvResult0point05();
		}
		else if (a < 3.707) {
			csvResult0point02();
		}
		else if (a < 5.208) {
			csvResult0point01();
		}
		else if (a < 5.959) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 7:
		if (a < 1.415) {
			csvResult0point2no();
		}
		else if (a < 1.895) {
			csvResult0point2yes();
		}
		else if (a < 2.365) {
			csvResult0point1();
		}
		else if (a < 2.998) {
			csvResult0point05();
		}
		else if (a < 3.499) {
			csvResult0point02();
		}
		else if (a < 4.785) {
			csvResult0point01();
		}
		else if (a < 5.408) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 8:
		if (a < 1.397) {
			csvResult0point2no();
		}
		else if (a < 1.86) {
			csvResult0point2yes();
		}
		else if (a < 2.306) {
			csvResult0point1();
		}
		else if (a < 2.896) {
			csvResult0point05();
		}
		else if (a < 3.355) {
			csvResult0point02();
		}
		else if (a < 4.501) {
			csvResult0point01();
		}
		else if (a < 5.041) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 9:
		if (a < 1.383) {
			csvResult0point2no();
		}
		else if (a < 1.833) {
			csvResult0point2yes();
		}
		else if (a < 2.262) {
			csvResult0point1();
		}
		else if (a < 2.821) {
			csvResult0point05();
		}
		else if (a < 3.25) {
			csvResult0point02();
		}
		else if (a < 4.297) {
			csvResult0point01();
		}
		else if (a < 4.781) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 10:
		if (a < 1.372) {
			csvResult0point2no();
		}
		else if (a < 1.812) {
			csvResult0point2yes();
		}
		else if (a < 2.228) {
			csvResult0point1();
		}
		else if (a < 2.764) {
			csvResult0point05();
		}
		else if (a < 3.169) {
			csvResult0point02();
		}
		else if (a < 4.144) {
			csvResult0point01();
		}
		else if (a < 4.587) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 11:
		if (a < 1.363) {
			csvResult0point2no();
		}
		else if (a < 1.796) {
			csvResult0point2yes();
		}
		else if (a < 2.201) {
			csvResult0point1();
		}
		else if (a < 2.718) {
			csvResult0point05();
		}
		else if (a < 3.106) {
			csvResult0point02();
		}
		else if (a < 4.025) {
			csvResult0point01();
		}
		else if (a < 4.437) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 12:
		if (a < 1.356) {
			csvResult0point2no();
		}
		else if (a < 1.782) {
			csvResult0point2yes();
		}
		else if (a < 2.179) {
			csvResult0point1();
		}
		else if (a < 2.681) {
			csvResult0point05();
		}
		else if (a < 3.055) {
			csvResult0point02();
		}
		else if (a < 3.93) {
			csvResult0point01();
		}
		else if (a < 4.318) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 13:
		if (a < 1.35) {
			csvResult0point2no();
		}
		else if (a < 1.771) {
			csvResult0point2yes();
		}
		else if (a < 2.16) {
			csvResult0point1();
		}
		else if (a < 2.65) {
			csvResult0point05();
		}
		else if (a < 3.012) {
			csvResult0point02();
		}
		else if (a < 3.852) {
			csvResult0point01();
		}
		else if (a < 4.221) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 14:
		if (a < 1.345) {
			csvResult0point2no();
		}
		else if (a < 1.761) {
			csvResult0point2yes();
		}
		else if (a < 2.145) {
			csvResult0point1();
		}
		else if (a < 2.624) {
			csvResult0point05();
		}
		else if (a < 2.977) {
			csvResult0point02();
		}
		else if (a < 3.787) {
			csvResult0point01();
		}
		else if (a < 4.14) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 15:
		if (a < 1.341) {
			csvResult0point2no();
		}
		else if (a < 1.753) {
			csvResult0point2yes();
		}
		else if (a < 2.131) {
			csvResult0point1();
		}
		else if (a < 2.602) {
			csvResult0point05();
		}
		else if (a < 2.947) {
			csvResult0point02();
		}
		else if (a < 3.733) {
			csvResult0point01();
		}
		else if (a < 4.073) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 16:
		if (a < 1.337) {
			csvResult0point2no();
		}
		else if (a < 1.746) {
			csvResult0point2yes();
		}
		else if (a < 2.12) {
			csvResult0point1();
		}
		else if (a < 2.583) {
			csvResult0point05();
		}
		else if (a < 2.921) {
			csvResult0point02();
		}
		else if (a < 3.686) {
			csvResult0point01();
		}
		else if (a < 4.015) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 17:
		if (a < 1.333) {
			csvResult0point2no();
		}
		else if (a < 1.74) {
			csvResult0point2yes();
		}
		else if (a < 2.11) {
			csvResult0point1();
		}
		else if (a < 2.567) {
			csvResult0point05();
		}
		else if (a < 2.898) {
			csvResult0point02();
		}
		else if (a < 3.646) {
			csvResult0point01();
		}
		else if (a < 3.965) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 18:
		if (a < 1.328) {
			csvResult0point2no();
		}
		else if (a < 1.729) {
			csvResult0point2yes();
		}
		else if (a < 2.093) {
			csvResult0point1();
		}
		else if (a < 2.539) {
			csvResult0point05();
		}
		else if (a < 2.861) {
			csvResult0point02();
		}
		else if (a < 3.579) {
			csvResult0point01();
		}
		else if (a < 3.883) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 19:
		if (a < 1.328) {
			csvResult0point2no();
		}
		else if (a < 1.729) {
			csvResult0point2yes();
		}
		else if (a < 2.093) {
			csvResult0point1();
		}
		else if (a < 2.539) {
			csvResult0point05();
		}
		else if (a < 2.861) {
			csvResult0point02();
		}
		else if (a < 3.579) {
			csvResult0point01();
		}
		else if (a < 3.883) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 20:
		if (a < 1.325) {
			csvResult0point2no();
		}
		else if (a < 1.725) {
			csvResult0point2yes();
		}
		else if (a < 2.086) {
			csvResult0point1();
		}
		else if (a < 2.528) {
			csvResult0point05();
		}
		else if (a < 2.845) {
			csvResult0point02();
		}
		else if (a < 3.552) {
			csvResult0point01();
		}
		else if (a < 3.85) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 21:
		if (a < 1.323) {
			csvResult0point2no();
		}
		else if (a < 1.721) {
			csvResult0point2yes();
		}
		else if (a < 2.08) {
			csvResult0point1();
		}
		else if (a < 2.518) {
			csvResult0point05();
		}
		else if (a < 2.831) {
			csvResult0point02();
		}
		else if (a < 3.527) {
			csvResult0point01();
		}
		else if (a < 3.819) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 22:
		if (a < 1.321) {
			csvResult0point2no();
		}
		else if (a < 1.717) {
			csvResult0point2yes();
		}
		else if (a < 2.074) {
			csvResult0point1();
		}
		else if (a < 2.508) {
			csvResult0point05();
		}
		else if (a < 2.819) {
			csvResult0point02();
		}
		else if (a < 3.505) {
			csvResult0point01();
		}
		else if (a < 3.792) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 23:
		if (a < 1.319) {
			csvResult0point2no();
		}
		else if (a < 1.714) {
			csvResult0point2yes();
		}
		else if (a < 2.069) {
			csvResult0point1();
		}
		else if (a < 2.5) {
			csvResult0point05();
		}
		else if (a < 2.807) {
			csvResult0point02();
		}
		else if (a < 3.485) {
			csvResult0point01();
		}
		else if (a < 3.768) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 24:
		if (a < 1.318) {
			csvResult0point2no();
		}
		else if (a < 1.711) {
			csvResult0point2yes();
		}
		else if (a < 2.064) {
			csvResult0point1();
		}
		else if (a < 2.492) {
			csvResult0point05();
		}
		else if (a < 2.797) {
			csvResult0point02();
		}
		else if (a < 3.467) {
			csvResult0point01();
		}
		else if (a < 3.745) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 25:
		if (a < 1.316) {
			csvResult0point2no();
		}
		else if (a < 1.708) {
			csvResult0point2yes();
		}
		else if (a < 2.06) {
			csvResult0point1();
		}
		else if (a < 2.485) {
			csvResult0point05();
		}
		else if (a < 2.787) {
			csvResult0point02();
		}
		else if (a < 3.45) {
			csvResult0point01();
		}
		else if (a < 3.725) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 26:
		if (a < 1.315) {
			csvResult0point2no();
		}
		else if (a < 1.706) {
			csvResult0point2yes();
		}
		else if (a < 2.056) {
			csvResult0point1();
		}
		else if (a < 2.479) {
			csvResult0point05();
		}
		else if (a < 2.779) {
			csvResult0point02();
		}
		else if (a < 3.435) {
			csvResult0point01();
		}
		else if (a < 3.707) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 27:
		if (a < 1.314) {
			csvResult0point2no();
		}
		else if (a < 1.703) {
			csvResult0point2yes();
		}
		else if (a < 2.052) {
			csvResult0point1();
		}
		else if (a < 2.473) {
			csvResult0point05();
		}
		else if (a < 2.771) {
			csvResult0point02();
		}
		else if (a < 3.421) {
			csvResult0point01();
		}
		else if (a < 3.689) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 28:
		if (a < 1.313) {
			csvResult0point2no();
		}
		else if (a < 1.701) {
			csvResult0point2yes();
		}
		else if (a < 2.048) {
			csvResult0point1();
		}
		else if (a < 2.467) {
			csvResult0point05();
		}
		else if (a < 2.763) {
			csvResult0point02();
		}
		else if (a < 3.408) {
			csvResult0point01();
		}
		else if (a < 3.674) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 29:
		if (a < 1.311) {
			csvResult0point2no();
		}
		else if (a < 1.699) {
			csvResult0point2yes();
		}
		else if (a < 2.045) {
			csvResult0point1();
		}
		else if (a < 2.462) {
			csvResult0point05();
		}
		else if (a < 2.756) {
			csvResult0point02();
		}
		else if (a < 3.396) {
			csvResult0point01();
		}
		else if (a < 3.66) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	case 30:
		if (a < 1.31) {
			csvResult0point2no();
		}
		else if (a < 1.697) {
			csvResult0point2yes();
		}
		else if (a < 2.042) {
			csvResult0point1();
		}
		else if (a < 2.457) {
			csvResult0point05();
		}
		else if (a < 2.75) {
			csvResult0point02();
		}
		else if (a < 3.385) {
			csvResult0point01();
		}
		else if (a < 3.646) {
			csvResult0point002();
		}
		else {
			csvResult0point001();
		}
		break;
	default:
		if (30 < n || n <= 60) {
			if (a < 1.296) {
				csvResult0point2no();
			}
			else if (a < 1.671) {
				csvResult0point2yes();
			}
			else if (a < 2) {
				csvResult0point1();
			}
			else if (a < 2.39) {
				csvResult0point05();
			}
			else if (a < 2.66) {
				csvResult0point02();
			}
			else if (a < 3.232) {
				csvResult0point01();
			}
			else if (a < 3.46) {
				csvResult0point002();
			}
			else {
				csvResult0point001();
			}
			break;
		}
		else if (n < 60 || n <= 120) {
			if (a < 1.289) {
				csvResult0point2no();
			}
			else if (a < 1.658) {
				csvResult0point2yes();
			}
			else if (a < 1.98) {
				csvResult0point1();
			}
			else if (a < 2.358) {
				csvResult0point05();
			}
			else if (a < 2.617) {
				csvResult0point02();
			}
			else if (a < 3.16) {
				csvResult0point01();
			}
			else if (a < 3.373) {
				csvResult0point002();
			}
			else {
				csvResult0point001();
			}
			break;
		}
		else if (n > 120) {
			if (a < 1.282) {
				csvResult0point2no();
			}
			else if (a < 1.645) {
				csvResult0point2yes();
			}
			else if (a < 1.96) {
				csvResult0point1();
			}
			else if (a < 2.326) {
				csvResult0point05();
			}
			else if (a < 2.576) {
				csvResult0point02();
			}
			else if (a < 3.091) {
				csvResult0point01();
			}
			else if (a < 3.291) {
				csvResult0point002();
			}
			else {
				csvResult0point001();
			}
			break;
		}
		else {
			cout << "自由度有误" << endl;
			break;
		}
	}
}

void csvResult0point2no() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.2的alpha水平下，无法拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point2yes() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.2的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point1() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.1的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point05() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.05的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point02() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.02的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point01() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.01的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point002() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.002的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}

void csvResult0point001() {
	initgraph(1000, 600, 1);
	HWND hn = GetHWnd();
	SetWindowText(hn, L"可视化结果");
	outtextxy(10, 50, L"说明：可以拒绝均值之间没有差异的原假设的alpha水平越小，犯错的概率越小，意味着两组数据均值差异越大。如果在0.2的alpha水平下");
	outtextxy(10, 65, L"仍无法拒绝原假设，则两组数据均值之间的差异非常小。一般如果能在0.05的alpha水平下拒绝原假设则认为证据是有力的");
	outtextxy(280, 300, L"在0.001的alpha水平下，可以拒绝均值之间没有差异的原假设");
	system("pause");
	closegraph();
}