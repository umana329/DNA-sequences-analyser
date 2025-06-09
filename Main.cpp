#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<iomanip>
#include<algorithm>
#include<graphics.h>
#undef max
#include<locale>
#include<codecvt>
#include<cmath>
#include"calculate.h"
using namespace std;
void readfasta(const string&, vector<string>&, vector<string>&);
void theResult(const wstring, const wstring,int);
void readCsv(const string&, vector<vector<string>>&);

int main() {
	cout << "输入0读取fasta文件，输入1读取csv文件" << endl;
	double xuanze;


	while (1) {
		while (!(cin >> xuanze)) {//防止各种奇怪的输入

			cin.clear();
			cin.ignore(numeric_limits<streamsize>::max(), '\n');
			cout << "请输入0或1" << endl;

		}

		if (xuanze == 0) {//fasta
			vector<string> seq;
			vector<string> seqInf;
			string filename = "C:/file/", filename1;
			cout << "输入要读取的文件名(记得加上\".fa\")：";
			cin >> filename1;
			filename += filename1;

			readfasta(filename, seq, seqInf);

			for (size_t i = 0; i < seqInf.size(); i++) {

				cout << i + 1 << "  " << "Sequence Info:" << seqInf[i] << endl;
				cout << "Sequence:" << seq[i] << endl;
				cout << "Length:" << seq[i].size() << endl;
				cout << "Proportion of GC:" << setprecision(3) << GCpercentMin(seq[i]) * 100 << "%";
				if (GCpercentMin(seq[i]) != GCpercentMax(seq[i])) {

					cout << " ～ " << setprecision(3) << GCpercentMax(seq[i]) * 100 << "%";

				}

				cout << "\n" << endl;

			}
			while (1) {

				cout << "输入1以进行序列比对，输入0退出程序" << endl;

				double choice;

				while (!(cin >> choice)) {

					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
					cout << "请输入0或1" << endl << "输入1以进行序列比对，输入0退出程序" << endl;

				}

				if (choice == 0) {
					break;
				}

				else if (choice == 1) {
					double s1, s2;

					while (1) {
						cout << "输入要比对的2个序列的编号：";
						while (!(cin >> s1 >> s2)) {
							cin.clear();
							cin.ignore(numeric_limits<streamsize>::max(), '\n');
							cout << "输入有误，请重新输入" << endl << "输入要比对的2个序列的编号：";
						}

						if (s1 <= 0 || s2 <= 0 || s1 > seq.size() || s2 > seq.size() || (s1 - int(s1) != 0) || (s2 - int(s2) != 0)) {
							cout << "输入有误，请重新输入" << endl;
							continue;
						}

						break;

					}

					auto result = needlemanWunsch(seq[s1 - 1], seq[s2 - 1]);

					cout << "score " << result.first << endl << "比对结果详见窗口" << endl;

					const string& str1 = result.second.first;
					const string& str2 = result.second.second;
					wstring_convert<codecvt_utf8<wchar_t>>converter;
					wstring wide1 = converter.from_bytes(str1);
					wstring wide2 = converter.from_bytes(str2);

					initgraph(1600, 800, 1);
					HWND hn = GetHWnd();
					SetWindowText(hn, L"比对结果");
					settextstyle(16, 0, _T("Courier New"));//令输出的字符等宽
					theResult(wide1, wide2, 10);
					system("pause");
					closegraph();
				}

				else {
					cout << "请输入0或1" << endl;
				}

			}

			break;
		}


		else if (xuanze == 1) {//csv
			vector<vector<string>>csv;
			string filename = "C:/file/", filename1;
			cout << "输入要读取的文件名(记得加上\".csv\")：";
			cin >> filename1;
			filename += filename1;

			readCsv(filename, csv);

			vector<vector<double>>data;
			vector<vector<double>>dataN;

			for (int r = 1; r < csv.size(); r++) {

				vector<double>record;
				vector<double>recordN;
				for (int c = 2; c < csv[0].size(); c++) {
					record.push_back(turnIntoDouble(csv[r][c]));
					recordN.push_back(turnIntoDouble(csv[r][c]));
				}

				normalization(recordN);
				data.push_back(record);
				dataN.push_back(recordN);

			}

			cout << "原始数据：" << endl;
			for (const auto& row : csv) {

				for (int i = 0; i < row.size();i++) {
					const auto& column = row[i];

					if (i == 0) {
						cout << setw(20) << left << column;
					}

					else {
						cout << setw(10) << left << column;
					}

				}

				cout << endl;
			}

			cout << "归一化后的数据：" << endl;
			for (const auto& row : dataN) {

				for (int i = 0; i < row.size(); i++) {
					const auto& column = row[i];
					cout << setw(10) << left << column<<" ";

				}
				cout << endl;
			}

			while (1) {

				cout << "输入1以对某个基因进行配对样本T检验，输入0退出程序" << endl;

				double choice;

				while (!(cin >> choice)) {

					cin.clear();
					cin.ignore(numeric_limits<streamsize>::max(), '\n');
					cout << "请输入0或1" << endl << "输入1以对某个基因进行配对样本T检验，输入0退出程序" << endl;

				}

				if (choice == 0) {
					break;
				}

				else if (choice == 1) {
					double s0,s1,s2,s3,s4;

					while (1) {
						cout << "输入要检验的基因所对应行数以及2组数据的列数范围：(例如要检测第2行的基因，第一组是第3～5列，第二组是6～8列则输入\"2 3 5 6 8\",注意：列数以原始数据表为准，2组数据的数量应该一致)"<<endl;
						while (!(cin >> s0 >> s1 >> s2 >> s3 >> s4)) {
							cin.clear();
							cin.ignore(numeric_limits<streamsize>::max(), '\n');
							cout << "输入有误，请重新输入" << endl ;
						}

						if (s0<=0||s0>csv.size() || (s0 - int(s0) != 0) || s1 <= 2 || s2 <= 2 || s1 > csv[0].size() || s2 > csv[0].size() || (s1 - int(s1) != 0) || (s2 - int(s2) != 0) || s3 <= 2 || s4 <= 2 || s3 > csv[0].size() || s4 > csv[0].size() || (s3 - int(s3) != 0) || (s4 - int(s4) != 0)) {
							cout << "输入有误，请重新输入" << endl;
							continue;
						}

						break;

					}

					vector<double>score1;
					vector<double>score2;
					const int n = s2 - s1;//样本量
					double sumOfD = 0,sumOfDsquare=0;//差值之和、平方差之和
					for (int i = s1; i <= s2; i++) {
						score1.push_back(data[s0][i - 3]);
					}

					for (int i = s3; i <= s4; i++) {
						score2.push_back(data[s0][i - 3]);
					}

					for (int i = 0; i < score1.size(); i++) {
						sumOfD += score1[i] - score2[i];
						sumOfDsquare += pow(score1[i] - score2[i], 2);
					}

					double t = (sumOfD / n) / (pow(((sumOfDsquare - (pow(sumOfD, 2) / n)) / (n * (n - 1))), 0.5));
					checkT(t, n - 1);
					break;
				}

				else {
					cout << "请输入0或1" << endl;
				}

			}

			break;
		}

		else {
			cout << "请输入0或1" << endl;
		}
	}
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

		return;
	}


	void theResult(const wstring w1, const wstring w2,int y) {

		if (w1.size() <= 198) {
			settextcolor(WHITE);
			outtextxy(10, y, w1.c_str());
			settextcolor(RED);
			outtextxy(10, y+20, w2.c_str());
			return;
		}

		else {
			wstring s11, s12, s21, s22;//198个字符一行
			s11 = w1.substr(0, 198);
			s21 = w2.substr(0, 198);
			s12 = w1.substr(198);
			s22 = w2.substr(198);
			settextcolor(WHITE);
			outtextxy(10, y, s11.c_str());
			settextcolor(RED);
			outtextxy(10, y + 20, s21.c_str());
			theResult(s12, s22, y + 50);
			return;
		}

	}


	void readCsv(const string& filename, vector<vector<string>>&records){
		ifstream file(filename);
		string line;

		while (getline(file, line)) {
			stringstream ss(line);
			vector<string>record;
			string cell;

			while (getline(ss, cell, ',')) {
				record.push_back(cell);
			}

			if (!record.empty()) {
				records.push_back(record);
			}

		}
		return;
	}
