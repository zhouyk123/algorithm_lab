#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <chrono>
#include <cstdlib>
using namespace std;

int matchScore = 1; // 分数匹配
int mismatchScore = -21; // 分数不匹配
int gapPenalty = -22; // 间隙惩罚

void smithWaterman(string seq1, string seq2,vector<vector<int>> &dp,vector<int> &maxindex) {

    int m = seq1.size();
    int n = seq2.size();
    for (int i = 1; i <= m; ++i) {
        int score = dp[i][1];
        for (int j = 1; j <= n; ++j) {
            int match = dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? matchScore : mismatchScore);
            int deleteGap = dp[i - 1][j] + gapPenalty;
            int insertGap = dp[i][j - 1] + gapPenalty;

            dp[i][j] = max({0, match, deleteGap, insertGap});
            if(dp[i][j]>score){
                score=dp[i][j];
                maxindex[i]=j;
            }
        }
    }

    return;
}

pair<int,int> backtrace(vector<vector<int>>& dp,string seq1,string seq2, int i, int j) {
    int end=max(i-215,0);
    while (i > end && j > 0 ) {
        if (dp[i][j] == dp[i - 1][j - 1] + (seq1[i - 1] == seq2[j - 1] ? matchScore : mismatchScore)||dp[i][j]==0) {
            i--;
            j--;
        } else if (dp[i][j] == dp[i-1][j] + gapPenalty) {
            i--;
        } else {
            j--;
        }
    }
    return make_pair(i,j);
}


int main(int argc, char* argv[]) {
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <mismatchScore> <gapPenalty>" << endl;
        return 1;
    }

    mismatchScore = std::atoi(argv[1]);
    gapPenalty = std::atoi(argv[2]);

    ifstream  fin;
    fin.open("ref.txt",ios::in);
    string ref;
    fin>>ref;
    fin.close();
    fin.open("query.txt",ios::in);
    string tempt;
    fin>>tempt;
    fin.close(); 
    auto start = std::chrono::high_resolution_clock::now();
    int m = tempt.size();
    int n = ref.size();

    vector<pair<pair<int,int>,pair<int,int>>> result;
    vector<vector<int>> A(m + 1, vector<int>(n + 1, 0));
    vector<int> Amax(m + 1, 0);

    smithWaterman(tempt, ref,A, Amax);
    int index=0;
    pair<int,int> begin = backtrace(A,tempt,ref,m,n);
    int i = m;
    int j = Amax[m];
    while(i>=0&&j!=0){
        pair<int,int> begin = backtrace(A,tempt,ref,i,j);
        result.push_back(make_pair(begin,make_pair(i,j)));
        i=begin.first;
        j=Amax[begin.first];
    }
    for(i=1;i<result.size();i++){
        if(abs(result[i-1].first.second-result[i].second.second)<100){
            result[i-1].first.first=result[i].first.first;
            result[i-1].first.second=result[i].first.second;
            result.erase(result.begin()+i);
            i--;
        }
    }
    ofstream fout;
    fout.open("result.txt",ios::out);
    // for(i=0;i<result.size();i++){
    //     cout<<"("<<result[i].first.first<<", "<<result[i].second.first<<", ";
    //     cout<<result[i].first.second+160600000<<", "<<result[i].second.second+160600000<<"), ";
    // }
    // cout << endl << "[";
    // cout << "[";
    // for(i=0;i<result.size();i++){
    //     if (i!=0) cout << ", ";
    //     cout<<result[i].first.first<<", "<<result[i].second.first<<", ";
    //     cout<<result[i].first.second+160600000<<", "<<result[i].second.second+160600000;
    // }
    // cout << "]" << endl;
    fout << "[";
    for(i=0;i<result.size();i++){
        if (i!=0) fout << ", ";
        fout<<result[i].first.first<<", "<<result[i].second.first<<", ";
        // fout<<result[i].first.second+160600000<<", "<<result[i].second.second+160600000;
        fout<<result[i].first.second<<", "<<result[i].second.second;
    }
    fout << "]";
    fout.close();

    return 0;
}