#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <utility>
#include <set>
#include <cctype>
#include <queue>
#include <stack>
#include <cstdio>
#include <cstdlib>

#include <cmath>
#include <fstream>
#include <functional> 
#include <locale>
#include <ctime>
#include <iomanip>


using namespace std;

// trim from start (in place)
static inline void ltrim(string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
}

// trim from end (in place)
static inline void rtrim(string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
    ltrim(s);
    rtrim(s);
}

struct Result {
	string readName;
	int chromatidSequenceId;
	size_t startPos;
	size_t endPos;
	bool strand;
	float score;
	Result(string readName, int chromatidSequenceId, bool strand)
	: readName(readName), strand(strand), chromatidSequenceId(chromatidSequenceId), startPos(0), endPos(0), score(-1)
	{
	}
};

// utils

string createReverseRead(string read)
{
	string s;
	size_t len = read.size();
	s.resize(len);

	reverse_copy(read.begin(), read.end(), s.begin());

	for (size_t i = 0; i < len; ++i)
	{
		switch (s[i])
		{
			case 'T': s[i] = 'A'; break;
			case 'A': s[i] = 'T'; break;
			case 'C': s[i] = 'G'; break;
			case 'G': s[i] = 'C'; break;
		}
	}

	return s;
}

string toResultStr(Result& r)
{
	stringstream ss;

	ss << r.readName << "," << r.chromatidSequenceId << "," << (r.startPos + 1) << "," << (r.endPos + 1) << "," << (r.strand ? '+' : '-') << "," << fixed << setprecision(2) << r.score;

	return ss.str();
}

#define LEN_READ 150


class DNASequencing 
{

	vector<string> results;
  
	string chromatids;

	int currentChromatidSequenceId;

	vector< vector<size_t> > fastRef;

	size_t bits[256];

public:

	DNASequencing()
	: fastRef(65536)
	{
		bits['T'] = 0;
		bits['A'] = 1;
		bits['C'] = 2;
		bits['G'] = 3;
		bits['N'] = 0; // dummy
	}
      
	int passReferenceGenome(int chromatidSequenceId, vector<string> chromatidSequence) 
	{
		chromatids = "";
		for (int i = 0; i < chromatidSequence.size(); ++i)
		{
			trim(chromatidSequence[i]);
			chromatids += chromatidSequence[i];
		}

		currentChromatidSequenceId = chromatidSequenceId;

		cerr << "passReferenceGenome: id=" << chromatidSequenceId << endl;

		return 0;
	}

	int initTest(int) 
	{
		return 0;
	}

	size_t calcKey(string& c, size_t offset)
	{
		size_t pos = offset;
		size_t los = pos + LEN_READ - 1;

		size_t k1 = bits[c[pos + 0]] << 14;
		size_t k2 = bits[c[pos + 1]] << 12;
		size_t k3 = bits[c[pos + 2]] << 10;
		size_t k4 = bits[c[pos + 3]] <<  8;
		size_t k5 = bits[c[los - 3]] <<  6;
		size_t k6 = bits[c[los - 2]] <<  4;
		size_t k7 = bits[c[los - 1]] <<  2;
		size_t k8 = bits[c[los - 0]] <<  0;

		size_t key = k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8;
		return key;
	}

	int preProcessing() 
	{
		cerr << "preprocessing..." << endl;

		size_t lenDNA = chromatids.size();

		size_t max = 0;

		for (size_t pos = 0; pos < lenDNA - LEN_READ; ++pos)
		{
			if (chromatids[pos] == 'N')
			{
				continue;
			}

			size_t key = calcKey(chromatids, pos);

			fastRef[key].push_back(pos);

			size_t sz = fastRef[key].size();
			max = max > sz ? max : sz;
		}

		cerr << "done. max-len:" << max << endl;

		return 0;
	}

	vector<string> getAlignment(int n, double normA, double normS, vector<string> readNames, vector<string> reads) 
	{
		size_t cnt = 0;

		for (int i = 0; i < n; ++i)
		{
			// cerr << "read" << i << endl;

			string normalRead = reads[i];
			string reverseRead = createReverseRead(normalRead);

			Result normalResult(readNames[i], currentChromatidSequenceId, true);
			Result reverseResult(readNames[i], currentChromatidSequenceId, false);

			align(normalResult, normalRead);
			align(reverseResult, reverseRead);

			if (normalResult.score > reverseResult.score)
			{
				results.push_back(toResultStr(normalResult));
			}
			else
			{
				results.push_back(toResultStr(reverseResult));
			}

			if (++cnt % 1000 == 0)
			{
				cerr << ".";
			}
		}

		cerr << " done." << endl;

		return results;
	}

private:

	float calcMatchRate(string& chroma, size_t pos, size_t len, string& read)
	{
		size_t cnt = 0;

		size_t lim = 150/2;

		for (size_t i = 0; i < len; ++i)
		{
			if (chroma[pos + i] == read[i])
				++cnt;
			else  {
				if (--lim <= 0) {
					cnt = 0;
					break;
				}

			}
		}

		return (float)cnt / len;
	}

	// readは正順のみ
	void align(Result& result, string& r)
	{
		string& c = chromatids;
		size_t lenDNA = c.size();
		size_t len = r.size();
		size_t bestPos = -1;
		float bestRate = 0.0;

		size_t key = calcKey(r, 0);
		vector<size_t>& startPositions = fastRef[key];

		for (vector<size_t>::iterator it = startPositions.begin(); it != startPositions.end(); ++it)
		{
			size_t startPos = *it;

			float rate = calcMatchRate(c, startPos, len, r);
			if (rate > bestRate)
			{
				bestRate = rate;
				bestPos = startPos;
				// cerr << "best pos=" << startPos << " score=" << rate << endl;
			}
			if (rate == 1.0) {
				break;
			}
		}

		result.startPos = bestPos;
		result.endPos = bestPos + len;
		result.score = bestRate;
	}

};


/**********************************************************************************************

	 test bench

***********************************************************************************************/

bool loadChromatidSequence(vector<string>& v)
{
    const char *filename = "chromatid20.fa";

    std::ifstream ifs(filename);
    std::string str;
    if (ifs.fail())
    {
        std::cerr << "失敗" << std::endl;
        return false;
    }

    int lines = 0;

    getline(ifs, str); // skip header
    while (getline(ifs, str))
    {
        if (++lines % 10000 == 0)
            std::cerr << ".";

        v.push_back(str);
    }

    std::cerr << " loaded" << std::endl;
    return true;
}

bool loadReads(vector<string>& readNames, vector<string>& reads)
{
   const char *basename = "small5";
   string name1(basename);
   string name2(basename);

   name1.append(".fa1");
   name2.append(".fa2");

    std::ifstream ifs1(name1);
    std::ifstream ifs2(name2);
    if (ifs1.fail() || ifs2.fail())
    {
        std::cerr << "失敗" << std::endl;
        return false;
    }

    int lines = 0;
    std::string str1;
    std::string str2;

    while (getline(ifs1, str1) && getline(ifs2, str2))
    {
        if (++lines % 10000 == 0)
            std::cerr << ".";

        readNames.push_back(str1.substr(1));
        readNames.push_back(str2.substr(1));

	    getline(ifs1, str1) && getline(ifs2, str2);
	    reads.push_back(str1);
	    reads.push_back(str2);

    }

    std::cerr << " loaded " << reads.size() << "reads" << std::endl;
	return true;
}

int main() {

	clock_t begin = clock();

	vector<string> chromatidSequence;
	if (!loadChromatidSequence(chromatidSequence))
		return 0;

	vector<string> readNames;
	vector<string> reads;

	if (!loadReads(readNames, reads))
		return 0;

    DNASequencing dnaSequencing;

    dnaSequencing.passReferenceGenome(20, chromatidSequence);

    dnaSequencing.initTest(0);

    dnaSequencing.preProcessing();

	clock_t end = clock();
	cerr << "*** preprocess *** " << (end - begin) / (double) CLOCKS_PER_SEC << "(sec)" << endl;

	begin = clock();

	vector<string> results = dnaSequencing.getAlignment(reads.size(), 1.0, 1.0, readNames, reads);

	end = clock();

	cerr << "*** results *** " << (end - begin) / (double) CLOCKS_PER_SEC << "(sec)" << endl;
	cout << "*** results *** " << (end - begin) / (double) CLOCKS_PER_SEC << "(sec)" << endl;
	for (vector<string>::iterator it = results.begin(); it != results.end(); ++it)
	{
		cout << *it << endl;
	}

	return 1;
}

