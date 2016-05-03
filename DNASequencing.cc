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

#include <assert.h>
#include <cstdint>

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
size_t bits[256];


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

	ss << r.readName << "," << r.chromatidSequenceId << "," << (r.startPos + 1) << "," << (r.endPos + 1) 
	<< "," << (r.strand ? '+' : '-') << "," << fixed << setprecision(2) << r.score;

	return ss.str();
}

void setBitwiseRead(uint64_t bitwiseRead[], string& c, size_t len)
{
	for (size_t i = 0; i < len; ++i) 
	{
		size_t pos = i * 32;
		uint64_t v = 0;
		for (size_t dt = 0; dt < 32; ++dt)
		{
			v += bits[c[pos + dt]];
			if (dt != 31)
			{
				v <<= 2;
			}
		}

		bitwiseRead[i] = v;
	}
}
			
void setExtractedPart(uint64_t whole[], size_t len, uint64_t part[], size_t startPos)
{
	size_t base = startPos / 32;
	size_t offset = (startPos % 32) * 2;

	for (size_t i = 0; i < len; ++i) 
	{
		uint64_t cur = whole[base + i];
		uint64_t nxt = whole[base + i + 1];

		part[i] = ((cur << offset) & (~((1 << offset) - 1))) | (nxt >> (32 - offset));
	}
}

uint countBit(uint64_t val)
{
	// XXX: できれば64bitでやりたい
	uint32_t h = (uint32_t)(val >> 32);
	h = (h & 0x55555555) + ((h >> 1) & 0x55555555);
	h = (h & 0x33333333) + ((h >> 2) & 0x33333333);
	h = (h & 0x0f0f0f0f) + ((h >> 4) & 0x0f0f0f0f);
	h = (h & 0x00ff00ff) + ((h >> 8) & 0x00ff00ff);
	h = (h & 0x0000ffff) + ((h >>16) & 0x0000ffff);

	uint32_t l = (uint32_t)(val & UINT64_C(0x00000000ffffffff));
	l = (l & 0x55555555) + ((l >> 1) & 0x55555555);
	l = (l & 0x33333333) + ((l >> 2) & 0x33333333);
	l = (l & 0x0f0f0f0f) + ((l >> 4) & 0x0f0f0f0f);
	l = (l & 0x00ff00ff) + ((l >> 8) & 0x00ff00ff);
	l = (l & 0x0000ffff) + ((l >>16) & 0x0000ffff);

	return h + l;
}

uint countDiff(uint64_t v1[], uint64_t v2[], size_t len)
{
	uint cnt = 0;
	for (size_t i = 0; i < len; ++i)
	{
		cnt += countBit(v1[i] ^ v2[i]);
	}

	return cnt;
}



#define LEN_READ 150


class DNASequencing 
{

	vector<string> results;

	string chromatids;

	int currentChromatidSequenceId;

	vector< vector<size_t> > fastRef;

	uint64_t *bitwiseDNA;


public:

	DNASequencing()
	: fastRef(65536), bitwiseDNA(NULL)
	{
		bits['T'] = 0;
		bits['A'] = 1;
		bits['C'] = 2;
		bits['G'] = 3;
		bits['N'] = 0; // dummy
	}

	int passReferenceGenome(int chromatidSequenceId, vector<string> chromatidSequence) 
	{
		currentChromatidSequenceId = chromatidSequenceId;
		cerr << "passReferenceGenome: id=" << chromatidSequenceId << endl;

		chromatids = "";
		for (int i = 0; i < chromatidSequence.size(); ++i)
		{
			trim(chromatidSequence[i]);
			chromatids += chromatidSequence[i];
		}

		// align with 32chars
		size_t n = chromatids.size() % 32;
		if (n != 0)
		{
			chromatids += string(32 - n, 'N');
		}

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
		string& c = chromatids;

		size_t lenDNA = c.size();
		cerr << "lenDNA:" << lenDNA << endl;

		cerr << "preprocessing1..." << endl;

		// DNA has aligned.
		size_t bitwiseLen = lenDNA / 32;
		cerr << "bitwiseLen:" << bitwiseLen << endl;

		if (bitwiseDNA != NULL)
		{
			delete bitwiseDNA;
		}
		bitwiseDNA = new uint64_t[bitwiseLen];

		for (size_t i = 0; i < bitwiseLen; ++i) 
		{
			size_t pos = i * 32;
			uint64_t v = 0;
			for (size_t dt = 0; dt < 32; ++dt)
			{
				v += bits[c[pos + dt]];
				if (dt != 31)
				{
					v <<= 2;
				}
			}

			bitwiseDNA[i] = v;
		}

		cerr << "preprocessing2..." << endl;

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
			string normalRead = reads[i];
			string reverseRead = createReverseRead(normalRead);

			Result normalResult(readNames[i], currentChromatidSequenceId, true);
			Result reverseResult(readNames[i], currentChromatidSequenceId, false);

			align(normalResult, normalRead);
			align(reverseResult, reverseRead);

			results.push_back(toResultStr(normalResult.score > reverseResult.score ? normalResult : reverseResult));

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
		uint64_t bitwisePartialDNA[5];

		// 0. bitwiseReadを生成
		// align to x32
		uint64_t bitwiseRead[5];
		string alignedRead = r + string(10, 'N');
		setBitwiseRead(bitwiseRead, alignedRead, 5);

		size_t key = calcKey(r, 0);
		vector<size_t>& startPositions = fastRef[key];

		for (vector<size_t>::iterator it = startPositions.begin(); it != startPositions.end(); ++it)
		{
			size_t startPos = *it;

			// 1. bitwiseDNAをロード
			setExtractedPart(bitwiseDNA, 5, bitwisePartialDNA, startPos);

			// 3. xorと1-countでrate計算：あとは同じ
			uint diff = countDiff(bitwisePartialDNA, bitwiseRead, 5);
			cerr << "D" << hex << bitwisePartialDNA[0] << endl;
			cerr << "R" << hex << bitwiseRead[0] << endl;			
			cerr << "diff: " << diff << endl;
			float rate = (len - diff) / (float)len;

			// float rate = calcMatchRate(c, startPos, len, r);
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

		if (reads.size() == 20)
			break;

	}

	std::cerr << " loaded " << reads.size() << "reads" << std::endl;
	return true;
}

template<typename T, typename U>
void assertEqualHex(T expected, U value)
{
	cerr << ">> " << hex << value << endl << dec;
	assert(value == expected);
}

template<typename T, typename U>
void assertEqual(T expected, U value)
{
	cerr << ">> " << value << endl;
	assert(value == expected);
}

int main() {

	DNASequencing dnaSequencing;

	// unit test
	string test_srq = "CTAG";
	assertEqual(createReverseRead(test_srq), string("CTAG"));

	string test_srq2 = "TGTC";
	assertEqual(createReverseRead(test_srq2), string("GACA"));

	uint64_t whole[6];
	for (size_t i = 0; i < 6; ++i)
	{
		whole[i] = UINT64_C(0x00FF00FF00FF00FF);
	}
	uint64_t part[5];
	setExtractedPart(whole, 5, part, 0);
	assertEqualHex(part[0], UINT64_C(0x00FF00FF00FF00FF));

	setExtractedPart(whole, 5, part, 1);
	assertEqualHex(part[0],UINT64_C(0x03FC03FC03FC03FC));

	setExtractedPart(whole, 5, part, 32);
	assertEqualHex(part[0], UINT64_C(0x00FF00FF00FF00FF));

	setExtractedPart(whole, 5, part, 4);
	assertEqualHex(part[0], UINT64_C(0xFF00FF00FF00FF00));

	assertEqual(32, countBit(UINT64_C(0xFF00FF00FF00FF00)));

	uint64_t v1[1], v2[1];
	v1[0] = UINT64_C(0xFF00FF00FF00FF00);
	v2[0] = UINT64_C(0xFF00FF00FF01FF00);
	assertEqual(1, countDiff(v1, v2, 1));

	v2[0] = UINT64_C(0xFF10FF10FF01F000);
	assertEqual(7, countDiff(v1, v2, 1));

	string read1 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
	setBitwiseRead(part, read1, 1);
	assertEqual(0, part[0]);

	string read2 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAA";
	setBitwiseRead(part, read2, 1);
	assertEqual(5, part[0]);

	cerr << "unit test passed." << endl;

	// start evaluating
	clock_t begin = clock();

	vector<string> chromatidSequence;
	if (!loadChromatidSequence(chromatidSequence))
		return 0;

	vector<string> readNames;
	vector<string> reads;

	if (!loadReads(readNames, reads))
		return 0;

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

