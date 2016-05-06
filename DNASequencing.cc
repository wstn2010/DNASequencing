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

#define LOCAL_ONLY
//#define UNIT_TEST

#define LEN_READ 150
#define LEN_TARGET 32


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
	Result(string readName, bool strand)
	: readName(readName), startPos(0), endPos(0), strand(strand), score(-1)
	{
	}
};

#define LEN_READ 150
#define SZ_FAST_REF 16777216L

// utils
size_t bits[256];

void initBits()
{
  bits['T'] = 0; // 00
  bits['A'] = 1; // 01
  bits['C'] = 2; // 10
  bits['G'] = 3; // 11
  bits['N'] = 0; // dummy
}

inline size_t calcKey(string& c, size_t pos)
{
  size_t k1 = bits[c[pos + 0]] << 22;
  size_t k2 = bits[c[pos + 1]] << 20;
  size_t k3 = bits[c[pos + 2]] << 18;
  size_t k4 = bits[c[pos + 3]] << 16;
  size_t k5 = bits[c[pos + 4]] << 14;
  size_t k6 = bits[c[pos + 5]] << 12;
  size_t k7 = bits[c[pos + 6]] << 10;
  size_t k8 = bits[c[pos + 7]] <<  8;
  size_t k9 = bits[c[pos + 8]] <<  6;
  size_t ka = bits[c[pos + 9]] <<  4;
  size_t kb = bits[c[pos +10]] <<  2;
  size_t kc = bits[c[pos +11]] <<  0;

  return k1 + k2 + k3 + k4 + k5 + k6 + k7 + k8 + k9 + ka + kb + kc;
}

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
			
#define MASK_END UINT64_C(0xFFFFFFFFFFF00000)

// posは文字位置(!=ビット位置)
void setExtractedPart(uint64_t whole[], size_t len, uint64_t part[], size_t startPos)
{
	size_t base = startPos >> 5; // / 32;
	size_t offset = (startPos & 31) << 1; // (startPos % 32) * 2;

	for (size_t i = 0; i < len; ++i) 
	{
		uint64_t cur = whole[base + i];
		uint64_t nxt = whole[base + i + 1];

		part[i] = (cur << offset) | (nxt >> (64 - offset));
	}

	// お尻10文字は0fill: 20bit=2byte+4bit=0x00000
	part[len - 1] &= MASK_END;
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

#define MASK_L UINT64_C(0xAAAAAAAAAAAAAAAA)
#define MASK_H UINT64_C(0x5555555555555555)


// 2bitづつの区切りで、異なるときはそれぞれdiff=1としなくてはだめ
uint countDiff(uint64_t v1[], uint64_t v2[], size_t len)
{
	uint cnt = 0;
	for (size_t i = 0; i < len; ++i)
	{
		uint64_t val = v1[i] ^ v2[i];
		uint64_t adjusted = (val & MASK_H) | ((val & MASK_L) >> 1);

		cnt += countBit(adjusted);
	}

	return cnt;
}

// most left = 0
uint nlz(uint64_t x)
{
	x |= x >> 1;
	x |= x >> 2;
	x |= x >> 4;
	x |= x >> 8;
	x |= x >> 16;
	x |= x >> 32;

	return 64 - countBit(x);
}

// 最初の差異ビットの位置を返す：だいたいでよい
uint findStartDifferencePos(uint64_t v1[], uint64_t v2[], size_t len)
{
	uint cnt = 0;
	for (size_t i = 0; i < len; ++i)
	{
		uint64_t val = v1[i] ^ v2[i];
		uint64_t adjusted = (val & MASK_H) | ((val & MASK_L) >> 1);

		uint first = nlz(adjusted);
		if (first != 0)
		{
			return i * 64 + first;
		}
	}

	return -1;
}

// posからnビット分取り出す
// posは、v[]の通し位置: v[0]=0-63, v[1]=64-127,...
// 注：ベクタサイズ未指定
uint64_t extractRange(uint64_t whole[], size_t pos, size_t n)
{
	size_t base = pos >> 6; // / 64;
	size_t offset = pos & 63; //% 64;

	uint64_t cur = whole[base];
	uint64_t nxt = whole[base + 1];

	uint64_t result = (cur << offset) | (nxt >> (64 - offset));

	// mask
	size_t maskedBit = 64 - n;
	uint64_t mask = ~((UINT64_C(1) << maskedBit) - 1); 
	result &= mask;

	return result;
}

// posは、v[]の通し位置: v[0]=0-63, v[1]=64-127,...
// 注：ベクタサイズ決め打ち
void copy(uint64_t part[], uint64_t whole[], size_t pos, size_t n)
{
	size_t base = pos >> 6; // / 64;
	size_t offset = pos & 63; //% 64;

	size_t copiedBits = 0;
	size_t i = 0;
	for (; i < 5; ++i) 
	{
		uint64_t cur = whole[base + i];
		uint64_t nxt = whole[base + i + 1];

		part[i] = (cur << offset) | (nxt >> (64 - offset));
		copiedBits += 64;

		if (copiedBits > n)
		{
			// mask
			size_t maskedBit = copiedBits - n;
			uint64_t mask = ~((UINT64_C(1) << maskedBit) - 1);
			part[i++] &= mask;
			break;
		}
	}

	for (; i < 5; ++i)
	{
		part[i] = 0;
	}

}

// deletionを発見したら,その部分以外の差異数を返す
uint evaluateDeletion(uint64_t bitwisePartialDNA[], uint64_t bitwiseRead[], size_t n, uint pos)
{
	// cerr << "@dna : " << bitset<64>(bitwisePartialDNA[0]) << endl;
	// cerr << "@read: " << bitset<64>(bitwiseRead[0]) << endl;

	// 比較パターンを取り出す: READのposから、8文字(16bit): 後ろは0fill
	uint64_t target = extractRange(bitwiseRead, pos, LEN_TARGET);
	// cerr << "target:" << bitset<64>(target) << endl;

	// DNAのstartDifferencePos以降で一致個所を探す: 2bitずつ移動

	uint matchPos = 0;
	for (uint i = pos + 4; i < 200; i += 2)
	{
		uint64_t dna = extractRange(bitwisePartialDNA, i, LEN_TARGET);
		// cerr << "dna   :" << bitset<64>(dna) << endl;
		if ((dna ^ target) == 0)
		{
			matchPos = i;
			// cerr << "matchPos:" << dec << matchPos << endl;
			break;
		}

	}

	// 一致がないとか、ずれが4以下なら、でかい値を返す + 実行比較長が100以下
	if (matchPos == 0)
	{
		return 10000;
	}

	uint shift = matchPos - pos;
	// cerr << "shift:" << dec << shift << endl;

	// コピーし、ズラす。pos　＋ずらし分より前はマスク
	uint64_t copiedDNA[5], copiedRead[5];
	// copiedDNA: matchPosが先頭になるようにbitwiseReadから取り出し
	// copiedRead : posが先頭になるようにbitwiseReadから取り出し
	// 末尾matchpos分は、0マスク

	size_t effectiveReadBits = 300 - matchPos;
	// if (n == 1) // debug
	// 	effectiveReadBits = 64 - matchPos;

	copy(copiedDNA, bitwisePartialDNA, matchPos, effectiveReadBits);
	copy(copiedRead, bitwiseRead, pos, effectiveReadBits);
	// cerr << "+dna : " << bitset<64>(copiedDNA[0]) << endl;
	// cerr << "+read: " << bitset<64>(copiedRead[0]) << endl;

	// diffを求めて返す
	uint diff = countDiff(copiedDNA, copiedRead, n);

	return diff;		
}



uint evaluateInsertion(uint64_t bitwisePartialDNA[], uint64_t bitwiseRead[], size_t n, uint startDifferencePos)
{
	return 10000;		
}


//////////////////////////////////////////////////////////////////
// 
// class DNASequencing
//
//////////////////////////////////////////////////////////////////


class DNASequencing 
{

	vector<string> results;

	vector< vector<size_t> > fastRef;

	map< size_t, uint64_t *> dnaMap;

public:

	DNASequencing()
	: fastRef(SZ_FAST_REF)
	{
	  initBits();
	}

	int passReferenceGenome(int id, vector<string> chromatidSequence) 
	{
		cerr << "passReferenceGenome: id=" << id << endl;

		string chromatids = "";

		for (int i = 0; i < chromatidSequence.size(); ++i)
		{
			trim(chromatidSequence[i]);
			chromatids += chromatidSequence[i];
		}

		// align with 32chars
		size_t n = chromatids.size() & 31; //chromatids.size() % 32;
		if (n != 0)
		{
			chromatids += string(32 - n, 'N');
		}

		cerr << "creating bitwise DNA: id=" << id << endl;
		dnaMap[id] = createBitwiseDNA(chromatids);

		cerr << "creating fastRef: id=" << id << endl;
		setupFastRef(id, chromatids);

		return 0;
	}

	int initTest(int) 
	{
		return 0;
	}


	int preProcessing() 
	{
		return 0;
	}

	vector<string> getAlignment(int n, double normA, double normS, vector<string> readNames, vector<string> reads) 
	{
		cerr << "getAlignment: n=" << n << " normA=" << normA << " normS=" << normS << endl;

		results.clear();
		size_t cnt = 0;

		for (int i = 0; i < n; ++i)
		{
			string normalRead = reads[i];
			string reverseRead = createReverseRead(normalRead);

			Result normalResult(readNames[i], true);
			Result reverseResult(readNames[i], false);

			align(normalResult, normalRead);
			align(reverseResult, reverseRead);

			results.push_back(toResultStr(normalResult.score > reverseResult.score ? normalResult : reverseResult));

			if (++cnt % 10000 == 0)
			{
				cerr << ".";
			}
		}

		cerr << " done." << endl;

		return results;
	}

private:


	uint64_t *createBitwiseDNA(string& chromatids)
	{
		string& c = chromatids;

		size_t lenDNA = c.size();
		cerr << "lenDNA:" << lenDNA << endl;

		// DNA has aligned.
		size_t bitwiseLen = lenDNA / 32;
		cerr << "bitwiseLen:" << bitwiseLen << endl;

		uint64_t *bitwiseDNA = new uint64_t[bitwiseLen];

		for (size_t i = 0; i < bitwiseLen; ++i) 
		{
			size_t pos = i * 32;
			uint64_t v = 0;
			for (size_t dt = 0; dt < 32; ++dt)
			{
				v |= bits[c[pos + dt]];
				if (dt != 31)
				{
					v <<= 2;
				}
			}

			bitwiseDNA[i] = v;
		}

		return bitwiseDNA;
	}

	void setupFastRef(size_t id, string& chromatids)
	{
		size_t lenDNA = chromatids.size();

		for (size_t pos = 0; pos < lenDNA - LEN_READ; ++pos)
		{
			if (chromatids[pos] == 'N')
			{
				continue;
			}

			size_t key = calcKey(chromatids, pos);

			size_t val = (id << 32) | pos;
			fastRef[key].push_back(val);
		}
	}

	void align(Result& result, string& r)
	{
		uint64_t bitwisePartialDNA[5];

		// 0. bitwiseReadを生成
		// align to x32
		uint64_t bitwiseRead[5];
		string alignedRead = r + string(10, 'N');
		setBitwiseRead(bitwiseRead, alignedRead, 5);

		size_t key = calcKey(r, 0);
		vector<size_t>& startPositions = fastRef[key];

		if (startPositions.size() == 0)
		{
			// algorythm broken
			result.startPos = 1;
			result.endPos = 1 + LEN_READ;
			result.score = 0.0;
			result.chromatidSequenceId = 20;
			return;
		}

		int completedMatchCount = 0;
		int semiCompleteMatchCount = 0;

		size_t bestPos = -1;
		size_t bestId = 0;
		size_t bestDiff = 10000;

		for (vector<size_t>::iterator it = startPositions.begin(); it != startPositions.end(); ++it)
		{
			size_t startPos = *it;

			size_t id = startPos >> 32;

			startPos = startPos & 0x00000000ffffffff;

			uint64_t *bitwiseDNA = dnaMap[id];

			// 1. bitwiseDNAをロード
			setExtractedPart(bitwiseDNA, 5, bitwisePartialDNA, startPos);

			// 3. xorと1-countでrate計算：あとは同じ
			uint diff = countDiff(bitwisePartialDNA, bitwiseRead, 5);
			switch (diff)
			{
				case 0: ++completedMatchCount; break;
				case 1: ++semiCompleteMatchCount; break;
			}


			// search deletion
			// if (diff > 4) 
			// {
			// 	// 通しビット：24 < v < 300 - 80 までの値
			// 	uint startDifferencePos = findStartDifferencePos(bitwisePartialDNA, bitwiseRead, 5);
			// 	if (startDifferencePos != -1 && (24 + 4) < startDifferencePos && startDifferencePos < (300 - 80))
			// 	{
			// 		// 欠落一致トライ後のdiffを求める: pos - 1の1は、差異は2bitTACGの下位ビットに立つため
			// 		uint newDiff = evaluateDeletion(bitwisePartialDNA, bitwiseRead, 5, startDifferencePos - 1);
			// 		if (newDiff < diff)
			// 		{
			// 			diff = newDiff;

			// 			// if (diff <= 1)
			// 			// {
			// 			// 	result.startPos = startPos;
			// 			// 	result.endPos = startPos + LEN_READ;
			// 			// 	result.score = (1.0 - (diff / (float)LEN_READ)) * 0.9;
			// 			// 	result.chromatidSequenceId = id;
			// 			// 	return;
			// 			// }

			// 		} 
			// 		else
			// 		{
			// 			// 挿入一致トライ後のdiffを求める: pos - 1の1は、差異は2bitTACGの下位ビットに立つため
			// 			newDiff = evaluateInsertion(bitwisePartialDNA, bitwiseRead, 5, startDifferencePos - 1);
			// 			if (newDiff < diff)
			// 			{
			// 				// found an insertion
			// 				diff = newDiff;
			// 			}
			// 		}
			// 	}
			// }

			if (diff < bestDiff)
			{
				bestDiff = diff;				
				bestPos = startPos;
				bestId = id;
			}
		}

		if (completedMatchCount > 1)
		{
			result.startPos = bestPos;
			result.endPos = bestPos + LEN_READ;
			result.score = 1.0 / completedMatchCount;
			result.chromatidSequenceId = bestId;
		} 
		else if (semiCompleteMatchCount > 1)
		{
			result.startPos = bestPos;
			result.endPos = bestPos + LEN_READ;
			result.score = 1.0 / semiCompleteMatchCount;
			result.chromatidSequenceId = bestId;
		} 
		else
		{
			result.startPos = bestPos;
			result.endPos = bestPos + LEN_READ;
			result.score = (1.0 - (bestDiff / (float)LEN_READ)) / (bestDiff + 1);
			result.chromatidSequenceId = bestId;
		}

	}

};


#ifdef LOCAL_ONLY
/**
 * Constants from the problem statement
 */
const int MAX_POSITION_DIST = 300;
const double NORM_A_SMALL = -3.392;
const double NORM_A_MEDIUM = -3.962;
const double NORM_A_LARGE = -2.710;
const double MAX_AUC = 0.999999;

/**
 * Position: describe the position of a read within the genome
 */
struct Position {
	int rname;
	int from;
	int to;
	char strand;
};

/**
 * ReadResult: result of a read alignment
 */
struct ReadResult {
	double confidence;
	int r;
};

/**
 * Split a comma-separated string into a vector of string
 * @param row	the string to be split
 * @return	the vector of string
 */
vector<string> tokenize(const string& row) {
	vector<string> tokens;
	for(int i=0, pos=0, n=row.size(); i<n; ++i) {
		if(i==n-1 || row[i+1]==',') {
			string token = row.substr(pos, (i+1)-pos);
			tokens.push_back(token);
			pos = i+2;
		}
	}
	return tokens;
}

/**
 * Read a minisam file and build a map of ground truth
 * @param path	the path of the minisam file storing the ground truth 
 * @return a map[read_name] = read_Position
 */
map<string, Position> parse_truth(const string& path) {
	map<string, Position> res;
	ifstream ifs(path);
	string s;
	while(ifs >> s) {
		vector<string> tokens = tokenize(s);
		try {
			string qname = tokens[0];
			int chromatid = stoi(tokens[1]);
			int from = stoi(tokens[2]);
			int to = stoi(tokens[3]);
			char strand = tokens[4][0];
			res[qname] = Position{chromatid, from, to, strand};
		} catch(exception& e) {
			;
		}
	}
	return res;
}

/**
 * For each string of the results vector, build a read result {confidence, r}
 * @param truth		the map of ground truth position for each read
 * @param results	the vector of results as return by getAlignment
 * @return a vector of ReadResult, that is {confidence, r}
 */
vector<ReadResult> build_read_results(const map<string, Position>& truth, const vector<string>& results) {
	vector<ReadResult> read_results;
	int n = results.size();
	int correct = 0;
	for(int i=0; i<n; ++i) {
		vector<string> tokens = tokenize(results[i]);
		auto p = truth.find(tokens[0]);
		const Position& position = p->second;
		int r = 1;
		r = (stoi(tokens[1])==position.rname) ? r : 0;
		r = (tokens[4][0]==position.strand) ? r : 0;
		int start0 = stoi(tokens[2]);
		int start1 = position.from;
		r = (abs(start0-start1)<MAX_POSITION_DIST) ? r : 0;
		double confidence = stod(tokens[5]);
		read_results.push_back(ReadResult{confidence, r});
		correct += r;
	}
	cerr << "Number of correct answers: " << correct << '/' << n << " = " << (double)correct/(double)n << endl;
	return read_results;
}

/**
 * Compute the accuracy given the {confidence, r} pairs and the normalization facto
 * @param read_results	a vector of {confidence, r} results
 * @param norm_a		as described in the problem statement
 * @return	a double, the computed accuracy
 */
double compute_accuracy(vector<ReadResult>& read_results, double norm_a) {
	int n = read_results.size();
	sort(read_results.begin(), read_results.end(), 
		[](const ReadResult& lhs, const ReadResult& rhs){ return (lhs.confidence>rhs.confidence);});
	// merge results of equal confidence
	vector<int> cumul_si{read_results[0].r};
	vector<int> pos{0};
	for(int i=1; i<n; ++i) {
		if(read_results[i].confidence==read_results[i-1].confidence) {
			cumul_si.back() += read_results[i].r;
			pos.back() = i;
		} else {
			double cumul = cumul_si.back() + read_results[i].r;
			cumul_si.push_back(cumul);
			pos.push_back(i);
		}
	}
	// compute the AuC
	double auc = 0.0;
	double invn = 1.0 / (double)n;
	double invnp1 = 1.0 / (double)(n+1);
	double lfmultiplier = 1.0 / log(n+1);
	int m = cumul_si.size();
	for(int i=0; i<m; ++i) {
		double fi = 1.0 * (2+pos[i] - cumul_si[i])  * invnp1;
		double fi1 = (i==m-1) ? 1.0 : 1.0 * (2+pos[i+1] - cumul_si[i+1]) * invnp1;
		double lfi = lfmultiplier * log(fi);
		double lfi1 = lfmultiplier * log(fi1);
		auc += cumul_si[i] * (lfi1 - lfi) * invn;
	}
	cout << "auc = " << auc << endl;
	double tmp = log(1 - min(auc, MAX_AUC));
	cout << "log(1 - min(auc, MAX_AUC)) = " << tmp << endl; 
	cout << "NormA = " << norm_a << endl;
	double accuracy = tmp / norm_a;
	cout << "accuracy = " << accuracy << endl;
	return accuracy;
}

/**
 * Perform a single test
 * @param testDifficulty	define the test type (SMALL=0, MEDIUM=1, LARGE=2)
 * @return	alignments in format specified in the problem statement
 */
double time_cutoff;

vector<string> perform_test(int testDifficulty, double norm_a) {
	// test data path and description
	string fa1_path, fa2_path;
	vector<int> chr_ids;	
	if(testDifficulty==0) {
		fa1_path = "./data/small5.fa1";
		fa2_path = "./data/small5.fa2";
		chr_ids = vector<int>{20};
	} else if(testDifficulty==1) {
		fa1_path = "./data/medium5.fa1";
		fa2_path = "./data/medium5.fa2";
		chr_ids = vector<int>{1,11,20};
	} else if(testDifficulty==2) {
		fa1_path = "./data/large5.fa1";
		fa2_path = "./data/large5.fa2";
		for(int i=1; i<=24; ++i) chr_ids.push_back(i);		
	}	
	// call the MM DNASequencing methods
	DNASequencing dna_sequencing;
	dna_sequencing.initTest(testDifficulty);
	// load chromatid	
	for(int chromatid_seq_id: chr_ids) {
		vector<string> chromatid_seq;
		string path = "./data/chromatid" + to_string(chromatid_seq_id) + ".fa";
		ifstream ifs(path);
		string s;
		// skip header
		getline(ifs, s);
		cerr << "Skip header: " << s << endl;
		// pack all lines in chromatid_seq
		for(int i=0;getline(ifs, s); ++i) {
			if(s.back()=='\r') s.pop_back();
			chromatid_seq.push_back(s);
		}
		dna_sequencing.passReferenceGenome(chromatid_seq_id, chromatid_seq);		
	}
	dna_sequencing.preProcessing();
	// load reads
	vector<string> read_id, read_seq;
	{
		ifstream ifs1(fa1_path);
		ifstream ifs2(fa2_path);
		string s1, s2;
		while(getline(ifs1, s1) && getline(ifs2, s2)) {
			if(s1.back()=='\r') s1.pop_back();
			if(s2.back()=='\r') s2.pop_back();
			read_id.push_back(s1.substr(1, s1.size()-1));
			read_id.push_back(s2.substr(1, s2.size()-1));
			getline(ifs1, s1);
			getline(ifs2, s2);
			if(s1.back()=='\r') s1.pop_back();
			if(s2.back()=='\r') s2.pop_back();
			read_seq.push_back(s1);		
			read_seq.push_back(s2);
		}
	}
	int nreads = read_id.size();
	// compute alignments
	clock_t start_clock = clock();
	vector<string> results = dna_sequencing.getAlignment(nreads, norm_a, 0.5, read_id, read_seq);
	clock_t end_clock = clock();
	time_cutoff = (end_clock - start_clock) / (double) CLOCKS_PER_SEC;
	cerr << "time cutoff:" << time_cutoff << "(sec)" << endl; 

	return results;
}

/**
 * Main function: read the data, perform the DNA alignments and score results
 */
int main() {
	const int testDifficulty = 0;
	string minisam_path;
	double norm_a;
	double norm_s = 1.0; // 0.5
	double test_norm; // TestNorm = 1 000/1.05, 1 000 000/1.05, and 1 000 000/1.05
	double time_cut_off; // 16.1, 1102, 13730
	if(testDifficulty==0) {
		minisam_path = "./data/small5.minisam";
		norm_a = NORM_A_SMALL;
		test_norm = 1000/1.05;
		time_cut_off = 16.1;
	} else if(testDifficulty==1) {
		minisam_path = "./data/medium5.minisam";
		norm_a = NORM_A_MEDIUM;
		test_norm = 1000000/1.05;
		time_cut_off = 1102;
	} else if(testDifficulty==2) {
		minisam_path = "./data/large5.minisam";
		norm_a = NORM_A_LARGE;
		test_norm = 1000000/1.05;
		time_cut_off = 13730;
	}
	// perform test
	vector<string> results = perform_test(testDifficulty, norm_a);
	// load truth
	map<string, Position> truth = parse_truth(minisam_path);	
	vector<ReadResult> read_results = build_read_results(truth, results);
	// scoring
	double accuracy = compute_accuracy(read_results, norm_a);
	cerr << "accuracy = " << accuracy << endl; 
	double speed = 1.0 / norm_s * (1.0 - time_cutoff / time_cut_off);
	cerr << "speed = " << speed << endl;
	cerr << "score = " << (test_norm * accuracy * speed) << endl;

	cout << "*** results *** " << time_cutoff << "(sec)" << endl;
	for (vector<string>::iterator it = results.begin(); it != results.end(); ++it)
	{
		cout << *it << endl;
	}

	return 0;
}
#endif


/**********************************************************************************************

	 test bench

***********************************************************************************************/

#ifdef UNIT_TEST

#define assertEqualHex(T,U) { cerr << ">> " << hex << U << endl << dec; assert(T == U); }

#define assertEqual(T,U) { cerr << ">> " << U << endl; assert(T == U); } 

int main() 
{
	initBits();

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
	// setExtractedPart(whole, 5, part, 0);
	// assertEqualHex(part[0], UINT64_C(0x00FF00FF00FF00FF));

	setExtractedPart(whole, 5, part, 1);
	assertEqualHex(part[0],UINT64_C(0x03FC03FC03FC03FC));

	// setExtractedPart(whole, 5, part, 32);
	// assertEqualHex(part[0], UINT64_C(0x00FF00FF00FF00FF));

	setExtractedPart(whole, 5, part, 4);
	assertEqualHex(part[0], UINT64_C(0xFF00FF00FF00FF00));

	// test: countBit

	assertEqual(32, countBit(UINT64_C(0xFF00FF00FF00FF00)));

	// test: countDiff

	uint64_t v1[1], v2[1];
	v1[0] = UINT64_C(0xFF00FF00FF00FF00);
	v2[0] = UINT64_C(0xFF00FF00FF02FF00);
	assertEqual(1, countDiff(v1, v2, 1));

	v2[0] = UINT64_C(0xFF10FF10FF00F000);
	assertEqual(4, countDiff(v1, v2, 1));

	// test: setBitwiseRead

	string read1 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
	setBitwiseRead(part, read1, 1);
	assertEqual(0, part[0]);

	string read2 = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTAA";
	setBitwiseRead(part, read2, 1);
	assertEqual(5, part[0]);

	// test: nlz

	assertEqual(7, nlz(UINT64_C(0x0100000000000000)));

	assertEqual(1, nlz(UINT64_C(0x6000000000000000)));

	assertEqual(1, nlz(UINT64_C(0x60C0F0C0F0C0F0C0)));

	// test: findStartDifferencePos

	uint pos;
	v1[0] = UINT64_C(0x7F00FF00FF00FF00);
	v2[0] = UINT64_C(0xFF00FF00FF02FF00);
	pos = findStartDifferencePos(v1, v2, 1);
	assertEqual(1, pos);

	v1[0] = UINT64_C(0xFF00FF00FF00FF00);
	v2[0] = UINT64_C(0xFF00FF00FF02FF00);
	pos = findStartDifferencePos(v1, v2, 1);
	assertEqual(47, pos);

	// test: extractRange
	whole[0] = UINT64_C(0xFFFFFFFFFFFFFFFF);
	whole[1] = UINT64_C(0x0000000000000000);
	assertEqualHex(UINT64_C(0xFFFF000000000000), extractRange(whole, 32, 16));

	assertEqualHex(UINT64_C(0xFFFF000000000000), extractRange(whole, 48, 32));

	// test: copy
	for (size_t i = 0; i < 6; ++i)
	{
		whole[i] = UINT64_C(0x8811223344556677);
	}
	copy(part, whole, 48, 32);
	assertEqualHex(UINT64_C(0x6677881100000000), part[0]);
	assertEqualHex(UINT64_C(0x0000000000000000), part[1]);

	// test: evaluateDeletion
	v1[0] = UINT64_C(0x0011223344556677);
	v2[0] = UINT64_C(0x0011334455667788);
	pos = findStartDifferencePos(v1, v2, 1);
	assertEqual(19, pos);
	uint diff = evaluateDeletion(v1, v2, 1, pos - 1);
	assertEqual(0, diff);

	cerr << "unit test passed." << endl;

	return 1;
}

#endif
