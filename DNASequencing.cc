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
// #define TEST_BENCH

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
	: readName(readName), chromatidSequenceId(chromatidSequenceId), startPos(0), endPos(0), strand(strand), score(-1)
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
			
#define MASK_END UINT64_C(0xFFFFFFFFFFF00000)


void setExtractedPart(uint64_t whole[], size_t len, uint64_t part[], size_t startPos)
{
	size_t base = startPos / 32;
	size_t offset = (startPos % 32) * 2;

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
				v |= bits[c[pos + dt]];
				if (dt != 31)
				{
					v <<= 2;
				}
			}

			bitwiseDNA[i] = v;

			// verified: c.substr(pos, 32) = bitset<64>(v) 

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
			// if (i != 3) continue;

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

	// readは正順のみ
	void align(Result& result, string& r)
	{
		string& c = chromatids;
		size_t lenDNA = c.size();
		size_t len = r.size();
		size_t bestPos = -1;
		float bestRate = 0.0;
		uint64_t bitwisePartialDNA[5];
		uint64_t bestBitwisePartialDNA[5];


		// 0. bitwiseReadを生成
		// align to x32
		uint64_t bitwiseRead[5];
		string alignedRead = r + string(10, 'N');
		setBitwiseRead(bitwiseRead, alignedRead, 5);

		size_t key = calcKey(r, 0);
		vector<size_t>& startPositions = fastRef[key];

		int completedMatchCount = 0;
		int semiCompleteMatchCount = 0;

		for (vector<size_t>::iterator it = startPositions.begin(); it != startPositions.end(); ++it)
		{
			size_t startPos = *it;

			// 1. bitwiseDNAをロード
			setExtractedPart(bitwiseDNA, 5, bitwisePartialDNA, startPos);

			// 3. xorと1-countでrate計算：あとは同じ
			uint diff = countDiff(bitwisePartialDNA, bitwiseRead, 5);
			float rate = (len - diff) / (float)len;

			// float rate = calcMatchRate(c, startPos, len, r);
			if (rate > bestRate)
			{
				bestRate = rate;
				bestPos = startPos;
			}

			switch (diff)
			{
				case 0: ++completedMatchCount; break;
				case 1: ++semiCompleteMatchCount; break;
			}

		}

		// 完全一致が複数ある場合
		if (completedMatchCount > 1)
		{
			bestRate = 1.0 / completedMatchCount;
		} 
		else if (semiCompleteMatchCount > 1)
		{
			bestRate = 1.0 / semiCompleteMatchCount;
		} 

		result.startPos = bestPos;
		result.endPos = bestPos + len;
		result.score = bestRate;
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

#ifdef TEST_BENCH

bool loadChromatidSequence(vector<string>& v)
{
	const char *filename = "./data/chromatid20.fa";

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
	const char *basename = "./data/small5";
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

		// if (reads.size() == 20)
		// 	break;

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
	retry = false;
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

	assertEqual(32, countBit(UINT64_C(0xFF00FF00FF00FF00)));

	uint64_t v1[1], v2[1];
	v1[0] = UINT64_C(0xFF00FF00FF00FF00);
	v2[0] = UINT64_C(0xFF00FF00FF02FF00);
	assertEqual(1, countDiff(v1, v2, 1));

	v2[0] = UINT64_C(0xFF10FF10FF00F000);
	assertEqual(4, countDiff(v1, v2, 1));

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

#endif
