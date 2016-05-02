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


using namespace std;

// 
// Lightweight graph lib
// 
// graph
struct Node
{
	Node *next;
	size_t vertex;
	Node(size_t to) : vertex(to), next(NULL) {}
};

struct Info
{
	char base; // ATCG
	size_t start; // start pos on reference DNA
	size_t end; // end pos on reference DNA
	Info(char base, size_t start, size_t end) : base(base), start(start), end(end) {}
};

// vertex0 = root

class Graph 
{
private:
	vector<Node *> outwards; 
	vector<Node *> inwards; 
	vector<Info> info;

public:
	size_t root;

	Graph() 
	{
		root = addVertex('*', 0, 0);
	}

	size_t addVertex(char base, size_t start, size_t end)
	{
		size_t idx = outwards.size();
		outwards.push_back(NULL);
		inwards.push_back(NULL);
		info.push_back(Info(base, start, end));

		return idx;
	}

	size_t countVertex()
	{
		return outwards.size();
	}

	void addOutwardsVertex(size_t from, size_t to)
	{
		Node *next = new Node(to);
		Node *node = outwards[from];
		if (node == NULL)
		{
			outwards[from] = next;
		}
		else
		{
			while (node->next != NULL)
			{
				node = node->next;
			}
			node->next = next;
		}
	}

	void addInwardsVertex(size_t from, size_t to)
	{
		Node *prev = new Node(from);
		Node *node = inwards[to];
		if (node == NULL)
		{
			inwards[to] = prev;
		}
		else
		{
			while (node->next != NULL)
			{
				node = node->next;
			}
			node->next = prev;
		}
	}

	void addEdge(size_t from, size_t to)
	{
		assert(0 <= from && from < countVertex());
		assert(0 <= to && to < countVertex());

		addOutwardsVertex(from, to);
		addInwardsVertex(from, to);
	}

	void showOutwardEdge(size_t v)
	{
		assert(0 <= v && v < countVertex());
	
		cout << "outward edge of vertex:" << v << "(" << info[v].base << ") = ";

		Node *node = outwards[v];
		while (node)
		{
			cout << node->vertex << "(" << info[node->vertex].base << ") ";
			node = node->next;
		}

		cout << endl;
	}

	void showinwardEdge(size_t v)
	{
		assert(0 <= v && v < countVertex());
	
		cout << "inward edge of vertex:" << v << " = ";

		Node *node = inwards[v];
		while (node)
		{
			cout << node->vertex << " ";
			node = node->next;
		}

		cout << endl;
	}

};

//
// end 
//


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



class DNASequencing 
{

  vector<string> results;
  
  int currentChromatidSequenceId;
  string currentChromatid;

  Graph g;

public:
      
	int passReferenceGenome(int chromatidSequenceId, vector<string> chromatidSequence) 
	{
		currentChromatid = "";

		for (int i = 0; i < chromatidSequence.size(); ++i)
		{
			trim(chromatidSequence[i]);
			currentChromatid += chromatidSequence[i];
		}

		currentChromatidSequenceId = chromatidSequenceId;

		cerr << "passReferenceGenome: id=" << chromatidSequenceId << endl;

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
		for (int i = 0; i < n; ++i)
		{
			cerr << "read" << i << endl;

			string normalRead = reads[i];
			string reverseRead = createReverseRead(normalRead);

			Result normalResult(readNames[i], currentChromatidSequenceId, true);
			Result reverseResult(readNames[i], currentChromatidSequenceId, false);

			align(normalResult, normalRead);
			if (normalResult.score >= 0.99)
			{
				results.push_back(toResultStr(normalResult));
			}
			else
			{
				align(reverseResult, reverseRead);

				if (normalResult.score > reverseResult.score)
				{
					results.push_back(toResultStr(normalResult));
				}
				else
				{
					results.push_back(toResultStr(reverseResult));
				}
			}

		}

		return results;
	}

private:

	float calcMatchRate(string& chroma, size_t pos, size_t len, string& read)
	{
		size_t cnt = 0;

		size_t lim = 15;

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
		string& c = currentChromatid;
		size_t len = r.size();
		size_t bestPos = -1;
		float bestRate = 0.0;

		for (size_t pos = 0; pos < c.size() - len;)
		{
			if (c[pos] == r[0] && c[pos + 1] == r[1] && c[pos + 2] == r[2])
			{
				float rate = calcMatchRate(c, pos, len, r);
				if (rate > bestRate)
				{
					bestRate = rate;
					bestPos = pos;
					cerr << "best pos=" << pos << " score=" << rate << endl;
				}
				if (rate == 1.0) {
					break;
				}
			}
			
			pos += 1;

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

	    if (reads.size() >= 10)
	    	break;
    }

    std::cerr << " loaded " << reads.size() << "reads" << std::endl;
	return true;
}

int main() {

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

	clock_t begin = clock();

	vector<string> results = dnaSequencing.getAlignment(reads.size(), 1.0, 1.0, readNames, reads);

	clock_t end = clock();


	cerr << "*** results *** " << (end - begin) / 1000 << "(K clocks)" << endl;
	cout << "*** results *** " << (end - begin) / 1000 << "(K clocks)" << endl;
	for (vector<string>::iterator it = results.begin(); it != results.end(); ++it)
	{
		cout << *it << endl;
	}

	return 1;
}

