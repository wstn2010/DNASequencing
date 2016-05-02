/*
	Experimetal Graph Library

*/


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

	bool existOutwardEdge(size_t v, char base)
	{
		assert(0 <= v && v < countVertex());
	
		Node *node = outwards[v];
		while (node)
		{
			if (info[node->vertex].base == base)
			{
				return true;
			}
			node = node->next;
		}

		return false;
	}

	size_t findOutwardVertex(size_t v, char base, bool& found)
	{
		assert(0 <= v && v < countVertex());
	
		Node *node = outwards[v];
		while (node)
		{
			if (info[node->vertex].base == base)
			{
				found = true;
				return node->vertex;
			}
			node = node->next;
		}

		found = false;
		return 0;
	}
};

// int main()
// {
// 	Graph g;

// 	cerr << "#vertex = " << g.countVertex() << endl;

// 	size_t v1 = g.addVertex('A', 1, 2);
// 	size_t v2 = g.addVertex('G', 1, 3);
// 	g.addEdge(g.root, v1);
// 	g.addEdge(g.root, v2);
// 	g.addEdge(v2, v1);

// 	g.showOutwardEdge(g.root);
// 	g.showOutwardEdge(v1);
// 	g.showOutwardEdge(v2);

// 	g.showinwardEdge(v1);

// 	cout << "existOutwardEdge(g.root, 'A') = " << g.existOutwardEdge(g.root, 'A') << endl;
// 	cout << "existOutwardEdge(g.root, 'T') = " << g.existOutwardEdge(g.root, 'T') << endl;

// 	return 0;
// }

int main()
{
	Graph g;

	cerr << "building graph" << endl;

	string currentChromatid = "TGTGAAGCCAGGTA";
	size_t len = 3;

	for (size_t start = 0; start < currentChromatid.size() - len; ++start)
	{
		size_t curr = g.root;
		for (size_t i = 0; i < len; ++i)
		{
			size_t pos = start + i;
			char base = currentChromatid[pos];

			bool found;
			size_t v = g.findOutwardVertex(curr, base, found);
			if (found)
			{
				curr = v;				
			}
			else
			{
				size_t v = g.addVertex(base, start, pos);
				g.addEdge(curr, v);
				curr = v;
			}
		}
	}

	g.showOutwardEdge(0);
	g.showOutwardEdge(1);
	g.showOutwardEdge(2);
	g.showOutwardEdge(3);

	return 0;
}
