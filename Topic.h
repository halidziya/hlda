#pragma once
#include <vector>
#include <list>
#include "Doc.h"

using namespace std;

class Topic
{
public:
	Vector words; // Frequency matrix of topic
	list<Topic>::iterator parent; // Parent node
	list<Topic>::iterator copy;   // To facilitate copy of tree
	list<vector<Doc>::iterator> docs;

	int ndocs;
	int nwords; // Total number of words
	int id;		// Topic ID
	int level;
	int restriction;

	// Calculation buffer
	double loglikelihood; // Gamma ratio
	double ncrp;
	double cumlog;
	double prob;


	//void addDoc(Doc& d); // Prototype
	//void removeDoc(Doc& d);

	Topic(int id);
	Topic(list<Topic>::iterator parent,int id);
	~Topic();
};

