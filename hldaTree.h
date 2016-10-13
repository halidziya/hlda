#pragma once
#include <list>
#include "Topic.h"
#include "utils.h"
using namespace std;

class hldaTree
{
public:
	list<Topic> topics;
	int ntopics;
	int nextid;
	void addToPath(vector<Doc>::iterator d, list<Topic>::iterator topic,int heldout=0);
	void removeFromPath(vector<Doc>::iterator d,int heldout=0);


	void likelihoods(vector<Doc>::iterator d);
	void ncrp();
	void cum(vector<Doc>::iterator doc);
	list<Topic>::iterator sample();

	
	void newTopic(list<Topic>::iterator parent);
	void removeTopic(list<Topic>::iterator topic);

	hldaTree();
	~hldaTree();
};



// We have a list in each topic
// How to remove easy
// MH : Just discard the list in the node 
// Gibbs : Add to last place and remove from initial, all nodes sorted in same way