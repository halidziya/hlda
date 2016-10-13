#include "Topic.h"

Topic::Topic(int id)
{
	ndocs = 0;
	nwords = 0; // Total number of words
	this->id = id;		// Topic ID
	loglikelihood = 0; // Gamma ratio
	restriction = 1;
	level = 0;
	words = Vector(nterms);
	words.zero();
}


Topic::Topic(list<Topic>::iterator parent, int id) : Topic(id)
{
	this->parent = parent;
	level = parent->level + 1;
}


Topic::~Topic()
{

}
