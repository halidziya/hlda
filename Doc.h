#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm> 
#include <list>
#include "params.h"
#include "utils.h"
#include "Matrix.h"


using namespace std;

class Topic;
class Doc
{
public:

	vector<int> words;		  // Actual words
	//vector<int> heldout;	  // Held out portion
	vector<int> z;			  // Topic assignments path+LDA
	//vector<int> heldz;			  // Topic assignments path+LDA
	double like;			  // likelihood of words
	vector<int> indexset;	  // Index of words , sparse index

	Vector levelsum;		  // how many words in one level
	Vector emptyscore;		  // For empty topics hold the scores here
	Matrix bow;				  // Bag of Words, levels , full matrix


	int ndiffwords;			  // Number of different words in level
	int length;				  // Document length
	int validpath=0;
	//int nheld;
	vector<list<Topic>::iterator> path;	  // Document path in the tree
	vector<Doc>::iterator copy;		// To facilitate copy operation
	Vector prob;
	int id;

	Doc(string line,int id=0);
	void emptyScore();
	~Doc();

	void samplezfromPrior();
	void sampleunif();
	void sampleunifheld();
	void allLevel2();
	double sampleLevels(int calcprob=0);
	double sampleLevelsheld(int calcprob = 0);
	double immitateLevels(vector<Doc>::iterator& d);
};

