#pragma once
#include <vector>
#include "hldaTree.h"
#include <fstream>
#include <sstream>
#include <iostream>



using namespace std;


// State should be saved and overwrited 

class State
{
public:

	hldaTree	   tree;	 // hLDA Tree
	vector<Doc>    docs;     // Documents
	vector<Doc>    heldoutdocs;     // Documents
	int ndocs;				 // Number of documents in the corpus 

	State();				 // Default values
	State(string filename);	 // Data
	State(string datafile, string confile, string testfile);	 // Data and configuration
	State(string datafile,string confile);	 // Data and configuration
	~State();

	void readDocs(string filename, vector<Doc>& docs);		// Read data
	void readConf(string confile);		// Read config/parameters
	double likelihood();						// Likelihood
	double score;						// Calculated score

	void operator=(State& s);
	void write(string filename);
	double etaScore(); // Topic likelihoods
	double crpScore(); // nCRP prior
	double alphaScore(); // Level allocation dirichlet prior
	double heldlikelihood();
private:

	
};