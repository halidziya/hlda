#pragma once
#include "State.h"
#include <string>
#include <time.h>
#include <ppl.h>


using namespace std;
using namespace concurrency;

class Sampler
{
public:
	State state;									// holds the variables of gibbs sampler
	State beststate;
	State splitlaunch;
	State mergelaunch;
	//Parallel part
	State mergelaunches[8];
	State splitlaunches[8];
	State copies[8];

	Sampler();									// Default values
	Sampler(string filename);					// Data
	Sampler(string datafile, string confile);	// Data and configuration
	Sampler(string datafile, string testfile , string confile);	// Data and configuration

	double samplePath(vector<Doc>::iterator d, State& state, int init = 0,int heldout=0);
	double immitatePath(vector<Doc>::iterator d, vector<Doc>::iterator copy, State& state);
	double samplePath(vector<Doc>::iterator d, int init = 0,int heldout=0);
	double heldoutScore();
	~Sampler();
	
	void initialize();
	void gibbs();
	void mh();
	void parallelmh();



};

