#pragma once
#include <vector>
#include <map>
#include "Vector.h"

using namespace std;

extern int depth;				 // Depth of the tree
extern int nterms;
extern int niter;
extern int inititer;
extern int smiter;
extern int nheldout;
extern int lgammalimit;

extern Vector eta;		 // Topic Dirichlet Param
extern Vector etasum;	 // eta*nterms
extern Vector gamma;	 // nCRP param
extern Vector alpha;	 // Local LDA document Dirichlet param


class Table{
	map<double, Vector> data;
public:
	Vector& operator[](double s);
	double  get(double x, int y);
};

extern Table gammatable; // Holds the precomputed gamma function values. Be careful about floating point indexing