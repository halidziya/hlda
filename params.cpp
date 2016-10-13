#include "params.h"

int depth;				 // Depth of the tree
int nterms;
int niter;
int inititer;
int smiter;
int nheldout;
int lgammalimit;

Vector eta;		 // Topic Dirichlet Param
Vector etasum;	 // eta*nterms
Vector gamma;	 // nCRP param
Vector alpha;	 // Local LDA document Dirichlet param

Table gammatable;

double  Table::get(double x, int y)
{ 
	if (y < lgammalimit && (data.find(x) != data.end()))
		return data[x][y]; 
	else 
		return lgamma(x + y); 
}

Vector& Table::operator[](double s){
	return data[s]; 
};