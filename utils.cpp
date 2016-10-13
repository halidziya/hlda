#include "utils.h"
/*
using namespace std;
default_random_engine generator;


double urand()
{
	uniform_real_distribution<double> distribution(0.0, 1.0);
	return distribution(generator);
}*/

//#include <random>
//random_device rd;

double urand()
{
	return  randomMT() / (double(0xffffffff)+1); // ((double) rd() )/ rd.max();//  // randomMT() / (1.0 + (0xffffffff));  //  / ((double)(RAND_MAX + 1.0));
}