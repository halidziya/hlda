#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ostream>
#include "DebugUtils.h"
#include <thread>
#include "Buffer.h"
#include "utils.h"
#include <initializer_list>
#define EPS 1e-10 // Floating point equality

using namespace std;


/*	Vector Types :
	TYPE 0 : Abstract	vector , just a pointer
	TYPE 1 : Standart	vector , allocates deallocates memory automatically
	TYPE 2 : Buffer		vector , calculations are used on top.

*/

class Vector{
public:
	double* data;
	int n;
	int type;
	Vector(std::initializer_list<double> l);
	Vector(int size,int real=1);
	~Vector(void);
	
	Vector (const Vector& v);
	Vector (Vector&& v);
	Vector();
	

	void print();
	void zero();
	void resize(int size);
	Vector copy();		// Returns real vector
	Vector unique();	// Returns unique values 


	inline double& operator[](const int i){ return data[i]; };
	int		operator()(int i);	 // Return integer
	double operator*(Vector& v); // Inner product
	Vector operator-(Vector& v); // Subtraction
	void operator=(const Vector& v); // Assignment 
	Vector operator*(double scalar); // Scaling
	Vector operator+(Vector& v);	// Summation
	void operator<<(Vector& v);	// Add on it
	void operator>>(Vector& v);	// Subtract from it
	Vector operator/(double scalar); // Scaling
	double maximum();
	double sum();
	int sample();
	friend ostream& operator<<(ostream& os, const Vector& v);
};


// BE CAREFUL ONCE BUFFER IS INITIALIZED IT OPERATES IN 'd' DIMENSIONS
#define BUFF_SIZE 10
extern MultiBuffer<Vector> buffer;
extern MultiBuffer<Vector> absbuffer;

void init_buffer(int nthreads,int d);