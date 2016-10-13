#pragma once
#include "Vector.h"
class Matrix :
	public Vector
{
public:
	int m;
	int r;

	Matrix(void);
	Matrix(int d);
	Matrix(int r,int m,int real=1);			// refcount
	Matrix(Matrix&& mat);			// Move
	Matrix(const Matrix& mat);			// Copy
									// 1 is standart , 2  is persistent , 0 is abstract
	~Matrix(void);
	
	void readMatrix(char* filename);
	void writeMatrix(char* filename);
	void readBin(char* filename);
	void writeBin(char* filename);
	

	Vector operator[](int i);   // Get row
	Vector operator()(int i);   // Get specific item

	friend ostream& operator<<(ostream& os, const Matrix& v);
};

