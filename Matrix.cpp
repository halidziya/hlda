#include "Matrix.h"
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

Matrix::Matrix(void)
{
	r = 0;
	m = 0;
}



Matrix::~Matrix(void)
{
	// Vector destructor takes care of it
}

Matrix::Matrix(int r,int m,int real) : Vector(r*m,real)
{
	this->r = r;
	this->m = m;
}


Matrix::Matrix(int r) : Vector(r*r)
{
	this->r = r;
	this->m = r;
}

void Matrix::readMatrix(char* filename)
{
	ifstream file(filename);
	int i,j;
	double val;
	if (!file.is_open())
		printf("Could not operfile...\n");
	else
	{
		printf("Reading %s...\n",filename);
		file >> r >> m;
		data = (double*) malloc(sizeof(double)*r*m);
		type = 1;
		for(i=0;i<r;i++)
			for(j=0;j<m;j++)
				file >> data[i*m+j];
		file.close();
	}
}

void Matrix::writeMatrix(char* filename)
{
	ofstream file(filename);
	int i,j;
	file << r << " " << m << endl;
	for (i=0;i<r;i++)
	{
		for(j=0;j<m;j++)
			file << data[i*m+j] << " ";
		file << endl;
	}
	file.close();
}


void Matrix::readBin(char* filename)
{
	FILE* file;
	fopen_s(&file,filename,"rb");
	int i,j;
	double val;
	if (file==NULL)
		printf("Could not open file...\n");
	else
	{
		printf("Reading %s...\n",filename);
		fread(&r,sizeof(int),1,file);
		fread(&m,sizeof(int),1,file);
		n = r*m;
		type = 1;
		data  = (double*) malloc(sizeof(double)*n);
		fread(data,sizeof(double),n,file);
		fclose(file);
	}
}

void Matrix::writeBin(char* filename)
{
	FILE* file;
	fopen_s(&file,filename,"wb");
	int i,j;
	fwrite(&r,sizeof(int),1,file);
	fwrite(&m,sizeof(int),1,file);

	n = r*m;

	fwrite(data,sizeof(double),n,file);
	fclose(file);
}


/* Get Row */
Vector Matrix::operator[](int i){ 
	absbuffer.get().data = data + m*i;
	return absbuffer.next();
}

Vector Matrix::operator()(int i)
{
	return Vector::operator()(i);
}






Matrix::Matrix(Matrix&& mat) : Vector(mat)
{
	m = mat.m;
	r = mat.r;
}

Matrix::Matrix(const Matrix& mat) : Vector(mat) , r(mat.r) , m(mat.m)
{
}


ostream& operator<<(ostream& os, const Matrix& v)
{   // It does not save type , it save as real vector always
	os.write((char*) &v.r,sizeof(int)); 
	os.write((char*) &v.m,sizeof(int)); 
	os.write((char*) v.data,sizeof(double)*v.n); 
	// os.write((char*) &v.triangle,sizeof(int));  // Not used greatly
	return os;
}