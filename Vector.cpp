#include "Vector.h"
#include "Matrix.h"
#include <new>

MultiBuffer<Vector> buffer;
MultiBuffer<Vector> absbuffer;

void init_buffer(int nthreads,int d)
{
	nthreads = nthreads + 1; // And one master thread
	new (&buffer) MultiBuffer<Vector>(nthreads, BUFF_SIZE, d, 2);
	new (&absbuffer) MultiBuffer<Vector>(nthreads, BUFF_SIZE, d, 0);
}

Vector::Vector(int size,int real):n(size)
{
	if (real && size > 0) // Buffer or real vector
		data = (double*) malloc(sizeof(double)*size);
	if (size == 0)
		type = 0; // Abstract
	else
		this->type = real;

	n = size;
}

Vector::Vector(std::initializer_list<double> l) : Vector(l.size()) {
	std::initializer_list<double>::iterator it=l.begin();
	for (int i = 0; i < n; i++)
	{
		data[i] = *it;
		it++;
	}
}

Vector::Vector():n(0){type=1;data=NULL;}

Vector::~Vector()
{
	if (type==1)
		free(data);
}

void Vector::zero()
{
	memset(data,0,sizeof(double)*n);
}


double Vector::operator*(Vector& v) // Dot product
{
	int i;
	double res=0;
	for(i=0;i<n;i++)
		res += data[i] * v.data[i];

	return res;
}


/* Return integer value */
int Vector::operator()(int i){
		return (int) (data[i]+0.5);
}


void Vector::print()
{
	int i;
	printf("\n[");
	for(i=0;i<n;i++)
		printf("%f ",data[i]);
	printf("]\n");
}

Vector::Vector(Vector&& v)
{
	n = v.n;
	data  = v.data;
	v.data = NULL;
	v.n = 0;
	v.type = 0; // No longer valid vector
	type = 1; // Owns the resource
}

Vector::Vector(const Vector& v)
{
	n = v.n;
	type = v.type;
	switch (type) // Real or Buffer
	{
	case 0:
		data = v.data;
		break;
	case 1: // Deep copy
		data  = (double*) malloc(sizeof(double)*n);
		memcpy(data,v.data,sizeof(double)*n);
		break;
	case 2: // Buffer is persistent keep a abstract vector to point
		type = 0;
		data = v.data;
		break;
	}
}



void Vector::operator=(const Vector& v)
{
	type = 1;
	switch (v.type) // Real or Buffer
	{
	case 0:
	case 1:
	case 2:
		if (!data) // if size is same
			data = (double*) malloc(sizeof(double)*v.n);
		else if (n!=v.n)
			data = (double*) realloc(data,v.n*sizeof(double));
		memcpy(data,v.data,sizeof(double)*v.n);
		break;
	}
	n = v.n;
}



Vector  Vector::copy()
{
	int i;
	Vector res(n);
	memcpy(res.data,data,n*sizeof(double));
	return res;
}


void Vector::resize(int size)
{
	if (n != size)
	{
		n = size;
		data = (double*)realloc(data, sizeof(double)*size);
	}
}

Vector Vector::operator*(double scalar)
{
	Vector& r = buffer.get();
	r.n = n;
	for(int i=0;i<n;i++)
		r.data[i]= data[i] * scalar;
	return buffer.next();
}

Vector Vector::operator-(Vector& v)
{
	Vector& r = buffer.get();
	r.n = n;
	for(int i=0;i<n;i++)
		r.data[i]= data[i] - v.data[i];
	return buffer.next();
}

Vector Vector::operator+(Vector& v)
{
	Vector& r = buffer.get();
	r.n = n;
	for(int i=0;i<n;i++)
		r.data[i]= data[i] + v.data[i];
	return buffer.next();
}


void Vector::operator<<(Vector& v)
{
	for (int i = 0; i<n; i++)
		data[i] +=  v.data[i];
}

void Vector::operator>>(Vector& v)
{
	for (int i = 0; i<n; i++)
		data[i] = data[i] - v.data[i];
}

Vector Vector::operator/(double scalar)
{
	Vector& r = buffer.get();
	r.n = n;
	double divval = 1.0/scalar;
	for(int i=0;i<n;i++)
		r.data[i]= data[i] * divval;
	return buffer.next();
}


double Vector::maximum()
{
	double val  = data[0];
	for(int i=1;i<n;i++)
		if (val<data[i])
			val = data[i];

	return val;
}

double Vector::sum()
{
	double val = 0;
	for (int i = 0; i < n; i++)
		val += data[i];
	return val;
}

int Vector::sample() // Assuming value is normalized
{
	double s = sum();
	double u = urand()*s;
	int i = 0;
	while (u > data[i])
	{
		u = u - data[i];
		i++;
	}
	return i;
}

ostream& operator<<(ostream& os, const Vector& v)
{   // It does not save type , it save as real vector always
	int one = 1;
	os.write((char*) &v.n,sizeof(int)); 
	os.write((char*) &one,sizeof(int)); 
	os.write((char*) v.data,sizeof(double)*v.n); 

	return os;
}