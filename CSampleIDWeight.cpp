#include <cstring>
#include <iostream>
#include <fstream>
#include "CSampleIDWeight.h"

using namespace std; 

CSampleIDWeight::CSampleIDWeight() : dim(0), data(NULL), id(0), weight(0.0), log_prob(0.0)
{}

CSampleIDWeight::CSampleIDWeight(const double *x, int _dim, int _id, double _weight)
	: dim(_dim), data(NULL), id(_id), weight(_weight), log_prob(0.0)
{
	if (dim > 0)
	{
		data = new double[dim]; 
		memcpy(data, x, sizeof(double)*dim); 
	}
}

CSampleIDWeight::CSampleIDWeight(const CSampleIDWeight &right)
	: dim(right.dim), data(NULL), id(right.id), weight(right.weight), log_prob(right.log_prob)
{
	if (dim > 0)
	{
		data = new double[dim]; 
		memcpy(data, right.data, sizeof(double)*dim);
	}
}

CSampleIDWeight::~CSampleIDWeight()
{
	if (dim > 0)
		delete [] data; 
}

void CSampleIDWeight::SetDataDimension(int _dim)
{
	if (dim != _dim)
	{
		if (dim > 0)
			delete [] data; 
		dim = _dim; 
		if (dim > 0)
			data = new double[dim];
		else 
			data = NULL;  
	}
}

CSampleIDWeight & CSampleIDWeight:: operator = (const CSampleIDWeight &right)
{	
	SetDataDimension(right.dim); 
	memcpy(data, right.data, sizeof(double)*dim); 
	id = right.id; 
	weight = right.weight;	
	log_prob = right.log_prob; 
	return *this;
}

void CSampleIDWeight::Add(const CSampleIDWeight &right)
{
	for (int i=0; i<dim; i++)
		data[i] = data[i] + right.data[i];  
}

void CSampleIDWeight::Subtract (const CSampleIDWeight &right)
{
	for (int i=0; i<dim; i++)
		data[i] = data[i] - right.data[i]; 
}


void CSampleIDWeight:: PartialCopyFrom(const CSampleIDWeight &right, int offset, int length)
{
	memcpy(data+offset, right.data+offset, sizeof(double)*length); 
	if (offset ==0 && length == dim && length == right.dim)
	{
		weight = right. weight; 
		log_prob = right.log_prob; 
	}
}

void CSampleIDWeight:: PartialCopyFrom(int offset1, const CSampleIDWeight &right, int offset2, int length)
{
	memcpy(data+offset1, right.data+offset2, sizeof(double)*length); 
	if (offset1 == 0 && offset2 == 0 && length == dim && length == right.dim)
	{
		weight = right.weight;
		log_prob = right.log_prob; 
	}
}

istream & read (istream & input_stream, CSampleIDWeight *x)
{
	int dim; 
	input_stream.read((char*)&(dim), sizeof(int)); 
	x->SetDataDimension(dim);
	input_stream.read((char*)&(x->data[0]), sizeof(double)*x->dim); 
	input_stream.read((char*)&(x->id), sizeof(int)); 
	input_stream.read((char*)&(x->weight), sizeof(double)); 
	return input_stream; 
}

ostream& write(ostream & output_stream, const CSampleIDWeight *x)
{
	output_stream.write((char*)&(x->dim), sizeof(int)); 
	output_stream.write((char*)&(x->data[0]), sizeof(double)*x->dim); 
	output_stream.write((char*)&(x->id), sizeof(int)); 
	output_stream.write((char*)&(x->weight), sizeof(double)); 
	return output_stream; 
}

void CSampleIDWeight::CopyData(double *_x, int _dim, int &_id, double &_weight) const
{
	memcpy(_x, data, sizeof(double)*dim); 
	_id = id; 
	_weight = weight; 
}

void CSampleIDWeight::CopyData(double *_x, int _dim) const
{
	memcpy(_x, data, sizeof(double)*dim); 
}

int CSampleIDWeight::GetSize_Data() const
{	
	return sizeof(int)+sizeof(double)*dim+sizeof(int)+sizeof(double); 
}

istream& operator >> (istream &inputStr, CSampleIDWeight &sample)
{
	int dim; 
	inputStr >> dim; 
	sample.SetDataDimension(dim); 
	for (int i=0; i<sample.dim; i++)
		inputStr >> sample.data[i]; 
	inputStr >> sample.id; 
	inputStr >> sample.weight; 
	return inputStr; 
}

ostream& operator << (ostream &outputStr, const CSampleIDWeight &sample)
{
	outputStr << sample.dim << "\t"; 
	for (int i=0; i<sample.dim; i++)
		outputStr << sample.data[i] << "\t"; 
	outputStr << sample.id << "\t"; 
	outputStr << sample.weight << endl; 

	return outputStr; 
}
