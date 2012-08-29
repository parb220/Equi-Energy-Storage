#include <cstring>
#include <iostream>
#include <fstream>
#include "CSampleIDWeight.h"

using namespace std; 

void CSampleIDWeight::SetDataDimension(int _dim)
{
	if (dim != _dim)
	{
		if (dim > 0)
			delete [] data; 
		dim = _dim; 
		data = new double[dim]; 
	}
}

CSampleIDWeight::CSampleIDWeight()
{
	dim = 0; 
	data = NULL; 
	id = 0; 
	weight = 0; 
	if_weight_set = false; 
}

CSampleIDWeight::CSampleIDWeight(const double *x, int _dim, int _id, double _weight)
{
	dim = _dim; 
	data = new double[dim]; 
	memcpy(data, x, sizeof(double)*dim); 
	id = _id; 
	weight = _weight; 
	if_weight_set = true; 
}

/*CSampleIDWeight::CSampleIDWeight(const vector <double> &x, int _id, double _weight)
{
	data = x; 
	id = _id; 
	weight = _weight; 
}*/

CSampleIDWeight::CSampleIDWeight(const CSampleIDWeight &right)
{
	dim = right.dim; 
	data = new double[dim]; 
	memcpy(data, right.data, sizeof(double)*dim);
	id = right.id; 
	weight = right.weight;
	if_weight_set = right.if_weight_set; 
}

CSampleIDWeight::~CSampleIDWeight()
{
	/* if (data != NULL)
		delete [] data; */
	if (dim > 0)
		delete [] data; 
}

CSampleIDWeight & CSampleIDWeight:: operator = (const CSampleIDWeight &right)
{	
	SetDataDimension(right.dim); 
	memcpy(data, right.data, sizeof(double)*dim); 
	id = right.id; 
	weight = right.weight;	
	if_weight_set = right.if_weight_set
	return *this;
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

	/*for (int i=0; i<CSampleIDWeight::dim; i++)
		input_stream >> x.data[i]; 
	input_stream >> x.id; 
	input_stream >> x.weight; 
	return input_stream;*/
} 

ostream& write(ostream & output_stream, const CSampleIDWeight *x)
{
	output_stream.write((char*)&(x->dim), sizeof(int)); 
	output_stream.write((char*)&(x->data[0]), sizeof(double)*x->dim); 
	output_stream.write((char*)&(x->id), sizeof(int)); 
	output_stream.write((char*)&(x->weight), sizeof(double)); 
	return output_stream; 

	/*for (int i=0; i<CSampleIDWeight::dim; i++)
		output_stream << x.data[i] << "\t"; 
	output_stream << x.id << "\t"; 
	output_stream << x.weight << "\n";
	return output_stream;*/
}

void CSampleIDWeight::GetData(double *_x, int _dim, int &_id, double &_weight)
{
	memcpy(_x, data, sizeof(double)*dim); 
	/*for (int i=0; i<dim; i++)
		_x[i] = data[i]; */
	_id = id; 
	_weight = weight; 
}

/*void CSampleIDWeight::GetData(vector < double > &_x, int &_id, double &_weight)
{
	_x = data; 
	_id = id; 
	_weight = weight; 
}*/

int CSampleIDWeight::GetSize_Data()
{	
	return sizeof(int)+sizeof(double)*dim+sizeof(int)+sizeof(double); 
}

istream& operator >> (istream &inputStr, CSampleIDWeight &sample)
{
	for (int i=0; i<sample.dim; i++)
		inputStr >> sample.data[i]; 
	inputStr >> sample.id; 
	inputStr >> sample.weight; 
	return inputStr; 
}

ostream& operator << (ostream &outputStr, const CSampleIDWeight &sample)
{
	for (int i=0; i<sample.dim; i++)
		outputStr << sample.data[i] << "\t"; 
	outputStr << sample.id << "\t"; 
	outputStr << sample.weight << endl; 

	return outputStr; 
}
