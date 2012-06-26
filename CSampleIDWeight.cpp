#include <cstring>
#include <iostream>
#include <fstream>
#include "CSampleIDWeight.h"

using namespace std; 
int CSampleIDWeight::dim; 

void CSampleIDWeight::SetDataDimension(int _dim)
{
	dim = _dim; 
}

CSampleIDWeight::CSampleIDWeight()
{
	data = new double[dim]; 
	//data.resize(dim); 
	id = 0; 
	weight = 1.0; 
}

CSampleIDWeight::CSampleIDWeight(const double *x, int _id, double _weight)
{
	data = new double[dim]; 
	memcpy(data, x, sizeof(data)*dim); 
	// data.resize(dim); 
	for (int i=0; i<dim; i++)
		data[i] = x[i]; 
	id = _id; 
	weight = _weight; 
}

/*CSampleIDWeight::CSampleIDWeight(const vector <double> &x, int _id, double _weight)
{
	data = x; 
	id = _id; 
	weight = _weight; 
}*/

CSampleIDWeight::CSampleIDWeight(const CSampleIDWeight &right)
{
	data = new double[dim]; 
	memcpy(data, right.data, sizeof(data)*dim);
	// data = right.data; 
	id = right.id; 
	weight = right.weight;
}

CSampleIDWeight::~CSampleIDWeight()
{
	/* if (data != NULL)
		delete [] data; */
}

CSampleIDWeight & CSampleIDWeight:: operator = (const CSampleIDWeight &right)
{	
	if (data == NULL)
		data = new double[dim]; 
	memcpy(data, right.data, sizeof(data)*dim); 
	// data = right.data; 
	id = right.id; 
	weight = right.weight;	
	return *this;
}

istream & read (istream & input_stream, CSampleIDWeight *x)
{
	input_stream.read((char*)&(x->data[0]), sizeof(double)*CSampleIDWeight::dim); 
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
	output_stream.write((char*)&(x->data[0]), sizeof(double)*CSampleIDWeight::dim); 
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
	return sizeof(double)*dim+sizeof(int)+sizeof(double); 
}
