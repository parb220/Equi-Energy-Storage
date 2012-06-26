#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "CPutGetBin.h"

CPutGetBin::CPutGetBin(int _id, int n_TotalSamples, int _capacityPut, int _capacityGet, string _grandPrefix)
{
	id = _id; 
	nSamplesGeneratedByFar = n_TotalSamples; 
	
	capacityPut = _capacityPut; 
	nPutUsed = 0;
	dataPut.resize(capacityPut); 	
	
	capacityGet = _capacityGet; 
	nGetUsed = 0;
 	dataGet.resize(capacityGet); 

	filename_prefix = _grandPrefix; 
} 

CPutGetBin::~CPutGetBin()
{
}

void CPutGetBin::SetCapacity_Put(int _capacityPut)
{
	capacityPut = _capacityPut; 
	if ((int)(dataPut.size()) != capacityPut)
		dataPut.resize(capacityPut);
}

void CPutGetBin::SetCapacity_Get(int _capacityGet)
{
	capacityGet = _capacityGet; 
	if ((int)(dataGet.size()) != capacityGet)
		dataGet.resize(capacityGet); 
}

int CPutGetBin::DepositSample(const CSampleIDWeight &sample)
{
	int index =  nPutUsed; 
	dataPut[index] = sample; 
	nSamplesGeneratedByFar ++; 
	nPutUsed ++; 
	if (nPutUsed == capacityPut)
	{
		Dump(GetNumberDataFile()); 
		nPutUsed = 0; 
	}
	
	return nSamplesGeneratedByFar; 
}

int CPutGetBin::DepositSample(const double *x, int x_d, int x_index, double x_weight)
{
	int index = nPutUsed; 
	dataPut[index] = CSampleIDWeight(x, x_index, x_weight); 
	nSamplesGeneratedByFar ++; 
	nPutUsed ++; 
	if (nPutUsed == capacityPut)
	{
		Dump(GetNumberDataFile()); 
		nPutUsed = 0; 
	}
	return nSamplesGeneratedByFar; 
}

bool CPutGetBin::Dump(int n)
{
	stringstream convert;
	convert << "." << id << "." << n;	 
	string file_name = filename_prefix + convert.str(); 
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
		return false; 
	for (int i=0; i<nPutUsed; i++)
		write(oFile, &(dataPut[i])); 
	oFile.close();  
	return true; 
}

CSampleIDWeight CPutGetBin::DrawSample(const gsl_rng *r)
{
	int index; 
	if (GetNumberDataFile() <= 0) 
	/* when data have not been dumped to files
 	will directly get a data from dataPut
 	*/
		index = gsl_rng_uniform_int(r, nSamplesGeneratedByFar); 
	else 
	{
		index = gsl_rng_uniform_int(r, capacityGet); 
		nGetUsed ++; 
		if (nGetUsed == capacityGet)
		{
			Fetch(r);
			nGetUsed = 0; 
		}
	}
	return dataGet[index]; 
}

void CPutGetBin::DrawSample(double *x, int dim, int &id, double &weight, const gsl_rng *r)
{
	int index; 
	if (GetNumberDataFile() <= 0)
		index = gsl_rng_uniform_int(r, nSamplesGeneratedByFar);
	else 
	{
		index = gsl_rng_uniform_int(r, capacityGet);  
		nGetUsed ++; 
		if (nGetUsed == capacityGet)
		{
			Fetch(r);
			nGetUsed = 0; 
		} 
	}
	dataGet[index].GetData(x, dim, id, weight); 
}

bool CPutGetBin::Fetch(const gsl_rng *r)
{
	if (capacityGet > nSamplesGeneratedByFar)
		return false; 
	
	int select; 
	vector < vector <int > > select_per_file(GetNumberDataFile()+1); 
	// First GetNumberDataFile(): files
	// Last: cache
	for (int i=0; i<capacityGet; i++)
	{
		select=gsl_rng_uniform_int(r, nSamplesGeneratedByFar); 
		select_per_file[select/capacityPut].push_back(select%capacityPut); 
	}
	int counter =0; 
	for (int i=0; i<GetNumberDataFile(); i++)	// Read data from file
	{
		if (!select_per_file[i].empty())
		{
			if (!ReadFromOneFile(i+1, counter, select_per_file[i]) )
				return false; 
		}
	}
	if (!select_per_file[GetNumberDataFile()].empty() ) // Read data from cache
	{
		for (int n=0; n<(int)(select_per_file[GetNumberDataFile()].size()); n++)
		{
			dataGet[counter] = dataPut[select_per_file[GetNumberDataFile()][n]]; 
			counter ++; 
		}
	}
	return true; 
}

bool CPutGetBin::ReadFromOneFile(int i, int &counter, const vector <int> &index) 
{	
	fstream iFile; 
	string file_name; 
	stringstream convert;

	convert.str(std::string());
	convert << "." << id << "." << i;
        file_name = filename_prefix + convert.str();
        iFile.open(file_name.c_str(), ios::in|ios::binary);
        if (!iFile)
        	return false;
        for (int n=0; n<(int)(index.size()); n++)
        {
        	iFile.seekg(CSampleIDWeight::GetSize_Data()*index[n], ios_base::beg);
                read(iFile, &(dataGet[counter]));
                counter ++;
        }
      	iFile.close();
	return true; 
}

void CPutGetBin::finalize()
{
	if (nPutUsed > 0)
		Dump(GetNumberDataFile()+1);  
}
 
int CPutGetBin::GetNumberDataFile(bool flag) const
{
	if (flag)
		return ceil((double)nSamplesGeneratedByFar/capacityPut);  
	else 
		return GetNumberDataFile(); 
}
