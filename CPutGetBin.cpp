#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include "CPutGetBin.h"

CPutGetBin::CPutGetBin(int _id, int n_TotalSamples, int _capacityPut, int _capacityGet, string _grandPrefix)
{
	id = _id; 
	nSamplesGeneratedByFar = n_TotalSamples; 
	
	capacityPut = _capacityPut; 
	nPutUsed = 0;
	dataPut.resize(capacityPut); 	
	
	capacityGet = _capacityGet; 
	nGetUsed = capacityGet;
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
	convert << id << "." << n;	 
	string file_name = filename_prefix + convert.str(); 
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
		return false; 
	for (int i=0; i<nPutUsed; i++)
		write(oFile, &(dataPut[i])); 
	oFile.close();  
	return true; 
}

bool CPutGetBin::DrawSample(const gsl_rng *r, CSampleIDWeight &sample)
{
	int index; 
	if (nSamplesGeneratedByFar <= 0)
		return false; 
	if (GetNumberDataFile() <= 0) 
	/* when data have not been dumped to files
 	will directly get a data from dataPut
 	*/
	{
		index = gsl_rng_uniform_int(r, nSamplesGeneratedByFar); 
		sample = dataPut[index]; 
		return true; 
	}
	else 
	{
		if (nGetUsed == capacityGet)
		{
			Fetch(r);
			nGetUsed = 0; 
		}
		index = gsl_rng_uniform_int(r, capacityGet); 
		nGetUsed ++; 
		sample = dataGet[index]; 
		return true; 
	}
}

bool CPutGetBin::DrawSample(double *x, int dim, int &id, double &weight, const gsl_rng *r)
{
	int index; 
	if (nSamplesGeneratedByFar <= 0)
		return false; 
	if (GetNumberDataFile() <= 0)
	{
		index = gsl_rng_uniform_int(r, nSamplesGeneratedByFar);
		dataPut[index].GetData(x,dim,id, weight); 
		return true; 
	}
	else 
	{
		if (nGetUsed == capacityGet)
		{
			Fetch(r);
			nGetUsed = 0; 
		} 
		index = gsl_rng_uniform_int(r, capacityGet);  
		nGetUsed ++; 
		dataGet[index].GetData(x, dim, id, weight); 
		return true; 
	}
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
	convert << id << "." << i;
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

vector <CSampleIDWeight> CPutGetBin::RetrieveSamplesSequentially(bool if_clear_old_bin)
{
	vector <CSampleIDWeight> samples; 
	samples.clear(); 
	if (nPutUsed == 0 && GetNumberDataFile() > 0)
	{
		int index = GetNumberDataFile(); 
	        fstream iFile;
        	string file_name;
        	stringstream convert;

        	convert.str(std::string());
        	convert << id << "." << index;
        	file_name = filename_prefix + convert.str();
        	iFile.open(file_name.c_str(), ios::in|ios::binary);
        	if (!iFile)
                	return samples;
		while (nPutUsed < capacityPut)
		{
			read(iFile, &(dataGet[nPutUsed]));	
			nPutUsed ++;
		}
		iFile.close();
		remove(file_name.c_str()); 
	}
	if (nPutUsed)
	{
		samples.resize(nPutUsed); 
		for (int i=0; i<nPutUsed; i++)
			samples[i] = dataPut[i]; 
		nSamplesGeneratedByFar -= nPutUsed;
		nPutUsed = 0; 
	}
	return samples; 
}

void CPutGetBin::ChangeFileName(int new_id)
{
	stringstream convert; 
	string old_file, new_file; 
	for (int i=1; i<=GetNumberDataFile(); i++)
	{
		convert.str(std::string()); 
		convert << id << "." << i; 
		old_file = filename_prefix + convert.str(); 

		convert.str(std::string()); 
		convert << new_id << "." << i; 
		new_file = filename_prefix + convert.str(); 
		
		rename(old_file.c_str(), new_file.c_str()); 
	}
}

void CPutGetBin::ClearDepositDrawHistory()
{
	nSamplesGeneratedByFar = 0; 
 	nGetUsed = capacityGet;
	nPutUsed = 0; 
}
