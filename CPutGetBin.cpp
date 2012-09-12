#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include <glob.h>
#include "CPutGetBin.h"

CPutGetBin::CPutGetBin(int _id, int _nDumpFile, int _capacityPut, int _capacityGet, string _grandPrefix)
{
	id = _id; 
	nDumpFile = _nDumpFile; 
	
	capacityPut = _capacityPut; 
	nPutUsed = 0;
	// dataPut.resize(capacityPut); 	
	
	capacityGet = _capacityGet; 
	nGetUsed = capacityGet;
 	// dataGet.resize(capacityGet); 

	filename_prefix = _grandPrefix; 
} 

CPutGetBin::~CPutGetBin()
{
}

void CPutGetBin::SetCapacity_Put(int _capacityPut)
{
	capacityPut = _capacityPut; 
	// if ((int)(dataPut.size()) != capacityPut)
	//	dataPut.resize(capacityPut);
}

void CPutGetBin::SetCapacity_Get(int _capacityGet)
{
	capacityGet = _capacityGet; 
	//if ((int)(dataGet.size()) != capacityGet)
	//	dataGet.resize(capacityGet); 
}

int CPutGetBin::DepositSample(const CSampleIDWeight &sample)
{
	int index =  nPutUsed;
	if ((int)(dataPut.size()) < capacityPut)
		dataPut.push_back(sample); 
	else
		dataPut[index] = sample; 
	nPutUsed ++; 
	if (nPutUsed == capacityPut)
	{
		Dump(); 
		nDumpFile ++; 
		nPutUsed = 0; 
	}
	
	return nDumpFile*capacityPut+nPutUsed; 
}

int CPutGetBin::DepositSample(const double *x, int x_d, int x_index, double x_weight)
{
	int index = nPutUsed; 
	if ((int)(dataPut.size()) < capacityPut)
		dataPut.push_back(CSampleIDWeight(x, x_d, x_index, x_weight)); 
	else 
		dataPut[index] = CSampleIDWeight(x, x_d, x_index, x_weight); 
	nPutUsed ++; 
	if (nPutUsed == capacityPut)
	{
		Dump(); 
		nDumpFile++; 
		nPutUsed = 0; 
	}
	return nDumpFile*capacityPut+nPutUsed;
}

bool CPutGetBin::Dump()
{
	string file_name = GetFileNameForDump(); 
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
		return false; 
	for (int i=0; i<nPutUsed; i++)
		write(oFile, &(dataPut[i])); 
	oFile.flush(); 
	oFile.close();  
	return true; 
}

bool CPutGetBin::DrawSample(const gsl_rng *r, CSampleIDWeight &sample)
{
	int index; 
	if (nDumpFile*capacityPut+nPutUsed<= 0)
		return false; 
	else if (nPutUsed > 0) 
	/* when data have not been dumped to files
 	will directly get a data from dataPut
 	*/
	{
		index = gsl_rng_uniform_int(r, nPutUsed); 
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
	if (nDumpFile+capacityPut+nPutUsed<= 0)
		return false; 
	else if (nPutUsed<= 0)
	{
		index = gsl_rng_uniform_int(r, nPutUsed);
		dataPut[index].CopyData(x,dim,id, weight); 
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
		dataGet[index].CopyData(x, dim, id, weight); 
		return true; 
	}
}

bool CPutGetBin::Fetch(const gsl_rng *r)
{
	if (capacityGet > nDumpFile*capacityPut+nPutUsed)
		return false;

	if ((int)(dataGet.size()) < capacityGet)
		dataGet.resize(capacityGet);  
	
	int select; 
	vector < vector <int > > select_per_file(nDumpFile+1); 
	// First GetNumberDataFile(): files
	// Last: cache
	for (int i=0; i<capacityGet; i++)
	{
		select=gsl_rng_uniform_int(r, nDumpFile*capacityPut+nPutUsed); 
		select_per_file[select/capacityPut].push_back(select%capacityPut); 
	}
	int counter =0; 
	for (int i=0; i<nDumpFile; i++)	// Read data from file
	{
		if (!select_per_file[i].empty())
		{
			if (!ReadFromOneFile(i, counter, select_per_file[i]) )
				return false; 
		}
	}
	if (!select_per_file[nDumpFile].empty() ) // Read data from cache
	{
		for (int n=0; n<(int)(select_per_file[nDumpFile].size()); n++)
		{
			dataGet[counter] = dataPut[select_per_file[nDumpFile][n]]; 
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
	// determine the dimension of CSampleIDWeight
	CSampleIDWeight temp_data; 
	read(iFile, &temp_data); 

        for (int n=0; n<(int)(index.size()); n++)
        {
        	iFile.seekg(temp_data.GetSize_Data()*index[n], ios_base::beg);
                read(iFile, &(dataGet[counter]));
                counter ++;
        }
      	iFile.close();
	return true; 
}

bool CPutGetBin::restore()
{
	nDumpFile = GetNumberFileForDump(); 
	if (nDumpFile > 0)
	{
		fstream iFile;
       	 	string file_name;
        	stringstream convert;

        	convert.str("");
        	convert << id << "." << nDumpFile-1; 
        	file_name = filename_prefix + convert.str();
		iFile.open(file_name.c_str(), ios::in|ios::binary);
        	if (!iFile)
                	return false;
		
		// To determine size of each record
		CSampleIDWeight temp; 
		read(iFile, &temp); 

		iFile.seekg(0, ios::beg); 
		iFile.seekg(0, ios::end); 
		int lenFile = iFile.tellg(); 
		iFile.seekg(0, ios::beg); 

		nPutUsed = lenFile/temp.GetSize_Data(); 

		if (nPutUsed < capacityPut)
		{
			if ((int)dataPut.size() < nPutUsed)
				dataPut.resize(nPutUsed); 
			for (int n=0; n<nPutUsed; n++)
				read(iFile, &(dataPut[n])); 
			nDumpFile --; 
		}
		else 
			nPutUsed =0; 
		iFile.close(); 
	}
	return true; 
}

void CPutGetBin::finalize()
{
	if (nPutUsed > 0)
	{
		Dump();
		nDumpFile++;
		nPutUsed =0; 
	}  
}
 
void CPutGetBin::DisregardHistorySamples()
{
	string file_name; 
	stringstream convert; 
	
	for (int iFile =0; iFile<nDumpFile; iFile++)
	{
		convert.str(std::string());
                convert << id << "." << iFile; 
		file_name = filename_prefix + convert.str();	
		remove(file_name.c_str());
	}
}

void CPutGetBin::ClearDepositDrawHistory()
{
	nDumpFile = 0; 
 	nGetUsed = capacityGet;
	nPutUsed = 0; 
}

string CPutGetBin::GetFileNameForDump() const
{
	stringstream convert;
        convert << id << "." << nDumpFile;
        string file_name = filename_prefix + convert.str();
	return file_name; 
}

bool CPutGetBin::if_fetchable() 
{
	if (nDumpFile*capacityPut+nPutUsed > 0)
		return true; 
	else 
		return false; 
}

int CPutGetBin::GetNumberFileForDump() const
{
	stringstream convert; 
	convert.str(""); 
	convert << id << ".*"; // << "." << cluster_node; 

	string filename_pattern = filename_prefix + convert.str(); 
	
	glob_t glob_result; 
	glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result); 
	int final_result = (int)(glob_result.gl_pathc); 
	globfree(&glob_result); 
	return final_result; 
}
