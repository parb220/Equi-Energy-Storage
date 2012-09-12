#include <fstream>
#include <sstream>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <cstdio>
#include <glob.h>
#include "CPutGetBin.h"

CPutGetBin::CPutGetBin(int _id, int _nDumpFile, int _capacityPut, int _capacityGet, string _grandPrefix, int _suffix )
{
	suffix = _suffix; 
	id = _id; 
	nDumpFile = _nDumpFile; 
	
	capacityPut = _capacityPut; 
	nPutUsed = 0;
	// dataPut.resize(capacityPut); 	
	
	capacityGet = _capacityGet; 
	nGetUsed = capacityGet;
 	// dataGet.resize(capacityGet); 

	filename_prefix = _grandPrefix; 
	fetch_status = false; 
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
	vector <string> filename_fetch = GetFileNameForFetch(); 
	int nFetchFile = (int)(filename_fetch.size()); 
	if (nFetchFile*capacityPut+nPutUsed<= 0)
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
			Fetch(r, filename_fetch);
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
	vector <string> filename_fetch = GetFileNameForFetch(); 
	int nFetchFile = (int)(filename_fetch.size()); 
	if (nFetchFile*capacityPut+nPutUsed<= 0)
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
			Fetch(r, filename_fetch);
			nGetUsed = 0; 
		} 
		index = gsl_rng_uniform_int(r, capacityGet);  
		nGetUsed ++; 
		dataGet[index].CopyData(x, dim, id, weight); 
		return true; 
	}
}

bool CPutGetBin::Fetch(const gsl_rng *r, const vector<string> & filename_fetch)
{
	int nFetchFile = (int)(filename_fetch.size()); 
	if (capacityGet > nFetchFile*capacityPut+nPutUsed)
		return false;

	if ((int)(dataGet.size()) < capacityGet)
		dataGet.resize(capacityGet);  
	
	int select; 
	vector < vector <int > > select_per_file(nFetchFile+1); 
	// First nFetchFile(): files
	// Last: cache
	for (int i=0; i<capacityGet; i++)
	{
		select=gsl_rng_uniform_int(r, nFetchFile*capacityPut+nPutUsed); 
		select_per_file[select/capacityPut].push_back(select%capacityPut); 
	}
	int counter =0; 
	for (int i=0; i<nFetchFile; i++)	// Read data from file
	{
		if (!select_per_file[i].empty())
		{
			if (!ReadFromOneFile(filename_fetch[i], counter, select_per_file[i]) )
				return false; 
		}
	}
	if (!select_per_file[nFetchFile].empty() ) // Read data from cache
	{
		for (int n=0; n<(int)(select_per_file[nFetchFile].size()); n++)
		{
			dataGet[counter] = dataPut[select_per_file[nFetchFile][n]]; 
			counter ++;
		}
	}
	return true; 
}

bool CPutGetBin::ReadFromOneFile(string file_name, int &counter, const vector <int> &index) 
{	
	fstream iFile; 
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
        	convert << id << "." << nDumpFile-1 << "." << suffix; 
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
	vector <string> filename_fetch = GetFileNameForFetch(); 
	int nFetchFile = (int)(filename_fetch.size()); 
		
	for (int iFile =0; iFile<nFetchFile; iFile++)
		remove(filename_fetch[iFile].c_str());
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
        convert << id << "." << nDumpFile << "." << suffix;
        string file_name = filename_prefix + convert.str();
	return file_name; 
}

bool CPutGetBin::if_fetchable() 
{
	if (!fetch_status)
	{
		int nFetchFile = GetNumberFileForFetch(); 
		if (nFetchFile*capacityPut+nPutUsed > 0)
			fetch_status = true; 
		else 
			fetch_status = false; 
	} 
	return fetch_status; 
}

int CPutGetBin::GetNumberFileForDump() const
{
	stringstream convert; 
	convert.str(""); 
	convert << id << ".*." << suffix; // << "." << cluster_node; 

	string filename_pattern = filename_prefix + convert.str(); 
	
	glob_t glob_result; 
	glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result); 
	int final_result = (int)(glob_result.gl_pathc); 
	globfree(&glob_result); 
	return final_result; 
}

vector <string> CPutGetBin::GetFileNameForFetch() const
{
	stringstream convert;
        convert.str("");
        convert << id << ".*.*"; // << ".*" // all cluster node; 

        string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	vector <string> filename_fetch; 
	if (glob_result.gl_pathc > 0)
	{
		filename_fetch.resize(glob_result.gl_pathc); 
		for (int i=0; i<(int)(glob_result.gl_pathc); i++)
			filename_fetch[i] = string(glob_result.gl_pathv[i]); 
	}
	else 
		filename_fetch.clear(); 
        globfree(&glob_result);
        return filename_fetch; 
}

int CPutGetBin::GetNumberFileForFetch() const
{
	vector <string> file_fetch = GetFileNameForFetch();
	return (int)(file_fetch.size()); 
}
