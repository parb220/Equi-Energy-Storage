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

bool CPutGetBin::Dump(string _filename)
{
	string file_name; 
	if (_filename.empty())
		file_name = GetFileNameForDump(); 
	else 
		file_name = _filename; 
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

void CPutGetBin::consolidate()
{
	vector < string> file_consolidate = GetFileNameForConsolidate(); 
	if (file_consolidate.size() >= 2)
	{
		vector <CSampleIDWeight> sample_consolidate; 
		vector <CSampleIDWeight> temp_sample; 
		for (int i=0; i<(int)(file_consolidate.size()); i++)
		{
			temp_sample = ReadSampleFromFile(file_consolidate[i]); 
			sample_consolidate.insert(sample_consolidate.end(), temp_sample.begin(), temp_sample.end()); 
			remove(file_consolidate[i].c_str()); 
		}
		int nComplete = (int)(sample_consolidate.size())/capacityPut; 
		int nRemaining = (int)(sample_consolidate.size())%capacityPut; 
		for (int iComplete=0; iComplete<nComplete; iComplete ++)
		{
			dataPut.resize(capacityPut); 
			for (int j=0; j<capacityPut; j++)
				dataPut[j] = sample_consolidate[j+iComplete*capacityPut]; 
			nPutUsed = capacityPut; 
			Dump(file_consolidate[iComplete]); 
		}
		if (nRemaining > 0)
		{
			if ((int)(dataPut.size()) < nRemaining)
				dataPut.resize(nRemaining); 
			for (int j=0; j<nRemaining; j++)
				dataPut[j] = sample_consolidate[j+nComplete*capacityPut]; 
			nPutUsed = nRemaining; 
			Dump(file_consolidate[nComplete]); 
		}
	}
}

void CPutGetBin::restore()
{
	nDumpFile = GetNumberFileForDump(); 
	if (nDumpFile > 0)
	{
       	 	string file_name;
        	stringstream convert;

        	convert.str(string());
        	convert << id << "." << nDumpFile-1 << "." << suffix; 
        	file_name = filename_prefix + convert.str();
		dataPut = ReadSampleFromFile(file_name); 
		nPutUsed = (int)(dataPut.size()); 
 
		if (nPutUsed >0 && nPutUsed < capacityPut)
			nDumpFile --; 
		else if (nPutUsed == capacityPut)
			nPutUsed =0; 
	}
}

void CPutGetBin::RestoreForFetch()
{
	vector <string> filename_partial = GetFileNameForConsolidate(); 
	// After consolidation, there should be at most 1 partial file
	if (!filename_partial.empty())
	{
		vector<CSampleIDWeight> tempSample = ReadSampleFromFile(filename_partial[0]); 
		if (!tempSample.empty())
		{
			nPutUsed = (int)(tempSample.size()); 
			if (dataPut.size() < tempSample.size())
				dataPut.resize(capacityPut); 
			for (int j=0; j<(int)(tempSample.size()); j++)
				dataPut[j] = tempSample[j]; 
		}
	}
}


vector < CSampleIDWeight> CPutGetBin::ReadSampleFromFile(string file_name) const
{
	int nRecord = NumberRecord(file_name); 
	if (nRecord == 0)
		return vector<CSampleIDWeight> (0);
	fstream iFile;
        iFile.open(file_name.c_str(), ios::in|ios::binary);
        if (!iFile)
        	return vector<CSampleIDWeight> (0); 

        vector <CSampleIDWeight> sample(nRecord); 
        for(int n=0; n<nRecord; n++)
		read(iFile, &(sample[n]));
        iFile.close();
	return sample; 
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
	convert.str(string()); 
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
        convert.str(string());
        convert << id << ".*.*"; // << ".*" // all cluster node; 

        string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	vector <string> filename_fetch; 
	if (glob_result.gl_pathc > 0)
	{
		for (int i=0; i<(int)(glob_result.gl_pathc); i++)
			if (NumberRecord(string(glob_result.gl_pathv[i])) == capacityPut)
				filename_fetch.push_back(string(glob_result.gl_pathv[i])); 
	}
	else 
		filename_fetch.clear(); 
        globfree(&glob_result);
        return filename_fetch; 
}

vector <string> CPutGetBin::GetFileNameForConsolidate() const
{
	stringstream convert;
        convert.str(string());
        convert << id << ".*.*"; // << ".*" // all cluster node; 

        string filename_pattern = filename_prefix + convert.str();

        glob_t glob_result;
        glob(filename_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
	vector <string> filename_consolidate; 
	if (glob_result.gl_pathc > 0)
	{
		for (int i=0; i<(int)(glob_result.gl_pathc); i++)
			if (NumberRecord(string(glob_result.gl_pathv[i])) < capacityPut)
				filename_consolidate.push_back(string(glob_result.gl_pathv[i])); 
	}
	else 
		filename_consolidate.clear(); 
        globfree(&glob_result);
        return filename_consolidate; 
}
int CPutGetBin::GetNumberFileForFetch() const
{
	vector <string> file_fetch = GetFileNameForFetch();
	return (int)(file_fetch.size()); 
}

int CPutGetBin::NumberRecord(string file_name) const
{
	ifstream iFile; 
	iFile.open(file_name.c_str(), ios::in|ios::binary);
        if (!iFile)
		return -1;
	
	// To determine size of each record
	CSampleIDWeight temp; 
	read(iFile, &temp); 

	iFile.seekg(0, ios::beg); 
	iFile.seekg(0, ios::end); 
	int lenFile = iFile.tellg(); 
	int number_record = lenFile/temp.GetSize_Data(); 
	iFile.close(); 
	return number_record; 
}


