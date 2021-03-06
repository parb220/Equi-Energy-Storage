#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_permutation.h>
#include "CStorageHead.h"

using namespace std;

CStorageHead::CStorageHead(int _run_id, int _get_marker, int _put_marker, int _number_bins, string file_location, int _node_index)
{
	cluster_node = _node_index;
	run_id = _run_id; 
	get_marker = _get_marker; 
	put_marker = _put_marker; 
	number_bins = _number_bins; 
	filename_base = file_location; 

	stringstream str; 
	str << run_id << "/" << run_id << ".binary/"; 
	CPutGetBin temp(0, 0, put_marker, get_marker, filename_base+str.str(), cluster_node); 
	bin = vector <CPutGetBin> (number_bins, temp); 
	for (int i=0; i<number_bins; i++)
		bin[i].SetBinID(i, cluster_node);
}

CStorageHead::~CStorageHead()
{
}


int CStorageHead::DepositSample(int _bin_id, const CSampleIDWeight &sample)
{
	return bin[_bin_id].DepositSample(sample); 
}

int CStorageHead::DepositSample(int _bin_id, const double *x, int x_d, int x_ID, double x_weight)
{
	return bin[_bin_id].DepositSample(x, x_d, x_ID, x_weight); 
} 

bool CStorageHead::DrawLeastWeightSample(int _bin_id, CSampleIDWeight &sample) const
{
	return bin[_bin_id].DrawLeastWeightSample(sample); 
}

bool CStorageHead::DrawLeastWeightSample(int _bin_id, double *x, int x_d, int &x_ID, double &x_weight) const
{
	return bin[_bin_id].DrawLeastWeightSample(x, x_d, x_ID, x_weight); 
}

bool CStorageHead::DrawMostWeightSample(int _bin_id, CSampleIDWeight &sample) const
{
	return bin[_bin_id].DrawMostWeightSample(sample); 
}

bool CStorageHead::DrawMostWeightSample(int _bin_id, double *x, int x_d, int &x_ID, double &x_weight) const
{
	return bin[_bin_id].DrawMostWeightSample(x, x_d, x_ID, x_weight); 
}

bool CStorageHead::DrawSample(int _bin_id, const gsl_rng *r, CSampleIDWeight &sample)
{
	return bin[_bin_id].DrawSample(r, sample); 	
}

bool CStorageHead::DrawSample(int _bin_id, double *x, int x_d, int &x_ID, double &x_weight, const gsl_rng *r)
{
	return bin[_bin_id].DrawSample(x, x_d, x_ID, x_weight, r); 
}

bool CStorageHead::makedir()
{
	stringstream str; 
	str.str(string()); 
	str << run_id; 
	string dir = filename_base + str.str(); 
	errno = 0; // clear error number
	int status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
	if (status !=0 && errno != EEXIST)
		return false;
 
	str.str(string()); 
	str << run_id << "/" << run_id << ".binary" ; 
	dir = filename_base + str.str();
	errno = 0; 
	status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 
	if (status == 0 || errno == EEXIST)
		return true; 
	else 
		return false; 
}

void CStorageHead::finalize(int start_bin, int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (int i=start_bin; i<=end_bin; i++)
			bin[i].finalize(); 
	}
	else
	{
		for (int i=0; i<number_bins; i++)
			bin[i].finalize(); 
	}
}


void CStorageHead::DisregardHistorySamples(int start_bin, int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (int i=start_bin; i<=end_bin; i++)
			bin[i].DisregardHistorySamples();
	}
	else 
	{
		for (int i=0; i<number_bins; i++)
			bin[i].DisregardHistorySamples(); 
	} 
}

void CStorageHead::ClearDepositDrawHistory(int start_bin, int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (int i=start_bin; i<=end_bin; i++)
			bin[i].ClearDepositDrawHistory(); 
	}
	else 
	{
		for (int i=0; i<number_bins; i++)
			bin[i].ClearDepositDrawHistory(); 
	}
}

void CStorageHead::consolidate(int start_bin, int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (int i=start_bin; i<=end_bin; i++)
			bin[i].consolidate(); 
	}
	else 
	{
		for (int i=0; i<number_bins; i++)
			bin[i].consolidate(); 
	}
}

void CStorageHead::restore(int start_bin, int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
        {
                for (int i=start_bin; i<=end_bin; i++)
                        bin[i].restore();
        }
	else 
	{
		for (int i=0; i<number_bins; i++)
			bin[i].restore(); 
	}
}

bool CStorageHead::empty(int start_bin, int end_bin) const
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
	{
		for (int i=end_bin; i>=start_bin; i--)
			if (!bin[i].empty())
				return false; 
		return true; 
	}
	else 
	{
		for (int i=number_bins-1; i>=0; i--)
			if (!bin[i].empty())
				return false; 
		return true; 
	}
}

void CStorageHead::RestoreForFetch(int start_bin, int end_bin)
{
	if (start_bin >=0 && start_bin < number_bins && end_bin >= 0 && end_bin < number_bins && start_bin <= end_bin)
        {
                for (int i=start_bin; i<=end_bin; i++)
                        bin[i].RestoreForFetch();
        }
	else 
	{
		for (int i=0; i<number_bins; i++)
			bin[i].RestoreForFetch(); 
	}

}


void CParameterPackage::TraceStorageHead(const CStorageHead &storage)
{
       	// number_samples_generated_by_far.resize(storage.number_bins);
	number_files_fetch.resize(storage.number_bins); 
        for (int i=0; i<storage.number_bins; i++)
	{
                // number_samples_generated_by_far[i] = storage.bin[i].GetNumberFileForFetch();
		number_files_fetch[i] = storage.bin[i].GetNumberFileForFetch(); 
	}
}

bool CParameterPackage::LoadCurrentStateFromStorage(CStorageHead &storage, const gsl_rng *r, int level)
{
	x_current.resize(number_energy_level);
        for (int i=0; i<number_energy_level; i++)
                x_current[i].SetDataDimension(data_dimension);
	
	int startLevel, endLevel; 
	if (level >=0 && level<number_energy_level)
	{
		startLevel = level; 
		endLevel = level +1; 
	}	
	else 
	{
		startLevel = 0; 
		endLevel = number_energy_level; 
	}

	int bin_id, startBin, endBin;  
	bool if_success; 
	gsl_permutation *p; 

	for (int i=startLevel; i<endLevel; i++)
	{
		if_success = false; 
		startBin = (i+1)*number_energy_level; 
		endBin = number_bins-1; 
		
		p = gsl_permutation_alloc(endBin-startBin+1); 
		gsl_permutation_init(p); 
		gsl_ran_shuffle(r, p->data, endBin-startBin+1, sizeof(int)); 
		int index = 0; 
		while (!if_success && index < (endBin-startBin+1))
		{
			bin_id = startBin + gsl_permutation_get(p,index); 
			if_success = storage.DrawSample(bin_id, r, x_current[i]); 
			index ++; 
		}
		gsl_permutation_free(p); 
		if (!if_success)
		{
			startBin = (i+1)*number_energy_level-1; 
			endBin = 0; 
			p = gsl_permutation_alloc(startBin-endBin+1); 
			gsl_permutation_init(p); 
			gsl_ran_shuffle(r, p->data, startBin-endBin+1, sizeof(int)); 
			int index = 0; 
			while(!if_success && index < (startBin-endBin+1))
			{
				bin_id = startBin - gsl_permutation_get(p, index);  
				if_success = storage.DrawSample(bin_id, r, x_current[i]); 
				index ++; 
			}
			gsl_permutation_free(p); 
		}		
		if (!if_success)
			return false; 
	}
	return true; 
}
