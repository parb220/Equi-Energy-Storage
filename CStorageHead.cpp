#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include "CStorageHead.h"

using namespace std;

CStorageHead::CStorageHead(int _run_id, int _get_marker, int _put_marker, int _number_bins, string file_location)
{
	run_id = _run_id; 
	get_marker = _get_marker; 
	put_marker = _put_marker; 
	number_bins = _number_bins; 
	filename_base = file_location; 

	stringstream str; 
	str << run_id << ".binary/"; 
	CPutGetBin temp(0, 0, put_marker, get_marker, filename_base+str.str()); 
	bin = vector <CPutGetBin> (number_bins, temp); 
	for (int i=0; i<number_bins; i++)
		bin[i].SetBinID(i);
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

bool CStorageHead::DrawSample(int _bin_id, const gsl_rng *r, CSampleIDWeight &sample)
{
	return bin[_bin_id].DrawSample(r, sample); 	
}

bool CStorageHead::DrawSample(int _bin_id, double *x, int x_d, int &x_ID, double &x_weight, const gsl_rng *r)
{
	return bin[_bin_id].DrawSample(x, x_d, x_ID, x_weight, r); 
}

bool CStorageHead::empty(int _bin_id)
{
	if (bin[_bin_id].GetNumberSamplesGeneratedByFar() <= 0)
		return true; 
	else 
		return false; 
}

string CStorageHead::GetSummaryFileName() const
{
	stringstream convert; 
	convert << run_id; 
	return filename_base + convert.str() + ".summary"; 
}

ofstream & summary(ofstream &of, const CStorageHead &head)
{
	of << "Get Marker:\t" << head.get_marker << endl; 
	of << "Put Marker:\t" << head.put_marker << endl; 
	of << "Number of Bins: \t" << head.number_bins << endl; 
	for (int i=0; i<head.number_bins; i++)
	{
		of << head.bin[i].GetBinID() << ":\t" << head.bin[i].GetNumberSamplesGeneratedByFar() << "\t" << head.bin[i].GetNumberDataFile(true) << endl; 
	}	
	return of;	
}

bool CStorageHead::makedir()
{
	stringstream str; 
	str << run_id << ".binary" ; 
	string dir = filename_base + str.str(); 
	int status = mkdir(dir.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH); 
	if (status == 0 || status == EEXIST)
		return true; 
	else 
		return false; 
}

void CStorageHead::finalize()
{
	for (int i=0; i<number_bins; i++)
		bin[i].finalize(); 
}

vector <CSampleIDWeight> CStorageHead::RetrieveSamplesSequentially(bool if_clear_old_bin, int bin_id)
{
	vector <CSampleIDWeight> samples = bin[bin_id].RetrieveSamplesSequentially(if_clear_old_bin); 
	return samples; 
}

void CStorageHead::DisregardHistorySamples(int bin_id)
{
	bin[bin_id].DisregardHistorySamples(); 
}

void CStorageHead::ClearDepositDrawHistory(int bin_id)
{
	bin[bin_id].ClearDepositDrawHistory(); 
}

void CStorageHead::CreateTemporaryBin()
{
	if ((int)(bin.size()) < 2*number_bins)
	{
	 	int old_size = (int)(bin.size());
                bin.resize(2*number_bins);
                for (int i=old_size; i<(int)(bin.size()); i++)
                {
                        bin[i].SetBinID(i);
                        bin[i].SetCapacity_Put(bin[i-number_bins].GetCapacity_Put());
                        bin[i].SetCapacity_Get(bin[i-number_bins].GetCapacity_Get());
                        bin[i].SetFileNamePrefix(bin[i-number_bins].GetFileNamePrefix());
                        bin[i].ClearDepositDrawHistory();
                }
	}
}


int CStorageHead::DepositSample(bool if_new_bin, int bin_id, const double *x, int dX, int _id, double _weight)
{
	/* First number_bins: old
 	* Last number_bins: new */ 
	
	if (if_new_bin)
		return bin[bin_id+number_bins].DepositSample(x, dX, _id, _weight); 
	else
		return bin[bin_id].DepositSample(x, dX, _id, _weight); 
} 

int CStorageHead::DepositSample(bool if_new_bin, int bin_id, const CSampleIDWeight &sample)
{
	/* First number_bins: old 
 	*  others: new  */
	if (if_new_bin)	
		return bin[bin_id+number_bins].DepositSample(sample); 
	else 
		return bin[bin_id].DepositSample(sample); 
}

void CStorageHead::Consolidate(int bin_id)
{
	int old_bin_id = bin_id; 
	int new_bin_id = bin_id + number_bins;
	if ((int)(bin.size()) > new_bin_id && bin[new_bin_id].GetNumberSamplesGeneratedByFar())
	{
		bin[new_bin_id].ChangeFileName(old_bin_id); 
		bin[old_bin_id] = bin[new_bin_id]; 
		bin[old_bin_id].SetBinID(old_bin_id); 
	}
	else 
		bin[old_bin_id].ClearDepositDrawHistory(); 
}

void CStorageHead::ClearTemporaryBin()
{
	if ((int)(bin.size()) > number_bins)
		bin.erase(bin.begin()+number_bins, bin.end()); 
}
