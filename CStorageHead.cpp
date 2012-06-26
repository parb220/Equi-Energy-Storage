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
	str << run_id << "/"; 
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

CSampleIDWeight CStorageHead::DrawSample(int _bin_id, const gsl_rng *r)
{
	return bin[_bin_id].DrawSample(r); 	
}

void CStorageHead::DrawSample(int _bin_id, double *x, int x_d, int &x_ID, double &x_weight, const gsl_rng *r)
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
	str << run_id; 
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
