/* Test: 
 * CEquiEnergy
 * CModel
 * CMixtureModel
 * CTransitionModel
 * CSimpleGaussianModel
 * CTransitionModel_SimpleGaussian
 * CUniformModel
 * CBoundedModel
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "CSampleIDWeight.h"
#include "CGaussianModel.h"
#include "CTransitionModel.h"
#include "CTransitionModel_Gaussian.h"
#include "CStorageHead.h"

using namespace std;

bool Configure_GaussianModel_File(CGaussianModel &, const string); 

int main()
{
	/*
 	Initialize the target distribution as a Gaussian model;
	Mean and Covariance are stored in files
  	*/
	string filename_base = "gaussian_model."; 
	CGaussianModel target; 
	if (!Configure_GaussianModel_File(target, filename_base))
	{
		cout << "Error in configuring gaussian mixture model.\n"; 
		exit (-1);
	} 

	/*
 	Initialize the random_number_generator which will be used for all distributions to draw samples.
 	*/
	const gsl_rng_type *T; 
	gsl_rng *r;
	gsl_rng_env_setup(); 
	T = gsl_rng_default; 
	r = gsl_rng_alloc(T); 
	gsl_rng_set(r, (unsigned)time(NULL)); 

	/*
 	storage 
 	*/
	int run_id = time(NULL); 
	CStorageHead storage(run_id, 1000, 1000, 1, string("")); 
	
	/*
 	Transition model: CTransitionModel_Gaussian with dim = 2 and initially set as an identity matrix
 	*/
	CTransitionModel **mProposal; 
	int nBlock = target.GetDataDimension()/2; 
	vector <int > blockSize(nBlock); 
	mProposal = new CTransitionModel *[nBlock]; 
	for (int iBlock=0; iBlock<nBlock; iBlock ++)
	{
		mProposal[iBlock] = new CTransitionModel_Gaussian(2);
		blockSize[iBlock] = 2; 
	}


	/* Tune mProposal */
	double targetAcc = 0.25; 
	int LPeriod = 100; 
	int NPeriod = 50; 
	CSampleIDWeight mode; 
	target.GetMode(mode); 
	int dim_cum_sum = 0; 
	for (int iBlock=0; iBlock<nBlock; iBlock ++)
	{
		mProposal[iBlock]->Tune(targetAcc, LPeriod, NPeriod, r, &target, mode, dim_cum_sum, 2); 
		dim_cum_sum += 2; 
	}	

	// simulation run
	CSampleIDWeight x_current, x_new; 
	x_current.SetDataDimension(target.GetDataDimension()); 
	x_current = mode; 
	x_new.SetDataDimension(target.GetDataDimension()); 
	vector <bool> if_new_sample(nBlock); 
	int iBlock; 
	for (int n=0; n<10000; n++)
	{
		target.drawMH(x_new, mProposal, if_new_sample, r, x_current, nBlock, blockSize, 0);
		iBlock =0; 
		while (iBlock < nBlock && !if_new_sample[iBlock])
			iBlock++; 
		if (iBlock < nBlock)
			x_current = x_new; 
		storage.DepositSample(0, x_current); 
	} 	
	storage.finalize(); 

	for (int iBlock=0; iBlock<nBlock; iBlock++)
		delete mProposal[iBlock]; 
	delete []mProposal; 
	gsl_rng_free(r); 
}

bool Configure_GaussianModel_File(CGaussianModel &model, const string filename_base)
{
	/*mean */
	string filename = filename_base + "mean"; 
	ifstream inputFile; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false;
	int nDim; 
	inputFile >> nDim; 
	double *mean = new double [nDim];
	for (int i=0; i<nDim; i++)
		inputFile >> mean[nDim]; 
	inputFile.close(); 
	
	/*covariance */
	filename = filename_base + "covariance"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nDim; 
	double* covariance = new double [nDim*nDim];
	for (int i=0; i<nDim; i++)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> covariance[i*nDim+j]; 
	}
	inputFile.close(); 

	model = CGaussianModel(nDim, mean, covariance); 
	delete [] mean; 
	delete [] covariance; 
	return true;
}
