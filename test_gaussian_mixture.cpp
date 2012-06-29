/* Test: 
 * CEquiEnergy
 * CModel
 * CMixtureModel
 * CTransitionModel
 * CSimpleGaussianModel
 * CTransitionModel_SimpleGaussian
 * CUniformModel
 * CBoundedModel
 * CEES_Head
 * CEES_Node
 */

#include <fstream>
#include <iostream>
#include <cstring>
#include <ctime>
#include <cmath>
#include <gsl/gsl_rng.h>
#include "constant.h"
#include "equi_energy_setup_constant.h"
#include "CMixtureModel.h"
#include "CSimpleGaussianModel.h"
#include "CEquiEnergy.h"
#include "CUniformModel.h"
#include "CBoundedModel.h"
#include "CTransitionModel.h"
#include "CTransitionModel_SimpleGaussian.h"
#include "CEES_Node.h"
#include "CSampleIDWeight.h"
#include "CStorageHead.h"

using namespace std;

bool Configure_GaussianMixtureModel_File(CMixtureModel &, const string); 

int main()
{
	/*
 	Initialize the target distribution as a Gaussian mixture model;
	Mean Sigma and Weight are stored in files
  	*/
	string filename_base = "../equi_energy_generic/gaussian_mixture_model."; 
	CMixtureModel target; 
	if (!Configure_GaussianMixtureModel_File(target, filename_base))
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

	/* Initialize Storage  */
	int run_id = time(NULL); // by default, use current time as run_id; 
	int get_marker = 10000;	// Keep get_marker samples in memory for draw 
	int put_marker = 10000;	// Dump every put_marker samples
	int number_bins = NUMBER_ENERGY_LEVEL * NUMBER_ENERGY_LEVEL;   
	CSampleIDWeight::SetDataDimension(DATA_DIMENSION);	// Data dimension for storage (put and get) 
	CStorageHead storage(run_id, get_marker, put_marker, number_bins, string("/home/f1hxw01/equal_energy_hw/equi_energy_storage_data/")); 
	storage.makedir();
		
	/*
 	Initialize an equi_energy object for sampling
 	*/
	CEES_Node::SetEnergyLevelNumber(NUMBER_ENERGY_LEVEL); 		// Number of energy levels; 
	CEES_Node::SetEquiEnergyJumpProb(PEE);				// Probability for equal energy jump
	CEES_Node::SetPeriodBuildInitialRing(BUILD_INITIAL_ENERGY_SET_PERIOD);	// Period to build initial energy ring
	CEES_Node::SetDataDimension(DATA_DIMENSION); 		// Data dimension for simulation
	CEES_Node::SetDepositFreq(DEPOSIT_FREQUENCY); 		// Frequency of deposit
	CEES_Node::ultimate_target = &target;	

	/*
 	Set energy levels according to the geometric progression given H0 and H[K-1]
	Alternatively, could use SetEnergyLevels(double*, int) 
 	*/
	if (!CEES_Node::SetEnergyLevels_GeometricProgression(H0, HK_1))
	{
		cout << "Error in setting energy levels." << endl; 
		exit(-1); 
	}

	/*
 	Set temperatures for all levels, either according to the energy levels so that (H[i+1]-H[i])/T[i] is a constant, or use SetTemperatures(double*, int)
 	*/
	if (!CEES_Node::SetTemperatures_EnergyLevels(T0, TK_1, C) )
	{
		cout << "Error in setting temperature levels." << endl; 
		exit(-1);
	}

	/*  Generate K CEES_Node objects */
	CEES_Node *simulator_node = new CEES_Node[CEES_Node::GetEnergyLevelNumber()]; 
	
	double *sigma = new double[CEES_Node::GetDataDimension()];
	/* Uniform distribution in [0, 1]^2 used for initialization. */
	double *lB = new double [CEES_Node::GetDataDimension()]; 
	double *uB = new double [CEES_Node::GetDataDimension()]; 
	for (int i=0; i<CEES_Node::GetDataDimension(); i++)
	{
		lB[i] = 0.0; 
		uB[i] = 1.0; 
	}
	CModel *initial_model = new CUniformModel(CEES_Node::GetDataDimension(), lB, uB); 
	for (int i=0; i<CEES_Node::GetEnergyLevelNumber(); i++)
	{
		simulator_node[i].SetID_LocalTarget(i); 
		if (i < CEES_Node::GetEnergyLevelNumber() -1)
		{
			simulator_node[i].SetBurnInPeriod(0);
			simulator_node[i].SetHigherNodePointer(simulator_node+i+1);
		}
		else 
		{
			simulator_node[i].SetHigherNodePointer(NULL);
			simulator_node[i].SetBurnInPeriod(BURN_IN_PERIOD);  
		}
		/* Transition_SimpleGaussian with sigma=0.25sqrt(T) used for each energy level as the proposal model */	
                for (int j=0; j<CEES_Node::GetDataDimension(); j++)
                        sigma[j] = 0.25 * sqrt(simulator_node[i].GetTemperature());
                simulator_node[i].SetProposal(new CTransitionModel_SimpleGaussian(CEES_Node::GetDataDimension(), sigma));
		/*   Initialize       */
		simulator_node[i].Initialize(initial_model, r); 
	}
	delete initial_model;
	delete [] lB; 
	delete [] uB; 
	delete [] sigma; 

	/*
 	Burn-in and simulation
 	*/
	for (int n=0; n<(CEES_Node::GetEnergyLevelNumber()-1)*CEES_Node::GetPeriodBuildInitialRing()+SIMULATION_LENGTH; n++)
	{
		simulator_node[CEES_Node::GetEnergyLevelNumber()-1].draw(r, storage); 
		for (int i=CEES_Node::GetEnergyLevelNumber()-2; i>=0; i--)
		{
			if (simulator_node[i+1].EnergyRingBuildDone())
				simulator_node[i].draw(r, storage); 
		}
	}

	storage.finalize(); 		// save to hard-disk of those unsaved data
	string file = storage.GetSummaryFileName(); 
	ofstream oFile; 
	oFile.open(file.c_str());
	if (!oFile)
	{
		cout << "Error in writing the summary file.\n"; 
		exit(-1); 
	}
	summary(oFile, storage); 
	summary(oFile, simulator_node[CEES_Node::GetEnergyLevelNumber()-1]);   
	oFile.close(); 
	
	gsl_rng_free(r); 
}

bool Configure_GaussianMixtureModel_File(CMixtureModel &mixture_model, const string filename_base)
{
	/*weight */
	string filename = filename_base + "weight";
	int nComponent, nDim; 			// Number of components, dimension of variables
	bool equalComponent, equalDim;		// Whether parameters for different components are thes same; and whether parameters for different dimension are the same; 
	ifstream inputFile; 
	inputFile.open(filename.data());
	if (!inputFile)
		return false;
	inputFile >> nComponent; 
	double *weight = new double[nComponent]; 
	inputFile >> equalComponent; 
	if (!equalComponent)
	{
		for (int i=0; i<nComponent; i++)
			inputFile >> weight[i]; 
	} 
	else 
	{
		inputFile >> weight[0]; 
		for (int i=1; i<nComponent; i++)
			weight[i] = weight[0];
	}
	mixture_model.SetModelNumber(nComponent); 
	mixture_model.SetWeightParameter(weight, nComponent); 
	delete [] weight; 
	inputFile.close();

	/*sigma */
	filename = filename_base + "sigma"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** sigma = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		sigma[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> sigma[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> sigma[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				sigma[i][j] = sigma[i][0]; 
		} 
	}
	else
	{
		inputFile >> sigma[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				sigma[i][j] = sigma[0][0];
		}
	}

	inputFile.close(); 

	/*mean */
	filename = filename_base + "mean"; 
	inputFile.open(filename.data()); 
	if (!inputFile)
		return false; 
	inputFile >> nComponent >> nDim; 
	double** mean = new double* [nComponent];
	for (int i=0; i<nComponent; i++)
		mean[i] = new double[nDim]; 
	inputFile >> equalComponent >> equalDim; 
	if (!equalComponent && !equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				inputFile >> mean[i][j]; 
		}
	}
	else if (equalComponent && !equalDim)
	{
		for (int j=0; j<nDim; j++)
			inputFile >> mean[0][j]; 
		for (int i=1; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][j];
		}
	}
	else if (!equalComponent && equalDim)
	{
		for (int i=0; i<nComponent; i++)
		{
			inputFile >> sigma[i][0];
			for (int j=1; j<nDim; j++)
				mean[i][j] = mean[i][0]; 
		} 
	}
	else
	{
		inputFile >> mean[0][0]; 
		for (int i=0; i<nComponent; i++)
		{
			for (int j=0; j<nDim; j++)
				mean[i][j] = mean[0][0];
		}
	}

	inputFile.close(); 

	mixture_model.SetDataDimension(nDim); 
	for (int i=0; i<nComponent; i++)
	{
		mixture_model.Initialize(i, new CSimpleGaussianModel(nDim, mean[i], sigma[i])); 
		delete []mean[i]; 
		delete []sigma[i]; 
	}
	mixture_model.CalculateSetParameterNumber();
	return true;
}
