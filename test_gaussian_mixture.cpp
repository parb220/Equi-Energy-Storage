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
#include <cstdlib>
#include <unistd.h>
#include <gsl/gsl_rng.h>
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
bool TuneEnergyLevels_UpdateStorage(CEES_Node *, CStorageHead &);

void usage(int arc, char **argv)
{
	cout << "usage: " << argv[0] << " "; 
	cout << "-d <dimension> \n"; 
	cout << "-f <files of the target model> \n"; 
	cout << "-p <probability of equi-energy jump> \n"; 
	cout << "-h <energy bound of the highest energy level>\n"; 
	cout << "-l <simulation length>\n"; 
	cout << "? this message.\n";
}


int main(int argc, char **argv)
{
	// default setting, which is obtained from .h file
	string filename_base = "../equi_energy_generic/gaussian_mixture_model."; 
	int data_dimension = DATA_DIMENSION; 
	double pee = PEE;  
	double h_k_1 = HK_1; 
	int simulation_length =SIMULATION_LENGTH; 

	// parse command line
	int opt; 
	while ( (opt = getopt(argc, argv, "d:f:p:h:l:?")) != -1)
	{
		switch (opt) 
		{
			case 'f':
				filename_base = string(optarg); break; 
			case 'd':
				data_dimension = atoi(optarg); break; 	
			case 'p':
				pee = atof(optarg); break; 			
			case 'h':
				h_k_1 = atof(optarg); break; 
			case 'l':
				simulation_length = atoi(optarg); break; 
			case '?': 
			{
				usage(argc, argv); 
				exit(-1); 
			}	
			default:
			{ 
				usage(argc, argv); 
				return (-1); 
			}
		}
	}
	
	/*
 	Initialize the target distribution as a Gaussian mixture model;
	Mean Sigma and Weight are stored in files
  	*/
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
	CSampleIDWeight::SetDataDimension(data_dimension);	// Data dimension for storage (put and get) 
	CStorageHead storage(run_id, get_marker, put_marker, number_bins, string("/home/f1hxw01/equal_energy_hw/equi_energy_storage_data/")); 
	storage.makedir();
		
	/*
 	Initialize an equi_energy object for sampling
 	*/
	CEES_Node::SetEnergyLevelNumber(NUMBER_ENERGY_LEVEL); 		// Number of energy levels; 
	CEES_Node::SetEquiEnergyJumpProb(pee);				// Probability for equal energy jump
	CEES_Node::SetDataDimension(data_dimension); 		// Data dimension for simulation
	CEES_Node::ultimate_target = &target;			// Ultimate target distribution

	/*
 	Set energy levels according to the geometric progression given H0 and H[K-1]
	Alternatively, could use SetEnergyLevels(double*, int) 
 	*/
	CEES_Node::SetEnergyLevels_GeometricProgression(H0, h_k_1); 
	CEES_Node::InitializeMinMaxEnergy(ENERGY_TRACKING_NUMBER); // Initialize min_energy and if_tune_energy_level for tune_energy_level in the future; 
	
	/*
 	Set temperatures for all levels, either according to the energy levels, or use SetTemperatures(double*, int)
 	*/
	CEES_Node::SetTemperatures_EnergyLevels(T0, C, true); // (H[i+1]-H[i])/(T[i+1]-T[i]) is a constant

	/* Set target acceptance rate for each level's MH, (H(i+1)-H(i))/(targetAcc(i+1)-targetAcc(i)) = c */
	CEES_Node::SetTargetAcceptanceRate(MH_TARGET_ACC);
	// MH_TARGET_ACC:	target acceptance rate for the lowest energy level
	// 1.0:			target acceptance rate for the highest energy level 

	/* If MH block will be used */
	if (MH_BLOCK)	
		CEES_Node::SetBlockSize(NULL, CEES_Node::GetDataDimension()); 
	else 
		CEES_Node::SetBlockSize(NULL); 

	/*  Generate K CEES_Node objects */
	CEES_Node *simulator_node = new CEES_Node[CEES_Node::GetEnergyLevelNumber()]; 
	
	// Local target, pointer to the higher level node, MH Proposal distribution 
	double *sigma; 
 	for (int i=0; i<CEES_Node::GetEnergyLevelNumber(); i++)
	{
		simulator_node[i].SetID_LocalTarget(i); 	// Also set the local target distribution here
		if (i < CEES_Node::GetEnergyLevelNumber() -1)
			simulator_node[i].SetHigherNodePointer(simulator_node+i+1);
		else 
			simulator_node[i].SetHigherNodePointer(NULL);
		for (int iBlock = 0; iBlock <CEES_Node::GetNumberBlocks(); iBlock++)
		{
			sigma = new double[CEES_Node::GetBlockSize(iBlock)]; 
			for (int j=0; j<CEES_Node::GetBlockSize(iBlock); j++)
				sigma[j] = INITIAL_SIGMA * exp(log(simulator_node[i].GetTemperature()));
			simulator_node[i].SetProposal(new CTransitionModel_SimpleGaussian(CEES_Node::GetBlockSize(iBlock), sigma), iBlock); 
			delete [] sigma; 
		}
	}

	/*
 	Initialize --> BurnIn --> loop for several times (Tune MH StepSize --> Build Initial Energy Ring (Simulate for N samples) --> TuneEnergyLevel) --> simulate
 	*/
	double *lB = new double [CEES_Node::GetDataDimension()]; 
	double *uB = new double [CEES_Node::GetDataDimension()]; 
	for (int i=0; i<CEES_Node::GetDataDimension(); i++)
	{
		lB[i] = 0.0; 
		uB[i] = 1.0; 
	}
	CModel *initial_model = new CUniformModel(CEES_Node::GetDataDimension(), lB, uB); 
	for (int i=CEES_Node::GetEnergyLevelNumber()-1; i>=0; i--)
	{
		if (i == CEES_Node::GetEnergyLevelNumber()-1)
		/* Use a uniform distribution model for initialization. */
			simulator_node[i].Initialize(initial_model, r); 
		else if (!simulator_node[i].Initialize(storage, r)) 
			simulator_node[i].Initialize(initial_model, r); 
		cout << "Node " << i << " Burn in for " << BURN_IN_PERIOD << endl; 
		simulator_node[i].BurnIn(r, storage, BURN_IN_PERIOD, MULTIPLE_TRY_MH);
		cout << "Node " << i << " MH StepSize Tuning for " << MH_STEPSIZE_TUNING_MAX_TIME << " beginning for " << MH_TRACKING_LENGTH << endl;  
		simulator_node[i].MH_StepSize_Estimation(MH_TRACKING_LENGTH, MH_STEPSIZE_TUNING_MAX_TIME, r, MULTIPLE_TRY_MH);
		cout << "Node " << i << " simulate for " << ENERGY_LEVEL_TRACKING_WINDOW_LENGTH << endl; 
		simulator_node[i].Simulate(r, storage, ENERGY_LEVEL_TRACKING_WINDOW_LENGTH, DEPOSIT_FREQUENCY, MULTIPLE_TRY_MH);  
	}
	
	delete initial_model;
	delete [] lB; 
	delete [] uB; 

	// Tuning; 	
	int nEnergyLevelTuning = 0; 
	while (nEnergyLevelTuning < ENERGY_LEVEL_TUNING_MAX_TIME)
	{
		cout << "Energy level tuning " << nEnergyLevelTuning << endl; 
		TuneEnergyLevels_UpdateStorage(simulator_node, storage);
		for (int i=CEES_Node::GetEnergyLevelNumber()-1; i>=0; i--)
		{
               		simulator_node[i].MH_StepSize_Estimation(MH_TRACKING_LENGTH, MH_STEPSIZE_TUNING_MAX_TIME, r, MULTIPLE_TRY_MH);
                	simulator_node[i].Simulate(r, storage, ENERGY_LEVEL_TRACKING_WINDOW_LENGTH, DEPOSIT_FREQUENCY, MULTIPLE_TRY_MH);
		}
		nEnergyLevelTuning ++; 
	}
	
	cout << "Simulate for " << simulation_length<< endl; 
	// Runing through
	for (int i=CEES_Node::GetEnergyLevelNumber()-1; i>=0; i--)
		simulator_node[i].Simulate(r, storage, simulation_length, DEPOSIT_FREQUENCY, MULTIPLE_TRY_MH);

	storage.finalize(); 		// save to hard-disk of those unsaved data
	string file = storage.GetSummaryFileName(); 
	ofstream oFile; 
	oFile.open(file.c_str());
	if (!oFile)
	{
		cout << "Error in writing the summary file.\n"; 
		exit(-1); 
	}
	oFile << "Burn In:\t" << BURN_IN_PERIOD << endl; 
	oFile << "Tune Energy Level Window Length:\t" << ENERGY_LEVEL_TRACKING_WINDOW_LENGTH << endl; 
	oFile << "Tune Energy Level Number:\t" << ENERGY_LEVEL_TUNING_MAX_TIME << endl;
	oFile << "Deposit Frequency:\t"	<< DEPOSIT_FREQUENCY << endl; 
	oFile << "MH Target Probability:\t" << MH_TARGET_ACC << endl; 
	oFile << "MH Initial Window Length:\t" << MH_TRACKING_LENGTH << endl; 
	oFile << "MH Window Number:\t" << MH_STEPSIZE_TUNING_MAX_TIME << endl;   
	summary(oFile, storage); 
	summary(oFile, simulator_node);   
	for (int i=0; i<CEES_Node::GetEnergyLevelNumber(); i++)
		summary(oFile, simulator_node[i], i); 

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
