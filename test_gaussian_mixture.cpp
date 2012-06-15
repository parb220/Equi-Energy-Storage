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
#include "sdsm.h"

using namespace std;

bool Configure_GaussianMixtureModel_File(CMixtureModel &, const string); 
void handler(void *tag)
{
	cout << ((char*)tag)<< endl;
}

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

	/*
 	Initialize an equi_energy object for sampling
 	*/
	CEES_Node::SetEnergyLevelNumber(NUMBER_ENERGY_LEVEL); 		// Number of energy levels; 
	CEES_Node::SetEquiEnergyJumpProb(PEE);				// Probability for equal energy jump
	CEES_Node::SetBurnInPeriod(BURN_IN_PERIOD);			// Burn-in period
	CEES_Node::SetPeriodBuildInitialRing(BUILD_INITIAL_ENERGY_SET_PERIOD);	// Period to build initial energy ring
	CEES_Node::SetDataDimension(DATA_DIMENSION); 		// Data dimension 
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
		simulator_node[i].SetID(i);
		if (i < CEES_Node::GetEnergyLevelNumber() -1)
			simulator_node[i].SetHigherNodePointer(simulator_node+i+1);
		else 
			simulator_node[i].SetHigherNodePointer(NULL); 
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
 	SDSM initialization and set up
 	*/
	int run_id = time(NULL);
	bool if_restart_old_run = false; 
	int get_marker = 1000; 
	int put_marker = 1000; 
	int number_bins = CEES_Node::GetEnergyLevelNumber()*CEES_Node::GetEnergyLevelNumber();
	bool if_uniformly_weighted_items = true; 
	bool if_database_secondary_memory_available = true; 
	char db_serv_adrs[] = "localhost";
	bool if_cluster = false; 
	int cluster_task_id = 0; 	// only matters when running on clusters
	char tags[] = "equi-energy, sdsm, test for gaussian mixture";
	sdsm_init(run_id, if_restart_old_run, get_marker, put_marker, number_bins, if_uniformly_weighted_items, if_database_secondary_memory_available, db_serv_adrs, if_cluster, cluster_task_id, tags); 

	/* testing sdsm */
	string finish_remark("Finishing.");  
	sdsm_register_handler(handler, (void *)(finish_remark.data()));
	/* */
	
	/*
 	Burn-in and simulation
 	*/
	for (int n=0; n<(CEES_Node::GetEnergyLevelNumber()-1)*(CEES_Node::GetBurnInPeriod()+CEES_Node::GetPeriodBuildInitialRing())+CEES_Node::GetBurnInPeriod()+SIMULATION_LENGTH; n++)
	{
		for (int i=CEES_Node::GetEnergyLevelNumber()-1; i>=0; i--)
		{
			if (n >= (CEES_Node::GetEnergyLevelNumber()-1-i)*(CEES_Node::GetBurnInPeriod()+CEES_Node::GetPeriodBuildInitialRing()) )
				simulator_node[i].draw(r);
		}
	}

		
	sdsm_finalize(); 
	/*
	Output samples into files
 		
	filename_base = "gaussian_mixture_sample."; 
	char *char_buffer = new char[100]; 
	string filename; 
	for (int i=0; i<equi_energy_simulator.GetNumberEnergyLevels(); i++)
	{
		memset(char_buffer, 0, sizeof(char_buffer)); 
		sprintf(char_buffer, "%d", i); 
		filename = filename_base + char_buffer; 
		if (equi_energy_simulator.Output_Samples_EnergyLevel_File(i, filename) < 0)
		{
			cout << "Error in outputting the " << i << "-th energy level samples.\n";
			exit(-1);
		}
	}
	delete [] char_buffer; 
	*/
	/*
 	Release random number generator
 	*/
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
