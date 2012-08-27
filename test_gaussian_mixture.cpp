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
#include <sstream>
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
#include "CParameterPackage.h"

using namespace std;

bool Configure_GaussianMixtureModel_File(CMixtureModel &, const string); 
bool TuneEnergyLevels_UpdateStorage(CEES_Node *, CStorageHead &, CParameterPackage &);

void usage(int arc, char **argv)
{
	cout << "usage: " << argv[0] << endl; 
	cout << "-i <id>: id of simulation run\n"; 
	cout << "-y: to continue a previous simulation run (when -i is provided)\n"; 
	cout << "-d <dimension>: dimension of samples\n"; 
	cout << "-f <path>: path of the file containing the target distribution specifications\n"; 
	cout << "-p <probability>: probability of equi-energy jump\n"; 
	cout << "-h <energy>: energy bound of the highest energy level\n"; 
	cout << "-l <length>: length of simulation\n"; 
	cout << "-c <coefficient>: C factor to determine temperature bounds according to energy bounds\n";
	cout << "-t <number>: number of times max and min energy bounds are tracked and tuned\n"; 	
	cout << "-b <path>: directory to store samples\n"; 
	cout << "? this message.\n";
}


int main(int argc, char **argv)
{
 	// Initialize the random_number_generator
	const gsl_rng_type *T; 
	gsl_rng *r;
	gsl_rng_env_setup(); 
	T = gsl_rng_default; 
	r = gsl_rng_alloc(T); 
	gsl_rng_set(r, (unsigned)time(NULL)); 
	
	// default setting, which is obtained from .h file
	int _run_id = time(NULL); // by default, use current time as run_id; 
	bool if_continue = false; 
	string target_filename_base = "../equi_energy_generic/gaussian_mixture_model."; 
	string storage_filename_base = string("/home/f1hxw01/equal_energy_hw/equi_energy_storage_data/"); 
	int _data_dimension = DATA_DIMENSION; 
	double _pee = PEE;  
	double _h_k_1 = HK_1; 
	int _simulation_length =SIMULATION_LENGTH; 
	double _c_factor = C; 
	double _mh_target_acc = MH_TARGET_ACC; 
	double _energy_level_tuning_max_time =  ENERGY_LEVEL_TUNING_MAX_TIME; 

	// parse command line
	int opt; 
	while ( (opt = getopt(argc, argv, "i:yd:f:p:h:l:t:c:b:?")) != -1)
	{
		switch (opt) 
		{
			case 'i':
				_run_id = atoi(optarg); break; 
			case 'y':
				if_continue = true; break; 
			case 'f':
				target_filename_base = string(optarg); break; 
			case 'b':
				storage_filename_base = string(optarg); break; 
			case 'd':
				_data_dimension = atoi(optarg); break; 	
			case 'p':
				_pee = atof(optarg); break; 			
			case 'h':
				_h_k_1 = atof(optarg); break; 
			case 'l':
				_simulation_length = atoi(optarg); break; 
			case 'c':
				_c_factor = atof(optarg); 
			case 't': 
				_energy_level_tuning_max_time = atoi(optarg); break; 
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
	
	// Initialize parameters	
	CParameterPackage parameter; 
	stringstream convert; 
	string file_name; 
	if (if_continue)
	{
		convert.str(std::string()); 
		convert << _run_id << ".parameter"; 
		file_name = storage_filename_base + convert.str(); 
		parameter.LoadParameterFromFile(file_name); 
	}
	else
	{
		parameter.run_id = _run_id; 
		parameter.get_marker = 10000; 
		parameter.put_marker = 10000; 
		parameter.number_energy_level = NUMBER_ENERGY_LEVEL; 
		parameter.data_dimension = _data_dimension; 
		parameter.number_bins = parameter.number_energy_level * parameter.number_energy_level; 
		parameter.pee = _pee; 
		parameter.h0 = H0; 
		parameter.hk_1 = _h_k_1; 
		parameter.energy_tracking_number = ENERGY_TRACKING_NUMBER;
		parameter.t0 = T0; 
		parameter.c_factor = _c_factor;
		parameter.mh_target_acc = _mh_target_acc; 
		parameter.initial_sigma = INITIAL_SIGMA; 
		parameter.uniform_lb = 0.0; 
		parameter.uniform_ub = 1.0; 
		parameter.burn_in_period = BURN_IN_PERIOD; 
		parameter.multiple_try_mh = MULTIPLE_TRY_MH;
		parameter.mh_tracking_length = MH_TRACKING_LENGTH;
		parameter.mh_stepsize_tuning_max_time = MH_STEPSIZE_TUNING_MAX_TIME; 
		parameter.energy_level_tracking_window_length = ENERGY_LEVEL_TRACKING_WINDOW_LENGTH; 
		parameter.energy_level_tuning_max_time = _energy_level_tuning_max_time; 
		parameter.deposit_frequency = DEPOSIT_FREQUENCY; 

		if (MH_BLOCK)
			parameter.number_block = parameter.data_dimension; 
		else 
			parameter.number_block = 1; 
		parameter.SetBlock(); 
		parameter.SetEnergyBound(); 
		parameter.SetTemperature(); 
		parameter.SetMHProposalScale(); 
		parameter.SetCurrentState(r); 
	}
	parameter.simulation_length = _simulation_length; 
	
	/*
 	Initialize the target distribution as a Gaussian mixture model;
	Mean Sigma and Weight are stored in files
  	*/
	CMixtureModel target; 
	if (!Configure_GaussianMixtureModel_File(target, target_filename_base))
	{
		cout << "Error in configuring gaussian mixture model.\n"; 
		exit (-1);
	} 



	/* Initialize Storage  */
	CSampleIDWeight::SetDataDimension(parameter.data_dimension);	// Data dimension for storage (put and get)  
	CStorageHead storage(parameter.run_id, parameter.get_marker, parameter.put_marker, parameter.number_bins, storage_filename_base); 
	if (if_continue)
		storage.restore(parameter); 
	else 
		storage.makedir();

			
	/*
 	Initialize an equi_energy object for sampling
 	*/
	CEES_Node::SetEnergyLevelNumber(parameter.number_energy_level);	// Number of energy levels; 
	CEES_Node::SetEquiEnergyJumpProb(parameter.pee);	// Probability for equal energy jump
	CEES_Node::SetDataDimension(parameter.data_dimension); 	// Data dimension for simulation
	CEES_Node::ultimate_target = &target;			// Ultimate target distribution

	/*
 	Set energy levels according to the geometric progression given H0 and H[K-1]
	Alternatively, could use SetEnergyLevels(double*, int) 
 	*/
	CEES_Node::InitializeMinMaxEnergy(parameter.energy_tracking_number); // Initialize min_energy and if_tune_energy_level for tune_energy_level in the future; 
	
	double *temp_buffer_float = new double[parameter.number_energy_level];  
	parameter.GetEnergyBound(temp_buffer_float, parameter.number_energy_level); 
	CEES_Node::SetEnergyLevels(temp_buffer_float, parameter.number_energy_level);
	parameter.GetTemperature(temp_buffer_float, parameter.number_energy_level);  
	CEES_Node::SetTemperatures(temp_buffer_float, parameter.number_energy_level);  
	delete [] temp_buffer_float; 

	CEES_Node::SetTargetAcceptanceRate(parameter.mh_target_acc);

	int *temp_buffer_int = new int[parameter.number_block];
	parameter.GetBlockSize(temp_buffer_int, parameter.number_block);  	
	CEES_Node::SetBlockSize(temp_buffer_int, parameter.number_block); 
	delete [] temp_buffer_int; 

	/*  Generate K CEES_Node objects */
	CEES_Node *simulator_node = new CEES_Node[parameter.number_energy_level]; 
	
	// Local target, pointer to the higher level node, MH Proposal distribution 
	int dim_cum_sum; 
	temp_buffer_float = new double[parameter.data_dimension]; 
 	for (int i=0; i<parameter.number_energy_level; i++)
	{
		simulator_node[i].SetID_LocalTarget(i); 	// Also set the local target distribution here
		if (i < parameter.number_energy_level -1)
			simulator_node[i].SetHigherNodePointer(simulator_node+i+1);
		else 
			simulator_node[i].SetHigherNodePointer(NULL);

		parameter.GetMHProposalScale(i, temp_buffer_float, parameter.data_dimension); 
		dim_cum_sum =0; 
		for (int iBlock = 0; iBlock <parameter.number_block; iBlock++)
		{
			simulator_node[i].SetProposal(new CTransitionModel_SimpleGaussian(parameter.GetBlockSize(iBlock), temp_buffer_float+dim_cum_sum), iBlock);
			dim_cum_sum += parameter.GetBlockSize(iBlock); 
		}
	}

	if (!if_continue)
	{
	/*
 	Initialize --> BurnIn --> loop for several times (Tune MH StepSize --> Build Initial Energy Ring (Simulate for N samples) --> TuneEnergyLevel) --> simulate
 	*/
		for (int i=parameter.number_energy_level-1; i>=0; i--)
		{
			parameter.GetCurrentState(i, temp_buffer_float, parameter.data_dimension); 
			if (i == parameter.number_energy_level-1)
				simulator_node[i].Initialize(temp_buffer_float, parameter.data_dimension); 
			else if (!simulator_node[i].Initialize(storage, r)) 
				simulator_node[i].Initialize(temp_buffer_float, parameter.data_dimension); 
			cout << "Node " << i << " Burn in for " << parameter.burn_in_period << endl; 
			simulator_node[i].BurnIn(r, storage, parameter.burn_in_period, parameter.multiple_try_mh);
			cout << "Node " << i << " MH StepSize Tuning for " << parameter.mh_stepsize_tuning_max_time << " beginning for " << parameter.mh_tracking_length << endl;  
			// simulator_node[i].MH_StepSize_Regression(parameter.mh_tracking_length, parameter.mh_stepsize_tuning_max_time, r, parameter.multiple_try_mh);	// based on regression
			simulator_node[i].MH_StepSize_Tune(parameter.mh_tracking_length, parameter.mh_stepsize_tuning_max_time, r, parameter.multiple_try_mh);	// based on Dan's adaptive strategy 
			cout << "Node " << i << " simulate for " << parameter.energy_level_tracking_window_length << endl; 
			simulator_node[i].Simulate(r, storage, parameter.energy_level_tracking_window_length, parameter.deposit_frequency, parameter.multiple_try_mh);  
		}

		// Tuning; 	
		int nEnergyLevelTuning = 0; 
		while (nEnergyLevelTuning < parameter.energy_level_tuning_max_time) 
		{
			cout << "Energy level tuning " << nEnergyLevelTuning << endl; 
			TuneEnergyLevels_UpdateStorage(simulator_node, storage, parameter);
			for (int i=parameter.number_energy_level-1; i>=0; i--)
			{
               			simulator_node[i].MH_StepSize_Regression(parameter.mh_tracking_length, parameter.mh_stepsize_tuning_max_time, r, parameter.multiple_try_mh);	// based on regression
               			// simulator_node[i].MH_StepSize_Tune(20, 10, r, MULTIPLE_TRY_MH); 	// based on Dan's adaptive strategy
                		simulator_node[i].Simulate(r, storage, parameter.energy_level_tracking_window_length, parameter.deposit_frequency, parameter.multiple_try_mh);
			}
			nEnergyLevelTuning ++; 
		}
	}
	else 
	{
		for (int i=CEES_Node::GetEnergyLevelNumber()-1; i>=0; i--)
		{
			parameter.GetCurrentState(i, temp_buffer_float, parameter.data_dimension); 
			simulator_node[i].Initialize(temp_buffer_float, parameter.data_dimension);
		}
	}
	delete [] temp_buffer_float; 
	
	cout << "Simulate for " << parameter.simulation_length<< endl; 
	// Runing through
	for (int i=CEES_Node::GetEnergyLevelNumber()-1; i>=0; i--)
		simulator_node[i].Simulate(r, storage, parameter.simulation_length, parameter.deposit_frequency, parameter.multiple_try_mh);

	storage.finalize(); 		// save to hard-disk of those unsaved data
	
	// save parameters into a binary file
	parameter.TraceStorageHead(storage);
	for (int i=0; i<CEES_Node::GetEnergyLevelNumber(); i++)
		parameter.TraceSimulator(simulator_node[i]);  
	convert.str(std::string()); 
	convert << parameter.run_id << ".parameter"; 
	file_name = storage_filename_base + convert.str();  
	parameter.SaveParameterToFile(file_name);

	// save parameters into a text file
	convert.str(std::string()); 
	convert << parameter.run_id << ".summary"; 
	file_name = storage_filename_base + convert.str(); 
	parameter.WriteSummaryFile(file_name); 

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
