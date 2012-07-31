#include <algorithm>
#include <boost/thread/xtime.hpp>
#include <boost/thread/thread.hpp>
#include <cstring>
#include "CModel.h"
#include "CBoundedModel.h"
#include "CEES_Node.h"
#include "CStorageHead.h"

int CEES_Node::K;
vector <double> CEES_Node::H; 
vector <double> CEES_Node::T; 
double CEES_Node::pee; 
int CEES_Node::N; 
int CEES_Node::dataDim; 
int CEES_Node::depositFreq; 
CModel * CEES_Node::ultimate_target;

vector <double > CEES_Node::min_energy; 
vector <double > CEES_Node::max_energy; 
int CEES_Node:: min_max_energy_capacity; 
int CEES_Node::energylevel_tracking_window_length; 
int CEES_Node::energylevel_tune_number; 
vector <bool> if_tune_energy_level; 

int CEES_Node::nBlock; 
vector < int > CEES_Node::blockSize;

int CEES_Node::max_stepsize_tune_number;
int CEES_Node::stepsize_tracking_window_length;
double CEES_Node::stepsize_tune_lower_target_prob; 
double CEES_Node::stepsize_tune_upper_target_prob; 

CEES_Node::CEES_Node(int iEnergyLevel)
{
	id = iEnergyLevel; 
	nSamplesGenerated = 0; 

	nMHSamplesAccepted_Recent = vector<int>(nBlock, 0); 
	stepsize_tune_number = 0; 
	
	proposal = new CTransitionModel *[nBlock]; 
	proposal[0] = NULL;  
	next_level = NULL; 

	x_current = new double[dataDim]; 
	x_new = new double [dataDim]; 

	target = NULL; 
	ring_size = vector <int> (K, 0); 
}

CEES_Node::CEES_Node(int iEnergyLevel, CTransitionModel *transition, CEES_Node *pNext)
{
	id = iEnergyLevel; 
	nSamplesGenerated = 0; 

	nMHSamplesAccepted_Recent = vector<int>(nBlock, 0); 
	stepsize_tune_number = 0; 

	proposal = new CTransitionModel *[nBlock]; 
	proposal[0] = transition;
	next_level = pNext;
	
	x_current = new double [dataDim]; 
	x_new = new double [dataDim];

	target = new CBoundedModel(H[id], T[id], ultimate_target); 
	ring_size = vector< int> (K, 0);
}

CEES_Node::CEES_Node(int iEnergyLevel, CTransitionModel **transition, CEES_Node *pNext)
{
        id = iEnergyLevel;
        nSamplesGenerated = 0;

        nMHSamplesAccepted_Recent = vector<int>(nBlock, 0);
	stepsize_tune_number = 0; 

	proposal = new CTransitionModel *[nBlock]; 
	for (int iBlock =0; iBlock<nBlock; iBlock++)
        	proposal[iBlock] = transition[iBlock];
        next_level = pNext;

        x_current = new double [dataDim];
        x_new = new double [dataDim];

        target = new CBoundedModel(H[id], T[id], ultimate_target);
        ring_size = vector< int> (K, 0);
}

void CEES_Node::SetID_LocalTarget(int iEnergyLevel)
{
	id = iEnergyLevel; 
	AdjustLocalTarget();	// Set up local target according to H[id] and T[id] 
}

CEES_Node::~CEES_Node()
{
	if (sizeof(x_current))
		delete [] x_current; 

	if (sizeof(x_new))
		delete [] x_new; 

	if (target != NULL)
		delete target;

	for (int iBlock =0; iBlock <nBlock; iBlock++)
		delete proposal[iBlock]; 
	delete []proposal; 
}

int CEES_Node::BinID(double e) const
{
	int energy_index = GetRingIndex(e); 
	return id*K+energy_index; 
}


bool CEES_Node::BurnInDone()
{
	if (stepsize_tune_number >= max_stepsize_tune_number && energylevel_tune_number >= max_energylevel_tune_number && nSamplesGenerated >= B)
		return true; 
	else 
		return false;
}

void CEES_Node::TuneEnergyLevel_MHStepSize_IfAppropriate()
{
	if (energylevel_tune_number < max_energylevel_tune_number)	
	{
		UpdateMinMaxEnergy(energy_current);
		if (stepsize_tune_number >= max_stepsize_tune_number && nSampelsGenerated >= energylevel_tracking_window_length)
			if_tune_energy_level[id] = true; 	
		if (TimeToTuneEnergyLevel())
		{
			TuneEnergyLevels(); 
			RewindSampleCounter(); 
			RewindStepSizeTuneCounter(); 	
		}
	}
	if (stepsize_tune_number < max_stepsize_tune_number && nSamplesGenerated%stepsize_tracking_window_length ==0) 
	{
		MH_StepSize_Tune(); 
		RewindSampleCounter(); 
	}
}


void CEES_Node::Initialize(const gsl_rng *r)
{
	target->draw(x_current, dataDim, r); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current);

	nSamplesGenerated ++; 
	TuneEnergyLevel_MHStepSize_IfAppropriate(); 
}

void CEES_Node::Initialize(CModel * model, const gsl_rng *r)
{
	if (model == NULL)
		Initialize(r); 
	else
		model->draw(x_current, dataDim, r); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current); 
	nSamplesGenerated ++;
	TuneEnergyLevel_MHStepSize_IfAppropriate();
}

void CEES_Node::Initialize(const double *x, int x_d)
{
	memcpy(x_current, x, sizeof(double)*dataDim); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current); 
	nSamplesGenerated++;
	TuneEnergyLevel_MHStepSize_IfAppropriate();	
}

bool CEES_Node::SetEnergyLevels(double *e, int n)
{
	if (K > n)
		return false; 
	H = vector<double>(K);
	for (int i=0; i<K; i++)
		H[i] = e[i]; 
	return true;
}

bool CEES_Node::SetEnergyLevels_GeometricProgression(double H0, double HK_1)
{
        /*
 *     H[i] = H[i-1]+gamma^i
 *     gamma is determined by solving a polynomial equation 
 *     gamma+gamma^2+...+gamma^{K-1} = H[K-1]-H[0]; 
 *     */
        double *coefficients = new double [K];
        coefficients[0] = H0-HK_1;
        for (int i=1; i<K; i++)
                coefficients[i]=1;
        double *Z = new double [(K-1)*2];

        gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(K);
        gsl_poly_complex_solve(coefficients, K, w, Z);

        double gamma;
        bool continue_flag = true;
        for (int i=0; i<K-1 && continue_flag; i++)
        {
                if (Z[2*i]>0 && abs(Z[2*i+1]) <= DBL_EPSILON)
                {
                        gamma = Z[2*i];
                        continue_flag = false;
                }
        }
        delete [] Z;
        delete [] coefficients;
        if (continue_flag)
                return false;
        /* END: solving the polynomial equation*/

        H = vector<double>(K);
	H[0] = H0; 
	H[K-1] = HK_1;
        for (int i=1; i<K-1; i++)
                H[i] = H[i-1]+pow(gamma, i);
        return true;
}

bool CEES_Node::SetTemperatures(double* t, int n)
{
	if (K >n)
		return false; 
	T = vector<double>(K); 
	for (int i=0; i<K; i++)
		T[i] = t[i];
	return true;
} 

bool CEES_Node::SetTemperatures_EnergyLevels(double T0, double TK_1, double c)
{
	T = vector<double>(K);
	T[0] = T0; 
	T[K-1] = TK_1; 
	for (int i=1; i<K-1; i++)
		T[i] = (H[i+1]-H[i])/c;
	return true;
}

bool CEES_Node::SetTemperatures_EnergyLevels(double T0, double TK_1)
{
	T = vector<double > (K); 
	T[K-1] = TK_1; 
	T[0] = T0; 
	double gamma = (T[K-1]-T[0])/(H[K-1]-H[0]);
	for (int i=1; i<K-1; i++)
		T[i] = gamma * (H[i]-H[0]); 
	return true;  
	
}

int CEES_Node::GetRingIndex(double e) const
{
	for (int j=1; j<K; j++)
	{
		if (e < H[j])
			return j-1; 
	}	
	return K-1;
}

bool CEES_Node::EnergyRingBuildDone() const
{
	if (BurnInDone() && nSamplesGenerated >= B+N)
		return true; 
	else 
		return false;
}

void CEES_Node::draw(const gsl_rng *r, CStorageHead &storage, int mMH )
{
	int bin_id_next_level, bin_id, x_id; 
	bool new_sample_flag;
	double x_weight; 
	if (next_level != NULL)
		bin_id_next_level = next_level->BinID(ring_index_current); 
	if (next_level == NULL || storage.empty(bin_id_next_level)) // MH draw
	{
		target->draw(proposal[0], x_new, dataDim, x_current, r, new_sample_flag, mMH); 
		if (new_sample_flag)
		{
			memcpy(x_current, x_new, sizeof(x_new)*dataDim); 
			energy_current = OriginalEnergy(x_current, dataDim); 
			ring_index_current = GetRingIndex(energy_current); 

			if (energylevel_tune_number < max_energylevel_tune_number)
				UpdateMinMaxEnergy(energy_current); 
		}
		if (stepsize_tune_number < max_stepsize_tune_number && new_sample_flag)
			nMHSamplesAccepted_Recent[0] ++; 
	}
	else	// equi-energy draw with prob. of pee
	{
		double uniform_draw = gsl_rng_uniform(r); 
		if (uniform_draw <= pee && storage.DrawSample(bin_id_next_level, x_new, dataDim, x_id, x_weight, r))
		{ 
			/*double ratio=ProbabilityRatio(x_new, x_current, dataDim); 
			ratio = ratio * next_level->ProbabilityRatio(x_current, x_new, dataDim);  */
			// need to use LogProbRatio
			double ratio = LogProbRatio(x_new, x_current, dataDim); 
			ratio += next_level->LogProbRatio(x_current, x_new, dataDim); 
			double another_uniform_draw = gsl_rng_uniform(r); 
			if (log(another_uniform_draw) <= ratio)
				memcpy(x_current, x_new, sizeof(x_new)*dataDim);
		} 
		else	// MH draw 
		{
			target->draw(proposal[0], x_new, dataDim, x_current, r, new_sample_flag, mMH); 
			if (new_sample_flag)
			{
				memcpy(x_current, x_new, sizeof(x_new)*dataDim); 
				energy_current = OriginalEnergy(x_current, dataDim); 
				ring_index_current = GetRingIndex(energy_current); 
				if (energylevel_tune_number < max_energylevel_tune_number)
					UpdateMinMaxEnergy(energy_current); 
			}
                	if (stepsize_tune_number < max_stepsize_tune_number && new_sample_flag)
                        	nMHSamplesAccepted_Recent[0] ++;
		}
	}

	if (BurnInDone() && nSamplesGenerated%depositFreq==0)
	{
		bin_id = BinID(ring_index_current); 
		x_id = nSamplesGenerated; 
		x_weight = 1.0; 
		storage.DepositSample(bin_id, x_current, dataDim, x_id, x_weight); 
		ring_size[ring_index_current] ++;
	}
	nSamplesGenerated ++;
	TuneEnergyLevel_MHStepSize_IfAppropriate();
}

void CEES_Node::draw_block(const gsl_rng *r, CStorageHead &storage)
{
	int bin_id_next_level, bin_id, x_id; 
	vector <bool > new_sample_flag(nBlock, false); 
	double overall_new_sample_flag; 
	double x_weight; 
	if (next_level != NULL)
		bin_id_next_level = next_level->BinID(ring_index_current); 
	if (next_level == NULL || storage.empty(bin_id_next_level)) 
	{
		target->draw(proposal, x_new, dataDim, x_current, r, new_sample_flag, nBlock, blockSize); 
		// start: check each block to see if it has been updated 
		overall_new_sample_flag = false; 
		for (int iBlock =0; iBlock <nBlock; iBlock ++)
		{
			if(new_sample_flag[iBlock])
			{
				overall_new_sample_flag = true; 
				if (stepsize_tune_number < max_stepsize_tune_number)
					nMHSamplesAccepted_Recent[iBlock] ++;
			}
		}
		if (overall_new_sample_flag)
		{
			memcpy(x_current, x_new, sizeof(x_new)*dataDim); 
			energy_current = OriginalEnergy(x_current, dataDim);
			ring_index_current = GetRingIndex(energy_current); 
			if (energylevel_tune_number < max_energylevel_tune_number)
				UpdateMinMaxEnergy(energy_current); 
		}
		// end: check each block to see if it has been updated
	}
	else	// equi-energy draw with prob. of pee
	{
		double uniform_draw = gsl_rng_uniform(r); 
		if (uniform_draw <= pee && storage.DrawSample(bin_id_next_level, x_new, dataDim, x_id, x_weight, r))
		{ 
			/*double ratio=ProbabilityRatio(x_new, x_current, dataDim); 
			ratio = ratio * next_level->ProbabilityRatio(x_current, x_new, dataDim);  */
			// need to use LogProbRatio
			double ratio = LogProbRatio(x_new, x_current, dataDim); 
			ratio += next_level->LogProbRatio(x_current, x_new, dataDim); 
			double another_uniform_draw = gsl_rng_uniform(r); 
			if (log(another_uniform_draw) <= ratio)
				memcpy(x_current, x_new, sizeof(x_new)*dataDim);
		} 
		else	// Block MH
		{
			target->draw(proposal, x_new, dataDim, x_current, r, new_sample_flag, nBlock, blockSize); 
			overall_new_sample_flag = false;
                	for (int iBlock =0; iBlock <nBlock; iBlock ++)
                	{
                        	if(new_sample_flag[iBlock])
                        	{
                                	overall_new_sample_flag = true;
                                	if (stepsize_tune_number < max_stepsize_tune_number)
                                        	nMHSamplesAccepted_Recent[iBlock] ++;
                        	}
                	}
                	if (overall_new_sample_flag)
                	{
                        	memcpy(x_current, x_new, sizeof(x_new)*dataDim);
                        	energy_current = OriginalEnergy(x_current, dataDim);
                        	ring_index_current = GetRingIndex(energy_current);
				if (energylevel_tune_number < max_energylevel_tune_number)
                        		UpdateMinMaxEnergy(energy_current);
                	}
		}
	}

	if (BurnInDone() && nSamplesGenerated%depositFreq==0)
	{
		bin_id = BinID(ring_index_current); 
		x_id = nSamplesGenerated; 
		x_weight = 1.0; 
		storage.DepositSample(bin_id, x_current, dataDim, x_id, x_weight); 
		ring_size[ring_index_current] ++;
	}
	nSamplesGenerated ++;
	TuneEnergyLevel_MHStepSize_IfAppropriate();
}

ofstream & summary(ofstream &of, const CEES_Node *simulator)
{
	of << "Data Dimension:\t" << CEES_Node::dataDim << endl; 
	of << "Number of Energy Levels:\t" << CEES_Node::K << endl; 
	of << "Energy Thresholds:"; 
	for (int i=0; i<CEES_Node::K; i++)
		of << "\t" << CEES_Node::H[i];
	of << endl; 
	of << "Temperatures:"; 
	for (int i=0; i<CEES_Node::K; i++)
		of << "\t" << CEES_Node::T[i]; 
	of << endl; 
	of << "Prob Equi-Jump:\t" << CEES_Node::pee << endl; 
	of << "Build Initial Ring:\t" << CEES_Node::N << endl; 
	of << "Deposit Frequency:\t" << CEES_Node::depositFreq << endl; 
	return of;
}

void CEES_Node::MH_StepSize_Tune()
// Adapt from Dan's Adaptive Scaling
{
	double log_mid = log(stepsize_tune_target_prob); 
	double lower_bound = exp(log_mid/0.2); 
	double upper_bound = exp(log_mid/5.0);  
	
	double nSamplesGenerated = stepsize_tracking_window_length;
	double nSamplesAccepted; 
	double previous_ratio, current_ratio; 
	double low_scale, high_scale, scale; 
	double low_jump_ratio, high_jump_ratio;
	
	for (int iBlock=0; iBlock<nBlock; iBlock++)
	{
		nSamplesAccepted = nMHSamplesAccepted_Recent[iBlock];
		previous_ratio = (double)(nSamplesAccepted)/(double)(nSamplesGenerated); 
		if (previous_ratio < stepsize_tune_target_prob)
		{
			low_scale = 
		}
		{
		}	
	}

	/*double accept_ratio; 
	for (int iBlock = 0; iBlock < nBlock; iBlock ++)
	{
		accept_ratio = (double)(nMHSamplesAccepted_Recent[iBlock])/nMHSamplesGenerated_Recent; 
		if (accept_ratio < MH_lower_target_prob)	// should decrease stepsize
			proposal[iBlock]->step_size_tune(accept_ratio/MH_lower_target_prob < 0.5 ? 0.5 : accept_ratio/MH_lower_target_prob); 
		else if (accept_ratio > MH_upper_target_prob ) 	// should increase stepsize
			proposal[iBlock]->step_size_tune(accept_ratio/MH_upper_target_prob > 2 ? 2: accept_ratio/MH_upper_target_prob); 
	} */
}

void CEES_Node::InitializeMinMaxEnergy(int capacity)
{
	min_max_energy_capacity = capacity; 	
}

void CEES_Node::UpdateMinMaxEnergy(double _new_energy)
{
	vector<double>::iterator position; 
	position =upper_bound(min_energy.begin(), min_energy.end(), _new_energy); 
	min_energy.insert(position, _new_energy); 
	if ((int)(min_energy.size()) > min_max_energy_capacity)
		min_energy.pop_back(); 	
	
	position =lower_bound(max_energy.begin(), max_energy.end(), _new_energy); 
	max_energy.insert(position, _new_energy); 
	if ((int)(max_energy.size()) > min_max_energy_capacity)
		max_energy.erase(max_energy.begin()); 
}

void CEES_Node::AdjustLocalTarget()
{
	if (target)
		delete target; 
	target = new CBoundedModel(H[id], T[id], ultimate_target); 
}

void CEES_Node::SetBlockSize(const int* _bSize, int _nB)
{
	nBlock = _nB; 
	blockSize.resize(nBlock);
	if (nBlock == 1)
		blockSize[0] = dataDim;
	else if (nBlock == dataDim)
	{
		for (int i=0; i<nBlock; i++)
			blockSize[i] = 1; 
	}
	else
	{
		for (int i=0; i<nBlock; i++)
			blockSize[i] = _bSize[i]; 
	}
}

void CEES_Node::AssignSamplesGeneratedSoFar(CStorageHead &storage)
{
	// current sample
	ring_index_current = GetRingIndex(energy_current); 		
	vector <int> old_ring_size = ring_size; 
	ring_size = vector <int> (ring_size.size(), 0);  
	ring_size[ring_index_current] ++; 

	// samples from storage
	int storage_bin_id, new_storage_bin_id; 
	vector <CSampleIDWeight> samples; 
	double energy; 
	int ring_index, id;
	double  weight; 
	for (int bin_id = 0; bin_id < K; bin_id ++)
	{
		storage_bin_id = BinID(bin_id); 
		samples=storage.RetrieveSamplesSequentially(true, storage_bin_id); 
		while (!samples.empty())
		{
			for (int i=0; i<(int)(samples.size()); i++)
			{
				samples[i].GetData(x_new, dataDim, id, weight);
				energy = OriginalEnergy(x_new, dataDim); 
				ring_index = GetRingIndex(energy); 
				new_storage_bin_id = BinID(ring_index); 
				storage.DepositSample(true, new_storage_bin_id, x_new, dataDim, id, weight); 	
				ring_size[ring_index] ++;
			}
			samples=storage.RetrieveSamplesSequentially(true, storage_bin_id);
		}
	}
	for (int bin_id =0; bin_id <K; bin_id ++)
	{
		storage_bin_id = BinID(bin_id); 
		storage.Consolidate(storage_bin_id); 
	}
}

void CEES_Node:: SetMaxEnergyLevelTuneNumber(int _n) 
{ 
	max_energylevel_tune_number = _n; 
	energy_level_tune_number = 0;
	if_tune_energy_level = vector <bool> (K, false); 
}

bool CEES_Node::TimeToTuneEnergyLevel()
{
	for (int i=0; i<(int)(if_tune_energy_level.size()); i++)
	{
		if (!if_tune_energy_level[i])
			return false; 
	}
	return true; 
}
