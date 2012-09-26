#include <algorithm>
#include <cstring>
#include "CModel.h"
#include "CBoundedModel.h"
#include "CEES_Node.h"
#include "CStorageHead.h"
#include "MHAdaptive.h"
#include "CTransitionModel_SimpleGaussian.h"

int CEES_Node::K;
vector <double> CEES_Node::H; 
vector <double> CEES_Node::T; 
double CEES_Node::pee; 
int CEES_Node::dataDim; 
// CModel * CEES_Node::ultimate_target;
vector <double> CEES_Node::targetAcc;
int CEES_Node:: min_max_energy_capacity; 

int CEES_Node::nBlock; 
vector < int > CEES_Node::blockSize;

CEES_Node::CEES_Node(int iEnergyLevel, CTransitionModel *transition, CEES_Node *pNext)
{
	id = iEnergyLevel; 
	nSamplesGenerated = 0; 

	proposal = new CTransitionModel *[nBlock]; 
	proposal[0] = transition;
	next_level = pNext;
	
	x_current.SetDataDimension(dataDim); 
	x_new.SetDataDimension(dataDim); 

	// target = new CBoundedModel(H[id], T[id], ultimate_target); 
	// Because ultimate_target is not a static variable, 
	// it must be set explicitly.
	// Therefore, ultimate_target and target will be set outside the 
	// construction function.
	ultimate_target = NULL; 
	target = NULL; 
	
	ring_size = vector< int> (K, 0);
}

CEES_Node::CEES_Node(int iEnergyLevel, CTransitionModel **transition, CEES_Node *pNext)
{
        id = iEnergyLevel;
        nSamplesGenerated = 0;

	proposal = new CTransitionModel *[nBlock]; 
	for (int iBlock =0; iBlock<nBlock; iBlock++)
        	proposal[iBlock] = transition[iBlock];
        next_level = pNext;

        x_current.SetDataDimension(dataDim);
        x_new.SetDataDimension(dataDim); 

        // target = new CBoundedModel(H[id], T[id], ultimate_target);
        // Because ultimate_target is not a static variable, 
        // it must be set explicitly.
        // Therefore, ultimate_target and target will be set outside the 
        // construction function.
       
        ultimate_target = NULL; 
	target = NULL;  
        ring_size = vector< int> (K, 0);
}

void CEES_Node::SetID_LocalTarget(int iEnergyLevel)
{
	id = iEnergyLevel; 
	AdjustLocalTarget();	// Set up local target according to H[id] and T[id] 
}

CEES_Node::~CEES_Node()
{
	if (target)
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

void CEES_Node::Initialize(const gsl_rng *r, CModel *model)
{
	bool if_new_sample; 
	if (model == NULL)	// directly draw a sample from target
		target->draw(x_current, if_new_sample, r);  
	else
	{
		model->draw(x_current, if_new_sample, r); 
		target->log_prob(x_current); 
	}
	ring_index_current = GetRingIndex(x_current.GetWeight()); 
	UpdateMinMaxEnergy(x_current.GetWeight()); 
}

void CEES_Node::Initialize(const double *x, int x_d)
{
	x_current = CSampleIDWeight(x, x_d); 
	target->log_prob(x_current); 
	ring_index_current = GetRingIndex(x_current.GetWeight()); 
	UpdateMinMaxEnergy(x_current.GetWeight()); 
}

void CEES_Node::Initialize(const CSampleIDWeight &x)
{
	x_current = x; 
	target->log_prob(x_current); 
	ring_index_current = GetRingIndex(x_current.GetWeight()); 
	UpdateMinMaxEnergy(x_current.GetWeight()); 
}

bool CEES_Node::Initialize(CStorageHead &storage, const gsl_rng *r)
{
	// Initialize using samples from the next level; 
	if (next_level == NULL)
                return false;
        int bin_id_next_level; 
	// Try next levels' bins with the same or lower energies
	for (int try_id = id; try_id >= 0; try_id --)
	{
		bin_id_next_level = next_level->BinID(try_id); 
        	if (!storage.empty(bin_id_next_level) && storage.DrawSample(bin_id_next_level, r, x_current))
		{
			// x_current.weight will remain the same 
			// x_current.log_prob needs to be updated according to 
			// current level's H and T
			x_current.log_prob = -(x_current.GetWeight() > GetEnergy() ? x_current.GetWeight() : GetEnergy())/GetTemperature(); 
			ring_index_current = GetRingIndex(x_current.GetWeight());  
			UpdateMinMaxEnergy(x_current.GetWeight()); 
			return true;
		}
	}
	// If not successful, then try next level's bins with higher energies
	for (int try_id = id+1; try_id <K; try_id ++)
	{
		bin_id_next_level = next_level->BinID(try_id);
                if (!storage.empty(bin_id_next_level) && storage.DrawSample(bin_id_next_level, r, x_current))
                {
			// x_current.weight will remain the same 
			// x_current.log_prob needs to be updated according to
			// current level's H and  T
			x_current.log_prob = -(x_current.GetWeight() > GetEnergy() ? x_current.GetWeight() : GetEnergy())/GetTemperature(); 
                        ring_index_current = GetRingIndex(x_current.GetWeight());
			UpdateMinMaxEnergy(x_current.GetWeight()); 
                        return true;
                }
	}
	return false; 
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

bool CEES_Node::SetTemperatures(double* t, int n)
{
	if (K >n)
		return false; 
	T = vector<double>(K); 
	for (int i=0; i<K; i++)
		T[i] = t[i];
	return true;
} 

bool CEES_Node::SetTargetAcceptanceRate(double p0)
{
	targetAcc = vector<double>(K, p0); 
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


bool CEES_Node::MH_draw(const gsl_rng *r, int mMH)
{
	vector <bool> new_sample_flag(nBlock, false); 
	if (nBlock <= 1)
	{
		bool local_flag;
		target->draw(x_new, proposal[0], local_flag, r, x_current, mMH);
		new_sample_flag[0] = local_flag;
	}
	else 
		target->draw(x_new, proposal, new_sample_flag, r, x_current, nBlock, blockSize, mMH);
	bool overall_new_sample_flag = false;
	int iBlock = 0; 
	while (iBlock < nBlock && !overall_new_sample_flag)
	{ 
		if (new_sample_flag[iBlock])
			overall_new_sample_flag = true; 
		iBlock ++; 
	}
	if (overall_new_sample_flag)
	{
		x_current = x_new;  
		ring_index_current = GetRingIndex(x_current.GetWeight());
		UpdateMinMaxEnergy(x_current.GetWeight()); 
	}
	return overall_new_sample_flag; 
}

bool CEES_Node::EE_draw(const gsl_rng *r, CStorageHead &storage)
{
	if (next_level == NULL)
		return false; 
	int bin_id_next_level = next_level->BinID(ring_index_current); 
	if (storage.empty(bin_id_next_level) || !storage.DrawSample(bin_id_next_level, r, x_new))
		return false; 
	
	// if x_new is adopted
	// x_new.weight can be retained
	// x_new.log_prob needs to be updated to reflect current level's H and T

	double logRatio = LogProbRatio_Energy(x_new.GetWeight(), x_current.GetWeight()); 
	logRatio += next_level->LogProbRatio_Energy(x_current.GetWeight(), x_new.GetWeight()); 
	
	double uniform_draw = gsl_rng_uniform(r); 
	if (log(uniform_draw) <= logRatio)
	{
		x_current = x_new;  
		x_current.log_prob = -(x_current.GetWeight() > GetEnergy() ? x_current.GetWeight() : GetEnergy())/GetTemperature();	
		return true; 
	}
	else 
		return false; 
}

bool CEES_Node::draw(const gsl_rng *r, CStorageHead &storage, int mMH )
{
	bool new_sample_flag; 
	double uniform_draw; 
	if (next_level == NULL )	 // MH draw
		new_sample_flag = MH_draw(r, mMH); 
	else	// equi-energy draw with prob. of pee
	{
		uniform_draw = gsl_rng_uniform(r); 
		if (uniform_draw > pee ||  !(new_sample_flag = EE_draw(r, storage)) )
			new_sample_flag = MH_draw(r, mMH); 
	}
	return new_sample_flag; 
}

void CEES_Node::BurnIn(const gsl_rng *r, CStorageHead &storage, int B, int mMH)
// nSamplesGenerated does not grow during the burn-in period
{
	for (int i=0; i<B; i++)
		draw(r, storage, mMH); 
}

void CEES_Node::MH_StepSize_Tune(int initialPeriodL, int periodNumber, const gsl_rng *r, int mMH)
// Adapt from Dan's Adaptive Scaling
{
	// Save current state, because (1) tuning is based on a mode, and (2) after tuning is done, simulator will resume from the current state
	CSampleIDWeight x_current_saved = x_current; 

	ultimate_target->GetMode(x_current); 
	// At this moment x_current.weight and x_current.log_pron are wrt ultimate_target
	// x_current.weight will remain the same
	// x_current.log_prob needs to be updated to reflect this level's H and T
	x_current.log_prob = -(x_current.GetWeight() > GetEnergy() ? x_current.GetWeight() : GetEnergy())/GetTemperature(); 
	int dim_lum_sum=0; 
	for (int iBlock=0; iBlock<nBlock; iBlock++)
	{
		proposal[iBlock]->Tune(targetAcc[this->id], initialPeriodL, periodNumber, r, target, x_current, dim_lum_sum, blockSize[iBlock]); 
		dim_lum_sum += blockSize[iBlock];  
	}
	
	/*int nGenerated = initialPeriodL; // length of observation
	int nAccepted;  
	
	MHAdaptive *adaptive; 
	CTransitionModel_SimpleGaussian *individual_proposal = new CTransitionModel_SimpleGaussian(1); 
	  
	// Tune for a number of perioldNumber times
	int dim_lum_sum =0;
	bool new_sample_flag = false; 
	ultimate_target->GetMode(x_current, dataDim); 
	double log_prob_mode = target->log_prob(x_current, dataDim);  
	double log_prob_x; 
	for (int iBlock=0; iBlock<nBlock; iBlock++)
	{
		for (int iDim=0; iDim<blockSize[iBlock]; iDim++)
		{
			adaptive = new MHAdaptive(targetAcc[this->id], proposal[iBlock]->get_step_size(iDim)); 
			individual_proposal->set_step_size(proposal[iBlock]->get_step_size(iDim));
			for (int nPeriod=0; nPeriod<periodNumber; nPeriod++)
			{	
				// tuning starts from a mode of the target distribution
				ultimate_target->GetMode(x_current, dataDim); 
				nAccepted = 0;
				log_prob_x = log_prob_mode;  
				// draw samples 
				for (int iteration=0; iteration<nGenerated; iteration ++)
				{
					log_prob_x = target->draw_block(dim_lum_sum+iDim, 1, individual_proposal, x_new, dataDim, new_sample_flag, r, x_current,log_prob_x, mMH); 
					if (new_sample_flag)
					{
						memcpy(x_current, x_new, dataDim*sizeof(double));
						nAccepted ++; 
					}
				}
				// Update scale 
				if(adaptive->UpdateScale(nGenerated, nAccepted)) 
					individual_proposal->set_step_size(adaptive->GetScale()); 
			}
			proposal[iBlock]->set_step_size(adaptive->GetBestScale(), iDim); 
		}
		delete adaptive;
		dim_lum_sum += blockSize[iBlock]; 
	}
	
	delete individual_proposal; */
	// restore x_current
	x_current = x_current_saved; 
}

void CEES_Node::Simulate(const gsl_rng *r, CStorageHead &storage, int N, int depositFreq, int mMH)
{
	int bin_id;  
	for (int i=0; i<N; i++)
	{
		draw(r, storage, mMH); 
		if (i%depositFreq == 0)	// Deposit sample
		{
			bin_id = BinID(ring_index_current); 
			x_current.SetID(nSamplesGenerated); 
			storage.DepositSample(bin_id, x_current); 
			ring_size[ring_index_current] ++; 
		}
		nSamplesGenerated++;	
	}
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

void CEES_Node::DisregardHistorySamples(CStorageHead &storage)
{
	// keep current sample
	ring_index_current = GetRingIndex(x_current.GetWeight());
        vector <int> old_ring_size = ring_size;
        ring_size = vector <int> (ring_size.size(), 0);
        ring_size[ring_index_current] ++;
	
	// clear samples from storage
	int storage_bin_id; 
	for (int bin_id =0; bin_id <K; bin_id++)
	{
		storage_bin_id = BinID(bin_id); 
		storage.DisregardHistorySamples(storage_bin_id); 
		storage.ClearDepositDrawHistory(storage_bin_id); 
	}
}

void CParameterPackage::TraceSimulator(const CEES_Node &simulator)
{
	int id = simulator.id; 
	if (id == 0)
	{
		if ((int)h.size() < CEES_Node::K)
			h.resize(CEES_Node::K);
       		for (int i=0; i<CEES_Node::K; i++)
			h[i] = CEES_Node::H[i];

		if ((int)t.size() < CEES_Node::K)
			t.resize(CEES_Node::K);
		for (int i=0; i<CEES_Node::K; i++)
       			t[i] = CEES_Node::T[i];
	}
        
	if ((int)x_current.size() < CEES_Node::K)
       		x_current.resize(CEES_Node::K) ;
	//if ((int)energy_index_current.size() < CEES_Node::K)
       	//	energy_index_current.resize(CEES_Node::K);
	if ((int)scale.size() < CEES_Node::K)
       		scale.resize(CEES_Node::K); 

        // energy_index_current[id] = simulator.ring_index_current;
	if ((int)scale[id].size() < CEES_Node::dataDim)
		scale[id].resize(CEES_Node::dataDim); 
	int dim_cum_sum =0; 
       	for (int iBlock=0; iBlock<CEES_Node::nBlock; iBlock++)
       	{
               	for (int j=0; j<CEES_Node::blockSize[iBlock]; j++)
                       	scale[id][j+dim_cum_sum] = simulator.proposal[iBlock]->get_step_size(j);
		dim_cum_sum += CEES_Node::blockSize[iBlock]; 
       	}
	x_current[id] = simulator.x_current; 
}

double CEES_Node::LogProbRatio_Energy(double energy_x, double energy_y)
{
	double log_prob_x_bounded = -(energy_x > GetEnergy() ? energy_x : GetEnergy())/GetTemperature(); 
	double log_prob_y_bounded = -(energy_y > GetEnergy() ? energy_y : GetEnergy())/GetTemperature(); 
	return log_prob_x_bounded - log_prob_y_bounded; 
}
