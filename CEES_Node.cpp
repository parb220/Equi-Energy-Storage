#include <algorithm>
#include <cstring>
#include "CModel.h"
#include "CBoundedModel.h"
#include "CEES_Node.h"
#include "CStorageHead.h"
#include "MHAdaptive.h"

int CEES_Node::K;
vector <double> CEES_Node::H; 
vector <double> CEES_Node::T; 
double CEES_Node::pee; 
int CEES_Node::dataDim; 
CModel * CEES_Node::ultimate_target;
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
	
	x_current = new double [dataDim]; 
	x_new = new double [dataDim];

	target = new CBoundedModel(H[id], T[id], ultimate_target); 
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
	if (x_current)
		delete [] x_current; 

	if (x_new)
		delete [] x_new; 

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

void CEES_Node::Initialize(const gsl_rng *r)
{
	target->draw(x_current, dataDim, r); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current);
}

void CEES_Node::Initialize(CModel * model, const gsl_rng *r)
{
	if (model == NULL)
		Initialize(r); 
	else
		model->draw(x_current, dataDim, r); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current); 
}

void CEES_Node::Initialize(const double *x, int x_d)
{
	memcpy(x_current, x, sizeof(double)*dataDim); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current); 
}

bool CEES_Node::Initialize(CStorageHead &storage, const gsl_rng *r)
{
	if (next_level == NULL)
                return false;
        int bin_id_next_level; 
        int x_id;
        double x_weight;
	for (int try_id = id; try_id >= 0; try_id --)
	{
		bin_id_next_level = next_level->BinID(try_id); 
        	if (!storage.empty(bin_id_next_level))
		{
			storage.DrawSample(bin_id_next_level, x_new, dataDim, x_id, x_weight, r);
			memcpy(x_current, x_new, sizeof(double)*dataDim);
			energy_current = OriginalEnergy(x_current, dataDim); 
			ring_index_current = GetRingIndex(energy_current);  
			return true;
		}
	}
	for (int try_id = id+1; try_id <K; try_id ++)
	{
		bin_id_next_level = next_level->BinID(try_id);
                if (!storage.empty(bin_id_next_level))
                {
                        storage.DrawSample(bin_id_next_level, x_new, dataDim, x_id, x_weight, r);
                        memcpy(x_current, x_new, sizeof(double)*dataDim);
                        energy_current = OriginalEnergy(x_current, dataDim);
                        ring_index_current = GetRingIndex(energy_current);
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

bool CEES_Node::SetTemperatures_EnergyLevels(double T0, double c, bool flag)
{
	T = vector <double>(K); 
	T[0] = T0; 
	for (int i=1; i<K; i++)
		T[i] = T[i-1]+(H[i]-H[i-1])/c; 
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
		target->draw(proposal[0], x_new, dataDim, x_current, r, local_flag, mMH);
		new_sample_flag[0] = local_flag;
	}
	else 
		target->draw(proposal, x_new, dataDim, x_current, r, new_sample_flag, nBlock, blockSize, mMH);
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
		memcpy(x_current, x_new, sizeof(double)*dataDim); 
		energy_current = OriginalEnergy(x_current, dataDim); 
		ring_index_current = GetRingIndex(energy_current);
		UpdateMinMaxEnergy(energy_current); 
	}
	return overall_new_sample_flag; 
}

bool CEES_Node::EE_draw(const gsl_rng *r, CStorageHead &storage)
{
	if (next_level == NULL)
		return false; 
	int bin_id_next_level = next_level->BinID(ring_index_current); 
	if (storage.empty(bin_id_next_level))
		return false; 
	int x_id; 
	double x_weight; 
	storage.DrawSample(bin_id_next_level, x_new, dataDim, x_id, x_weight, r); 	

	double logRatio = LogProbRatio(x_new, x_current, dataDim);
	logRatio += next_level->LogProbRatio(x_current, x_new, dataDim);
	
	double uniform_draw = gsl_rng_uniform(r); 
	if (log(uniform_draw) <= logRatio)
	{
		memcpy(x_current, x_new, sizeof(double)*dataDim); 	
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
	int nPeriod; // number of periods of observation
	int nGenerated = initialPeriodL; // length of observation
	int nAccepted;  
	
	MHAdaptive **adaptive = new MHAdaptive *[nBlock];  

	// Initialize	
	for (int iBlock=0; iBlock<nBlock; iBlock++)
		adaptive[iBlock] = new MHAdaptive(targetAcc[this->id], proposal[iBlock]->get_step_size());
	
	// Tune for a number of perioldNumber times,
	nPeriod = 0; 	
	while(nPeriod < periodNumber)
	{
		// MH draw for  
		nAccepted = 0; 
		for (int iteration =0; iteration < nGenerated; iteration ++)
		{
			if (MH_draw(r, mMH))
				nAccepted ++; 
		}

		// Update Scale
		for (int iBlock=0; iBlock<nBlock; iBlock++)
		{
			if (adaptive[iBlock]->UpdateScale(nGenerated, nAccepted)) 
				proposal[iBlock]->set_step_size(adaptive[iBlock]->GetScale()); 
		}
		//nGenerated *=2; 
		nPeriod ++; 
	}
	for (int iBlock =0; iBlock < nBlock; iBlock++)
	{
		if (fabs(adaptive[iBlock]->GetBestScale()-adaptive[iBlock]->GetScale()) > 1.0e-6)
			proposal[iBlock]->set_step_size(adaptive[iBlock]->GetBestScale()); 
		delete adaptive[iBlock]; 
	}
	delete [] adaptive; 
}

void CEES_Node::MH_StepSize_Regression(int periodL, int periodNumber, const gsl_rng *r, int mMH)
{
	vector <AccStep> observation(periodNumber); 
	double log_increment; 
	double log_sn, log_min, log_max; 

	MHAdaptive *adaptive;  
	double *target_step_size = new double[nBlock]; 
       	for (int iBlock=0; iBlock<nBlock; iBlock++)
       	{
       		adaptive = new MHAdaptive(targetAcc[this->id], proposal[iBlock]->get_step_size());
		log_sn = log(proposal[iBlock]->get_step_size());
		log_min = log_sn - log(1.0e3) > log(1.0e-3) ? log_sn - log(1.0e3) : log(1.0e-3); 
		log_max = log_sn + log(1.0e3) < log(1.0e3) ? log_sn + log(1.0e3) : log(1.0e3);
		log_increment = (log_max-log_min)/(periodNumber-1); 
		for (int nPeriod=0; nPeriod<periodNumber; nPeriod ++)
		{
			observation[nPeriod].nGenerated = periodL;
			observation[nPeriod].nAccepted = 0; 
			observation[nPeriod].step = log_min+log_increment*nPeriod;   // log of step_size
			proposal[iBlock]->set_step_size(exp(observation[nPeriod].step)); 
			for (int iteration =0; iteration < periodL; iteration ++)
			{
				if (MH_draw(r, mMH))
					observation[nPeriod].nAccepted ++; 
			}
			observation[nPeriod].acc = (double)(observation[nPeriod].nAccepted)/(double)(observation[nPeriod].nGenerated);
		}
		adaptive->EstimateRegressionParameters(observation); 
		target_step_size[iBlock] = adaptive->GetStepSizeRegression(); 
		proposal[iBlock]->set_step_size(target_step_size[iBlock]); 
	       	delete adaptive;
	}
	delete [] target_step_size; 
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
			storage.DepositSample(bin_id, x_current, dataDim, nSamplesGenerated, 1.0); 
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
	ring_index_current = GetRingIndex(energy_current);
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
	if ((int)energy_index_current.size() < CEES_Node::K)
       		energy_index_current.resize(CEES_Node::K);
	if ((int)scale.size() < CEES_Node::K)
        	scale.resize(CEES_Node::K); 

	if ((int)x_current[id].size() < CEES_Node::dataDim)
		x_current[id].resize(CEES_Node::dataDim); 
	for (int j=0; j<CEES_Node::dataDim; j++)
		x_current[id][j] = simulator.x_current[j]; 
        energy_index_current[id] = simulator.ring_index_current;
	if ((int)scale[id].size() < CEES_Node::dataDim)
		scale[id].resize(CEES_Node::dataDim); 
	int dim_cum_sum =0; 
        for (int iBlock=0; iBlock<CEES_Node::nBlock; iBlock++)
        {
                for (int j=0; j<CEES_Node::blockSize[iBlock]; j++)
                        scale[id][j+dim_cum_sum] = simulator.proposal[iBlock]->get_step_size();
		dim_cum_sum += CEES_Node::blockSize[iBlock]; 
        }
}

