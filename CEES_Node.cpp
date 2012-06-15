#include <cstring>
#include "CModel.h"
#include "CBoundedModel.h"
#include "sdsm.h"
#include "CEES_Node.h"

int CEES_Node::K;
vector <double> CEES_Node::H; 
vector <double> CEES_Node::T; 
double CEES_Node::pee; 
int CEES_Node::B; 
int CEES_Node::N; 
int CEES_Node::dataDim; 
int CEES_Node::depositFreq; 
CModel * CEES_Node::ultimate_target;

CEES_Node::CEES_Node(int iEnergyLevel, CTransitionModel *transition, CEES_Node *higher)
{
	id = iEnergyLevel; 
	nSamplesGenerated = 0; 

	proposal = transition;
	higher_level = higher; 

	x_current = new double [dataDim]; 
	x_new = new double [dataDim];

	target = new CBoundedModel(H[id], T[id], ultimate_target); 
}

CEES_Node::~CEES_Node()
{
	if (sizeof(x_current))
		delete [] x_current; 

	if (sizeof(x_new))
		delete [] x_new; 

	if (target != NULL)
		delete target;

	if (proposal != NULL)
		delete proposal; 
}

int CEES_Node::BinID(double e) const
{
	int energy_index = GetRingIndex(e); 
	return id*K+energy_index; 
}

int CEES_Node::BinID(int energy_index) const
{
	return id*K+energy_index; 
}

bool CEES_Node::BurnInDone()
{
	if (nSamplesGenerated >= B)
		return true; 
	return false;
}

void CEES_Node::Initialize(CModel * model, const gsl_rng *r)
{
	model->draw(x_current, dataDim, r); 
	energy_current = OriginalEnergy(x_current, dataDim); 
	ring_index_current = GetRingIndex(energy_current); 
	nSamplesGenerated ++;
}

void CEES_Node::draw(const gsl_rng *r)
{
	int bin_id;
	if (higher_level != NULL)
		bin_id = higher_level->BinID(ring_index_current); // energy ring at the next level whose range of energy is the same that of x_current	
	else 
		bin_id = -1; 

	bool new_sample_flag; 
	sdsm_data_item *item_bin; 
	if (bin_id <0 || sdsm_db_get_item_count(bin_id) == 0) // MH draw
	{
		target->draw(proposal, x_new, dataDim, x_current, r, new_sample_flag); 
		if (new_sample_flag)
		{
			memcpy(x_current, x_new, sizeof(x_new)); 
			energy_current = OriginalEnergy(x_new, dataDim); 
			ring_index_current = GetRingIndex(energy_current); 
		}
	}
	else	// equi-energy draw
	{
		double uniform_draw = gsl_rng_uniform(r); 
		if (uniform_draw <= pee)
		{
			item_bin = sdsm_get(bin_id); // randomly pick a data point from bin_id; 	
			memcpy(x_new, item_bin->data, sizeof(item_bin->data));
			double ratio=ProbabilityRatio(x_new, x_current, dataDim); 
			ratio = ratio * higher_level->ProbabilityRatio(x_current, x_new, dataDim);  
			double another_uniform_draw = gsl_rng_uniform(r); 
			if (another_uniform_draw <= ratio)
				memcpy(x_current, x_new, sizeof(x_new));
			sdsm_release(item_bin);
		} 
		else 
		{
			target->draw(proposal, x_new, dataDim, x_current, r, new_sample_flag); 
			if (new_sample_flag)
			{
				memcpy(x_current, x_new, sizeof(x_new)); 
				energy_current = OriginalEnergy(x_new, dataDim); 
				ring_index_current = GetRingIndex(energy_current); 
			}
		}
	}

	if (BurnInDone() && (nSamplesGenerated-B)%depositFreq==0)
	{
		bin_id = BinID(ring_index_current); 
		if (sdsm_put(bin_id, x_current, sizeof(x_current), 1) < 0)
			cout << "Error in sdsm_put.\n"; 
	}
	nSamplesGenerated ++;
}

double CEES_Node::ProbabilityRatio(const double *x, const double *y, int dim)
{
	return target->probability(x, dim)/target->probability(y, dim); 
}

void CEES_Node::SetDataDimension(int d)
{
	dataDim = d; 
}

int CEES_Node::GetDataDimension()  
{
	return dataDim;
}

void CEES_Node::SetEnergyLevelNumber(int i)
{
	K = i; 
}

int CEES_Node::GetEnergyLevelNumber() 
{
	return K; 
}

void CEES_Node::SetEquiEnergyJumpProb(double p)
{
	pee = p;
}

int CEES_Node::GetEquiEnergyJumpProb() 
{
	return pee; 
}

void CEES_Node::SetBurnInPeriod(int b)
{
	B = b; 
}

int CEES_Node::GetBurnInPeriod() 
{
	return B; 
}

void CEES_Node::SetPeriodBuildInitialRing(int n)
{
	N = n; 
}

int CEES_Node::GetPeriodBuildInitialRing() 
{
	return N;
}

void CEES_Node::SetDepositFreq(int f)
{
	depositFreq = f;
}

int CEES_Node::GetDepositFreq() 
{
	return depositFreq;
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
                if (Z[2*i]>0 && abs(Z[2*i+1]) <= EPSILON)
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

int CEES_Node::GetRingIndex(double e) const
{
	for (int j=1; j<K; j++)
	{
		if (e < H[j])
			return j-1; 
	}	
	return K-1;
}
