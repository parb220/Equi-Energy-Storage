#include <cfloat>
#include <cmath>
#include <algorithm>
#include "MHAdaptive.h"

using namespace std;
MHAdaptive::MHAdaptive(double _targetProb, double _scale, double _a, double _b)
{
	// Target acceptance rate and its lower- and upper-bounds
	ResetTargetAcceptanceRate(_targetProb); 

        best_scale = scale = _scale;  
	previous_ratio = 0.0; 
        low_scale = high_scale = -1.0; 
        low_jump_ratio = high_jump_ratio = -1.0; 

	// used for regression
	a = _a; 
	b = _b; 
}

bool MHAdaptive::UpdateScale(int nGenerated, int nAccepted)
{
	old_scale = scale; 
	double new_scale, diff; 
	previous_ratio = (double)(nAccepted)/(double)(nGenerated); 
	if (previous_ratio < mid)	// need to decrease scale
	{
		low_scale = scale; 
		low_jump_ratio = previous_ratio; 
	}	
	else 	// need to increase scale
	{
		high_scale = scale; 
		high_jump_ratio = previous_ratio; 
	}

	if (low_jump_ratio < 0.0)	// when previous_ratio >= mid
	{
		best_scale = scale; 
		new_scale = (previous_ratio > upper_bound ? 2.0 : log_mid/log(previous_ratio))*high_scale; 
	} 
	else if (high_jump_ratio < 0.0) // when previous_rati < mid
	{
		best_scale = scale; 
		new_scale = (previous_ratio < lower_bound ? 0.5 : log_mid/log(previous_ratio))*low_scale; 
	}
	else 
	{
		diff = high_jump_ratio - low_jump_ratio; 
		new_scale = best_scale =  diff > 1.0e-6 ? ((mid-low_jump_ratio)*low_scale + (high_jump_ratio-mid)*high_scale)/diff : 0.5*(low_scale+high_scale); 
		low_jump_ratio = high_jump_ratio = -1.0; 
	}
	scale = new_scale;
	
	if (fabs(scale - old_scale) > 1.0e-6)
		return true; 
	else 
		return false; 	
}

bool comparator (const AccStep &l, const AccStep &r) { return l.acc < r.acc; }

void MHAdaptive::ResetTargetAcceptanceRate(double _targetProb)
{
        mid = _targetProb;
        log_mid = log(mid);
        lower_bound = exp(log_mid*2.0);
        upper_bound = exp(log_mid*0.5);
}
