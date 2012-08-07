#include <cfloat>
#include <cmath>
#include <algorithm>
#include "MHAdaptive.h"

using namespace std;
MHAdaptive::MHAdaptive(int _periodL, double _targetProb, double _a, double _b)
{
	// Target acceptance rate and its lower- and upper-bounds
	mid = _targetProb; 
	log_mid = log(mid); 
	lower_bound = exp(log_mid*2.0); 
	upper_bound = exp(log_mid*0.5); 

	period = _periodL; 
        best_scale = scale = 1.0; 
	previous_ratio = 0.0; 
        low_scale = high_scale = -1.0; 
        low_jump_ratio = high_jump_ratio = -1.0; 

	// used for regression
	a = _a; 
	b = _b; 
	A = 0.0; 
	B = 0.0; 
	C = 0.0; 
	sum_diff_p =0; 
	sum_log_diff_p = 0;
}

void MHAdaptive::UpdateScale(int nGenerated, int nAccepted)
{
	double new_scale, diff; 
	previous_ratio = (double)(nAccepted)/(double)(nGenerated); 
	if (previous_ratio < mid)
	{
		low_scale = scale; 
		low_jump_ratio = previous_ratio; 
	}	
	else 
	{
		high_scale = scale; 
		high_jump_ratio = previous_ratio; 
	}

	if (low_jump_ratio < 0.0)	// need to increase scale
	{
		best_scale = scale; 
		new_scale = (previous_ratio > upper_bound ? 2.0 : log_mid/log(previous_ratio))*high_scale; 
	} 
	else if (high_jump_ratio < 0.0) // need to decrease scale
	{
		best_scale = scale; 
		new_scale = (previous_ratio < lower_bound ? 0.5 : log_mid/log(previous_ratio))*low_scale; 
	}
	else 
	{
		diff = high_jump_ratio - low_jump_ratio; 
		new_scale = best_scale =  diff > 1.0e-6 ? ((mid-low_jump_ratio)*low_scale + (high_jump_ratio-mid)*high_scale)/diff : 0.5*(low_scale+high_scale); 
		period *= 2; 
		low_jump_ratio = high_jump_ratio = -1.0; 
	}
	scale = new_scale;	
}

void MHAdaptive::EstimateRegressionParameters(vector <AccStep> &observation)
{
	double p, x;  
	double denominator; 
	// logit^(-1)(x) = exp(x)/(1+exp(x))
	for (int i=0; i<(int)(observation.size()); i++)
	{
		x = a+b*observation[i].step; 
		p = exp(x)/(1.0+exp(x)); // p = logit^(-1)(a+b*log(s))
		A += observation[i].nGenerated*p*(1.0-p); 
		B += observation[i].nGenerated*observation[i].step*p*(1.0-p); 
		C += observation[i].nGenerated*observation[i].step*observation[i].step*p*(1.0-p); 
		sum_diff_p += (observation[i].nAccepted - observation[i].nGenerated*p); 
		sum_log_diff_p += observation[i].step * (observation[i].nAccepted - observation[i].nGenerated*p); 

		denominator = A*C-B*B; 
		if (abs(denominator) > DBL_EPSILON)	// need at least two periods to estimate a and b
		{
			a += (A*sum_diff_p - B*sum_log_diff_p)/denominator; 
			b -= (B*sum_diff_p + C*sum_log_diff_p)/denominator;
		}
	} 
}

bool comparator (const AccStep &l, const AccStep &r) { return l.acc < r.acc; }

double MHAdaptive::GetStepSize(vector <AccStep> &observation)
{
	// sorting acceptance rate and its company step_size
	stable_sort(observation.begin(), observation.end(), &comparator); 
	double small_step = observation[observation.size()/4].step;	// 25% percentile  
	double big_step = observation[observation.size()*3/4].step; 	// 75% percentile
	
	// position of lower_bound and upper_bound
	AccStep target_observation;
	target_observation.acc = log(mid)-log(1.0-mid);
	//AccStep lower_target_observation, upper_target_observation;  
	//lower_target_observation.acc = lower_bound; 
	//upper_target_observation.acc = upper_bound; 
	vector <AccStep>::iterator target_lower = std::lower_bound(observation.begin(), observation.end(), target_observation, &comparator); 
	vector <AccStep>::iterator target_upper = std::upper_bound(observation.begin(), observation.end(), target_observation, &comparator); 
	if (target_lower >= observation.end())
		return big_step;  
	else if (target_upper <= observation.begin())
		return small_step;  
	else 
	{
		EstimateRegressionParameters(observation); 
		return GetStepSizeRegression(); 
		// double average=0; 
		// for (int i=target_lower-observation.begin(); i<=target_upper-observation.begin(); i++)
		//	average += observation[i].step; 
		// average /= (target_upper-target_lower+1); 
		// return exp(average); 
	}
}

double MHAdaptive::GetStepSizeRegression() const
{
	double logit = log(mid) -log(1.0-mid); 
	return exp((logit-a)/b); 
}
