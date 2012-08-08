#include <cfloat>
#include <cmath>
#include <algorithm>
#include "MHAdaptive.h"

using namespace std;
MHAdaptive::MHAdaptive(int _periodL, double _targetProb, double _a, double _b)
{
	// Target acceptance rate and its lower- and upper-bounds
	ResetTargetAcceptanceRate(_targetProb); 

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

void MHAdaptive::ResetTargetAcceptanceRate(double _targetProb)
{
        mid = _targetProb;
        log_mid = log(mid);
        lower_bound = exp(log_mid*2.0);
        upper_bound = exp(log_mid*0.5);
}

double MHAdaptive::GetStepSize(vector <AccStep> &observation)
{
	// sorting logit(acceptance rate) 
	stable_sort(observation.begin(), observation.end(), &comparator); 
	double small_step = observation[observation.size()/10].step;	// 10% percentile  
	double big_step = observation[observation.size()*9/10].step; 	// 90% percentile
	
	// position of lower_bound and upper_bound
	AccStep target_observation, lower_observation, upper_observation; 
	vector <AccStep>::iterator target_lower, target_upper, lower_lower, upper_upper; 
	bool continue_flag = true;
	while (continue_flag && fabs(upper_bound-lower_bound)>DBL_EPSILON && lower_bound>DBL_EPSILON && upper_bound<1.0-DBL_EPSILON)
	{ 
		target_observation.acc = log(mid)-log(1.0-mid);
		lower_observation.acc = log(lower_bound) - log(1.0-lower_bound); 
		upper_observation.acc = log(upper_bound) - log(1.0-upper_bound); 
		target_lower = std::lower_bound(observation.begin(), observation.end(), target_observation, &comparator); 
		target_upper = std::upper_bound(observation.begin(), observation.end(), target_observation, &comparator); 
		lower_lower = std::lower_bound(observation.begin(), observation.end(), lower_observation, &comparator); 
		upper_upper = std::upper_bound(observation.begin(), observation.end(), upper_observation, &comparator);
		
		if (lower_lower >= observation.end())
			ResetTargetAcceptanceRate(lower_bound); 
		else if (upper_upper <= observation.begin())
			ResetTargetAcceptanceRate(upper_bound); 
		else if (target_lower >= observation.end())
			ResetTargetAcceptanceRate(exp((log(lower_bound)+log(mid))/2)); 
		else if (target_upper <= observation.begin())
			ResetTargetAcceptanceRate(exp((log(upper_bound)+log(mid))/2)); 
		else 
			continue_flag = false; 
	}
	if (lower_bound<=DBL_EPSILON)	// lower_bound of acceptance is still too high
		return exp(big_step); 
	else if (upper_bound>=1.0-DBL_EPSILON) 	// upper bound of acceptance is still too low
		return exp(small_step); 
	else 
	{
		int begin_ = target_lower > observation.begin() ? target_lower-observation.begin() : 0; 
		int end_ = target_upper < observation.end() ? target_upper-observation.begin(): (int)(observation.size()-1); 
		double average =0; 
		for (int i=begin_; i<=end_; i++)
			average += observation[i].step; 
		average = average/(end_-begin_+1); 
		return exp(average); 
	}
}

double MHAdaptive::GetStepSizeRegression(int flag) const
{
	double logit; 
	if (flag <0)
		logit = log(lower_bound) - log(1.0-lower_bound); 
	else if (flag > 0)
		logit = log(upper_bound) - log(1.0-upper_bound); 
	else 
		logit = log(mid) -log(1.0-mid); 
	return exp((logit-a)/b); 
}
