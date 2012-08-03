#include <cmath>
#include "MHAdaptive.h"

MHAdaptive::MHAdaptive(int _periodL, double _targetProb, double _a, double _b)
{
	// Target acceptance rate and its lower- and upper-bounds
	mid = _targetProb; 
	log_mid = log(mid); 
	lower_bound = exp(log_mid*5.0); 
	upper_bound = exp(log_mid*0.2); 

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
		new_scale = (previous_ratio > upper_bound ? 5.0 : log_mid/log(previous_ratio))*high_scale; 
	} 
	else if (high_jump_ratio < 0.0) // need to decrease scale
	{
		best_scale = scale; 
		new_scale = (previous_ratio < lower_bound ? 0.2 : log_mid/log(previous_ratio))*low_scale; 
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

void MHAdaptive::UpdateRegressionParameters(int nGenerated, int nAccepted, double current_scale, int iteration)
{
	double p; 
	// logit^(-1)(x) = exp(x)/(1+exp(x))
	double x = a+b*log(current_scale); 
	p = exp(x)/(1.0+exp(x)); // p = logit^(-1)(a+b*log(s))
	A += nGenerated*p*(1.0-p); 
	B += nGenerated*log(current_scale)*p*(1.0-p); 
	C += nGenerated*log(current_scale)*log(current_scale)*p*(1.0-p); 
	sum_diff_p += (nAccepted - nGenerated*p); 
	sum_log_diff_p += log(current_scale) * (nAccepted - nGenerated*p); 

	if (iteration>=1)	// need at least two periods to estimate a and b
	{
		a += (A*sum_diff_p - B*sum_log_diff_p)/(A*C-B*B); 
		b -= (B*sum_diff_p + C*sum_log_diff_p)/(A*C-B*B);
	}
	else 
	{
		a = a/2.0; 
		b = b/2.0;
	} 
}

double MHAdaptive::GetStepSizeRegression() const
{
	double logit = log(mid) -log(1.0-mid); 
	return exp((logit-a)/b); 
}
