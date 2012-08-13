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

void MHAdaptive::EstimateRegressionParameters(const vector < AccStep > &observation )
{
	double b_reserved = b; 
	double a_reserved = a; 
	
	double a_new = a; 
	double b_new = b; 
	double p, x;  
 	double denominator; 
	double sum_diff_p, sum_log_diff_p; 
	double A, B, C; 
	int count = 0; 
	bool denominator_continue = true; 
	while (denominator_continue && (count == 0 || fabs(a_new-a) > 1.0e-6 || fabs(b_new-b) > 1.0e-6))
	{
		a = a_new; 
		b = b_new; 
		A = B = C = sum_diff_p = sum_log_diff_p = 0.0; 
		for (int i=0; i<(int)(observation.size()); i++)
		{
			x = a+b*observation[i].step;  
			p = exp(x)/(1.0+exp(x)); // p = logit^(-1)(a+b*log(s))
			A += observation[i].nGenerated*p*(1.0-p); 
			B += observation[i].nGenerated*observation[i].step*p*(1.0-p); 
			C += observation[i].nGenerated*observation[i].step*observation[i].step*p*(1.0-p); 
			sum_diff_p += observation[i].nAccepted - observation[i].nGenerated*p; 
			sum_log_diff_p += observation[i].step * (observation[i].nAccepted - observation[i].nGenerated*p); 
		}

		denominator = A*C-B*B; 
		if (denominator > 1.0e-6)
		{
			a_new = a + (A*sum_diff_p - B*sum_log_diff_p)/denominator; 
			b_new = b - (B*sum_diff_p + C*sum_log_diff_p)/denominator;
			count ++;  
		}
		else
			denominator_continue = false; 
	}
	if ( denominator_continue )
	{
		a = a_new; 
		b = b_new;
		return; 
	}
	else
	{
		a = a_reserved; 
		b = b_reserved; 
	 	count = 0; 
		a_new = a; 
		while (count == 0 || fabs(a_new -a) > 1.0e-6)
		{
			a = a_new; 
			sum_diff_p = A = 0.0; 
	 		for (int i=0; i<(int)(observation.size()); i++)
                	{
                        	x = a+b*observation[i].step;
                        	p = exp(x)/(1.0+exp(x)); // p = logit^(-1)(a+b*log(s))
				sum_diff_p += observation[i].nAccepted - observation[i].nGenerated*p;
				A += observation[i].nGenerated*p*(1.0-p);
			}
			a_new = a + (sum_diff_p-(a+3.0)/25.0)/(A+1.0/25.0); 
			count ++;
		}
		a = a_new;
		return; 
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
