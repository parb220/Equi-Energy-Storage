#include <string>
#include <fstream>
#include <gsl/gsl_poly.h>
#include <cmath>
#include "CParameterPackage.h"
#include "CUniformModel.h"
using namespace std; 

CParameterPackage::CParameterPackage()
{
}

CParameterPackage::~CParameterPackage()
{
}		

bool CParameterPackage::SaveParameterToFile(string file_name) const
{
	fstream oFile(file_name.c_str(), ios::out | ios::binary); 
	if (!oFile)
		return false; 
	oFile.write((char*)(&number_cluster_node), sizeof(int)); 
	oFile.write((char*)(&run_id), sizeof(int)); 
	oFile.write((char*)(&get_marker), sizeof(int)); 
	oFile.write((char*)(&put_marker), sizeof(int)); 
	oFile.write((char*)(&number_energy_level), sizeof(int)); 
	oFile.write((char*)(&data_dimension), sizeof(int)); 
	oFile.write((char*)(&number_bins), sizeof(int)); 
	oFile.write((char*)(&pee), sizeof(double)); 
	oFile.write((char*)(&h0), sizeof(double)); 
	oFile.write((char*)(&hk_1), sizeof(double)); 
	oFile.write((char*)(&energy_tracking_number), sizeof(int)); 
	oFile.write((char*)(&t0), sizeof(double)); 
	oFile.write((char*)(&c_factor), sizeof(double));
	oFile.write((char*)(&mh_target_acc), sizeof(double)); 
	oFile.write((char*)(&number_block), sizeof(int));
	
	for (int i=0; i<number_block; i++)
		oFile.write((char*)(&(block_size[i])), sizeof(int));  
	
	oFile.write((char*)(&initial_sigma), sizeof(double)); 
	oFile.write((char*)(&uniform_lb), sizeof(double)); 
	oFile.write((char*)(&uniform_ub), sizeof(double)); 
	oFile.write((char*)(&burn_in_period), sizeof(int)); 
	oFile.write((char*)(&multiple_try_mh), sizeof(int)); 
	oFile.write((char*)(&mh_tracking_length), sizeof(int)); 
	oFile.write((char*)(&mh_stepsize_tuning_max_time), sizeof(int)); 
	oFile.write((char*)(&energy_level_tracking_window_length), sizeof(int)); 
	oFile.write((char*)(&energy_level_tuning_max_time), sizeof(int)); 
	oFile.write((char*)(&deposit_frequency), sizeof(int)); 
	oFile.write((char*)(&simulation_length), sizeof(int)); 

	for (int i=0; i<number_energy_level; i++)
		for (int j=0; j<data_dimension; j++)
			oFile.write((char*)(&(scale[i][j])), sizeof(double)); 
	
	// for (int i=0; i<number_bins; i++)
	//	oFile.write((char*)(&(number_samples_generated_by_far[i])), sizeof(int));
	
	for (int i=0; i<number_bins; i++)
		oFile.write((char*)(&(number_files_fetch[i])), sizeof(int)); 

	for (int i=0; i<number_energy_level; i++)
		oFile.write((char*)(&(h[i])), sizeof(double));

	for (int i=0; i<number_energy_level; i++)
	 	oFile.write((char*)(&(t[i])), sizeof(double)); 

	// for (int i=0; i<number_energy_level; i++)
	//	oFile.write((char*)(&(energy_index_current[i])), sizeof(int)); 
	oFile.close(); 
	return true; 
}

bool CParameterPackage::SaveCurrentStateToFile(string file_name) const
{
	fstream oFile(file_name.c_str(), ios::out | ios::binary);
        if (!oFile)
                return false;

	for (int i=0; i<number_energy_level; i++)
		write(oFile, &(x_current[i])); 
	oFile.close(); 
	return true;
}

bool CParameterPackage::LoadParameterFromFile(string file_name)
{
        fstream iFile(file_name.c_str(), ios::in | ios::binary);
        if (!iFile)
                return false;
	iFile.read((char*)(&number_cluster_node), sizeof(int)); 
        iFile.read((char*)(&run_id), sizeof(int));
        iFile.read((char*)(&get_marker), sizeof(int));
        iFile.read((char*)(&put_marker), sizeof(int));
        iFile.read((char*)(&number_energy_level), sizeof(int));
        iFile.read((char*)(&data_dimension), sizeof(int));
        iFile.read((char*)(&number_bins), sizeof(int));
        iFile.read((char*)(&pee), sizeof(double));
        iFile.read((char*)(&h0), sizeof(double));
        iFile.read((char*)(&hk_1), sizeof(double));
        iFile.read((char*)(&energy_tracking_number), sizeof(int));
        iFile.read((char*)(&t0), sizeof(double));
        iFile.read((char*)(&c_factor), sizeof(double));
        iFile.read((char*)(&mh_target_acc), sizeof(double));
	iFile.read((char*)(&number_block), sizeof(int)); 
	
	block_size.resize(number_block); 
	for (int i=0; i<number_block; i++)
		iFile.read((char*)(&(block_size[i])), sizeof(int));
        
	iFile.read((char*)(&initial_sigma), sizeof(double));
        iFile.read((char*)(&uniform_lb), sizeof(double));
        iFile.read((char*)(&uniform_ub), sizeof(double));
        iFile.read((char*)(&burn_in_period), sizeof(int));
        iFile.read((char*)(&multiple_try_mh), sizeof(int));
        iFile.read((char*)(&mh_tracking_length), sizeof(int));
        iFile.read((char*)(&mh_stepsize_tuning_max_time), sizeof(int));
        iFile.read((char*)(&energy_level_tracking_window_length), sizeof(int));
        iFile.read((char*)(&energy_level_tuning_max_time), sizeof(int));
        iFile.read((char*)(&deposit_frequency), sizeof(int));
        iFile.read((char*)(&simulation_length), sizeof(int));
	
	scale.resize(number_energy_level); 
	for (int i=0; i<number_energy_level; i++)
		scale[i].resize(data_dimension); 
	for (int i=0; i<number_energy_level; i++)
		for (int j=0; j<data_dimension; j++)
			iFile.read((char*)(&(scale[i][j])), sizeof(double));

	// number_samples_generated_by_far.resize(number_bins);
	// for (int i=0; i<number_bins; i++)
        //	iFile.read((char*)(&(number_samples_generated_by_far[i])), sizeof(int));

	number_files_fetch.resize(number_bins); 
	for (int i=0; i<number_bins; i++)
		iFile.read((char*)(&(number_files_fetch[i])), sizeof(int)); 

	h.resize(number_energy_level);
	for (int i=0; i<number_energy_level; i++) 
		iFile.read((char*)(&(h[i])), sizeof(double));
        
	t.resize(number_energy_level);
	for (int i=0; i<number_energy_level; i++) 
		iFile.read((char*)(&(t[i])), sizeof(double));

	iFile.close(); 
	return true; 
}

bool CParameterPackage::LoadCurrentStateFromFile(string file_name)
{
	x_current.resize(number_energy_level); 
	for (int i=0; i<number_energy_level; i++)
		x_current[i].SetDataDimension(data_dimension);
	fstream iFile(file_name.c_str(), ios::in | ios::binary);
        if (!iFile)
                return false;
	for (int i=0; i<number_energy_level; i++)
		read(iFile, &(x_current[i])); 	
	iFile.close(); 
	return true; 
}


bool CParameterPackage::WriteSummaryFile(string file_name) const
{
	ofstream oFile; 
	oFile.open(file_name.c_str(), ios::out); 
	if (!oFile) 
		return false; 
	oFile << "Burn In:\t" << burn_in_period << endl; 
	oFile << "Tune Energy Level Window Length:\t" << energy_level_tracking_window_length << endl; 
	oFile << "Tune Energy Level Number:\t" << energy_level_tuning_max_time << endl; 
	oFile << "Deposit Frequency:\t" << deposit_frequency << endl; 
	oFile << "MH Target Probability:\t" << mh_target_acc << endl; 
	oFile << "MH Initial Window Length:\t" << mh_tracking_length << endl; 
	oFile << "MH Window Number:\t" << mh_stepsize_tuning_max_time << endl; 
	oFile << "Get Marker:\t" << get_marker << endl; 
	oFile << "Put Marker:\t" << put_marker << endl; 
	oFile << "Number of Bins:\t" << number_bins << endl; 
	for (int i=0; i<number_bins; i++)
		// oFile << i << ":\t" << number_samples_generated_by_far[i] << "\t" << number_files_by_far[i] << endl; 
		oFile << i << ":\t" << number_files_fetch[i]*put_marker << "\t" << number_files_fetch[i] << endl; 
	oFile << "Data Dimension:\t" << data_dimension << endl; 
	oFile << "Number of Energy Levels:\t" << number_energy_level << endl; 
	oFile << "Energy Thresholds:"; 
	for (int i=0; i<number_energy_level; i++)
		oFile << "\t" << h[i]; 
	oFile << endl; 
	oFile << "Temperatures:"; 
	for (int i=0; i<number_energy_level; i++)
		oFile << "\t" << t[i]; 
	oFile << endl; 
	oFile << "Prob Equi-Jump:\t" << pee << endl;  
	for (int i=0; i<number_energy_level; i++)
	{
		oFile << "Step size " << i << ":"; 
		for (int j=0; j<data_dimension; j++)
			oFile << "\t" << scale[i][j]; 
		oFile << endl; 
	}
	oFile.close(); 
	return true; 
}

bool CParameterPackage::SetEnergyBound()
{
        /*
 *  *     H[i] = H[i-1]+gamma^i
 *   *     gamma is determined by solving a polynomial equation 
 *    *     gamma+gamma^2+...+gamma^{K-1} = H[K-1]-H[0]; 
 *     *     */

        double *coefficients = new double [number_energy_level];
        coefficients[0] = h0-hk_1;
        for (int i=1; i<number_energy_level; i++)
                coefficients[i]=1.0;
        double *Z = new double [(number_energy_level-1)*2];

        gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(number_energy_level);
        gsl_poly_complex_solve(coefficients, number_energy_level, w, Z);

        double gamma;
        bool continue_flag = true;
        for (int i=0; i<number_energy_level-1 && continue_flag; i++)
        {
                if (Z[2*i]>0 && abs(Z[2*i+1]) <= 1.0e-6)
                {
                        gamma = Z[2*i];
                        continue_flag = false;
                }
        }
        delete [] Z;
        delete [] coefficients;
        if (continue_flag)
                return false;

	h.resize(number_energy_level);  
	h[0] = h0; 
	h[number_energy_level-1] = hk_1; 
	for (int i=1; i<number_energy_level-1; i++)
		h[i] = h[i-1]+pow(gamma, i); 
	return true; 
}

bool CParameterPackage::SetTemperature()
{
	t.resize(number_energy_level); 
	t[0] = t0; 
	for (int i=1; i<number_energy_level; i++)
		t[i] = t[i-1]+(h[i]-h[i-1])/c_factor; 
	return true;
}

bool CParameterPackage::SetMHProposalScale()
{
	scale.resize(number_energy_level); 
	for (int i=0; i<number_energy_level; i++)
	{
		scale[i].resize(data_dimension); 
		for (int j=0; j<data_dimension; j++)
			scale[i][j] = initial_sigma*t[i]; 
	}
	return true; 
}

void CParameterPackage::SetCurrentState(const gsl_rng *r)
{
	double *lB = new double [data_dimension];
        double *uB = new double [data_dimension];
        for (int i=0; i<data_dimension; i++)
        {
                lB[i] = uniform_lb;
                uB[i] = uniform_ub;
        }
	CUniformModel initial_model(data_dimension, lB, uB); 
	bool if_new_sample; 
	x_current.resize(number_energy_level); 
	for (int i=0; i<number_energy_level; i++)
		initial_model.draw(x_current[i], if_new_sample, r); 
	delete [] lB; 
	delete [] uB; 
}

void CParameterPackage::SetBlock(int *_block_size)
{
	block_size.resize(number_block); 
	if (number_block == data_dimension)
	{
		for (int i=0; i<number_block; i++)
			block_size[i] = 1; 
	}
	else if (number_block == 1)
		block_size[0] = data_dimension; 
	else 
	{
		for (int i=0; i<number_block; i++)
			block_size[i] = _block_size[i]; 
	}
}

void CParameterPackage::GetMHProposalScale (int _id, double *_buffer, int _buffer_size) const 
{
	for (int j=0; j<data_dimension; j++)
		_buffer[j] = scale[_id][j]; 
}

void CParameterPackage::GetEnergyBound(double *_buffer, int _buffer_size) const
{
	for (int i=0; i<number_energy_level; i++)
		_buffer[i] = h[i]; 
}

void CParameterPackage::GetTemperature(double *_buffer, int _buffer_size) const
{
	for (int i=0; i<number_energy_level; i++)
		_buffer[i] = t[i]; 
}

void CParameterPackage::GetBlockSize(int *_buffer, int _buffer_size) const
{
	for (int i=0; i<number_block; i++)
		_buffer[i] = block_size[i]; 
}

