#include "CEES_Node.h"
#include "CStorageHead.h"
#include "CBoundedModel.h"
#include "CParameterPackage.h"

using namespace std;

bool TuneEnergyLevels_UpdateStorage(CEES_Node *simulator, CStorageHead &storage, CParameterPackage &parameter)
{
	double global_min = simulator[0].min_energy[0]; 
	double global_max = simulator[CEES_Node::K-1].max_energy[0]; 
	for (int i=1; i<CEES_Node::K; i++)
		global_min = global_min < simulator[i].min_energy[0] ? global_min : simulator[i].min_energy[0]; 
	double new_H0 = global_min < CEES_Node::H[0] ? global_min : CEES_Node::H[0];
	// double new_HK_1 = CEES_Node::max_energy[0] < 1.0e3 ?CEES_Node::max_energy[0]: 1.0e3;  
	double new_HK_1 = global_max < CEES_Node::H[CEES_Node::K-1] ? global_max : CEES_Node::H[CEES_Node::K-1]; 

	// Re-determine and adjust energy level and temperature levels
	if (new_H0 < CEES_Node::H[0] || new_HK_1 > CEES_Node::H[CEES_Node::K-1])
	{
		parameter.h0 = new_H0; 
		parameter.hk_1 = new_HK_1; 
		parameter.SetEnergyBound(); 
		parameter.SetTemperature(); 
		double *temp_buffer_float = new double[parameter.number_energy_level];
        	parameter.GetEnergyBound(temp_buffer_float, parameter.number_energy_level);
        	CEES_Node::SetEnergyLevels(temp_buffer_float, parameter.number_energy_level);
        	parameter.GetTemperature(temp_buffer_float, parameter.number_energy_level);
        	CEES_Node::SetTemperatures(temp_buffer_float, parameter.number_energy_level);
        	delete [] temp_buffer_float;

		// Re-adjust local target distribution and process samples that have been generated; 
		// storage.CreateTemporaryBin(); 
		for (int i=0; i<CEES_Node::K; i++)
			simulator[i].AdjustLocalTarget(); 
			// simulator[i].AssignSamplesGeneratedSoFar(storage); 
		storage.DisregardHistorySamples(); 
		return true; 
	}
	return false; 
} 
