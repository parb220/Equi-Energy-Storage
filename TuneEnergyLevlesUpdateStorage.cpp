#include "CEES_Node.h"
#include "CStorageHead.h"
#include "CBoundedModel.h"

using namespace std;

bool TuneEnergyLevels_UpdateStorage(CEES_Node *simulator, CStorageHead &storage, double c_factor, double mh_target_acc)
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
		CEES_Node::SetEnergyLevels_GeometricProgression(new_H0, new_HK_1);
		CEES_Node::SetTemperatures_EnergyLevels(CEES_Node::T[0], c_factor, true);
		CEES_Node::SetTargetAcceptanceRate(mh_target_acc); 
	
	// Re-adjust local target distribution and process samples that have been generated; 
		// storage.CreateTemporaryBin(); 
		for (int i=0; i<CEES_Node::K; i++)
		{
			simulator[i].AdjustLocalTarget(); 
			// simulator[i].AssignSamplesGeneratedSoFar(storage); 
			simulator[i].DisregardHistorySamples(storage); 
		}
		// storage.ClearTemporaryBin(); 
		return true; 
	}
	return false; 
} 
