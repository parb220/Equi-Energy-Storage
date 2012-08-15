#include "CEES_Node.h"
#include "CStorageHead.h"
#include "CBoundedModel.h"

using namespace std;

bool TuneEnergyLevels_UpdateStorage(CEES_Node *simulator, CStorageHead &storage, double c_factor, double mh_target_acc)
{
	double new_H0 = CEES_Node::min_energy[0] < CEES_Node::H[0] ? CEES_Node::min_energy[0] : CEES_Node::H[0];
	// double new_HK_1 = CEES_Node::max_energy[0] < 1.0e3 ?CEES_Node::max_energy[0]: 1.0e3;  
	double new_HK_1 = CEES_Node::max_energy[0] < CEES_Node::H[CEES_Node::K-1] ?CEES_Node::max_energy[0]: CEES_Node::H[CEES_Node::K-1]; 

	// Re-determine and adjust energy level and temperature levels
	if (new_H0 < CEES_Node::H[0] || new_HK_1 > CEES_Node::H[CEES_Node::K-1])
	{
		CEES_Node::SetEnergyLevels_GeometricProgression(new_H0, new_HK_1);
		CEES_Node::SetTemperatures_EnergyLevels(CEES_Node::T[0], c_factor, true);
		CEES_Node::SetTargetAcceptanceRate(mh_target_acc); 
	
	// Re-adjust local target distribution and process samples that have been generated; 
		storage.CreateTemporaryBin(); 
		for (int i=0; i<CEES_Node::K; i++)
		{
			simulator[i].AdjustLocalTarget(); 
			simulator[i].AssignSamplesGeneratedSoFar(storage); 
		}
		storage.ClearTemporaryBin(); 
		return true; 
	}
	return false; 
} 
