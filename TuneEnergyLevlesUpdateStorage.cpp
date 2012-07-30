#include "equi_energy_setup_constant.h"
#include "CEES_Node.h"
#include "CStorageHead.h"
#include "CBoundedModel.h"

using namespace std;

void TuneEnergyLevels_UpdateStorage(CEES_Node *simulator, CStorageHead &storage)
{
	double new_H0_average = 0; 
	for (int i=0; i<(int)(CEES_Node::min_energy.size()); i++)
		new_H0_average += CEES_Node::min_energy[i]; 
	new_H0_average = new_H0_average/(int)(CEES_Node::min_energy.size()); 
	double new_HK_1_average = 0;
	for (int i=0; i<(int)(CEES_Node::max_energy.size()); i++)
		new_HK_1_average += CEES_Node::max_energy[i]; 
	new_HK_1_average=new_HK_1_average/(int)(CEES_Node::max_energy.size()); 

	// Re-determine and adjust energy level and temperature levels
	if (new_H0_average < CEES_Node::H[0])
		CEES_Node::SetEnergyLevels_GeometricProgression(new_H0_average, HK_1);

	// if (!CEES_Node::SetTemperatures_EnergyLevels(T0, TK_1, C) )
	double new_TK_1 = new_HK_1_average * 100; 
	CEES_Node::SetTemperatures_EnergyLevels(T0, new_TK_1);
	
	// Re-adjust local target distribution and process samples that have been generated; 
	storage.CreateTemporaryBin(); 
	for (int i=0; i<CEES_Node::K; i++)
	{
		simulator[i].AdjustLocalTarget(); 
		simulator[i].AssignSamplesGeneratedSoFar(storage); 
	}
	storage.ClearTemporaryBin(); 
} 
