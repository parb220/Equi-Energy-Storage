#include "equi_energy_setup_constant.h"
#include "CEES_Node.h"
#include "CStorageHead.h"
#include "CBoundedModel.h"

using namespace std;

void TuneEnergyLevels_UpdateStorage(CEES_Node *simulator, CStorageHead &storage)
{
	CEES_Node::if_tune_energy_level = false; 
	double new_H0 = CEES_Node::min_energy; 

	// Re-determine and adjust energy level and temperature levels
	if (!CEES_Node::SetEnergyLevels_GeometricProgression(new_H0, HK_1))
        {
                cout << "Error in setting energy levels." << endl;
                exit(-1);
        }

	if (!CEES_Node::SetTemperatures_EnergyLevels(T0, TK_1, C) )
        {
                cout << "Error in setting temperature levels." << endl;
                exit(-1);
        }
	
	// Re-adjust local target distribution and process samples that have been generated; 
	for (int i=0; i<CEES_Node::K; i++)
	{
		simulator[i].AdjustLocalTarget(); 
		simulator[i].AssignSamplesGeneratedSoFar(storage); 
	}
	storage.ClearTemporaryBin(); 
} 
