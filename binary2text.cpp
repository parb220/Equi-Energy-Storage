#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include "CSampleIDWeight.h"

using namespace std; 

bool GetNumberFilesPerBin(string, vector <int> &); 

int main(int argc, char* argv[])
{
	if (argc < 5)
	{
		cout << argv[0] << " data_dimension input_filename(binary) output_filename(text) max_record_number"; 
		exit(-1); 
	}
	CSampleIDWeight::SetDataDimension(atoi(argv[1])); 
	int max_record_number = atoi(argv[4]); 
	vector <CSampleIDWeight> sample(max_record_number);
	fstream iFile; 
	iFile.open(argv[2], ios::in|ios::binary); 
	if (!iFile)
	{
		cout << "Error in loading " << argv[2] << endl; 
		exit(-1); 
	}
	int counter = 0; 
	while(!read(iFile, &(sample[counter])).eof())
		counter ++ ;
	

	iFile.close(); 

	ofstream oFile; 
	oFile.open(argv[3]); 
	if (!oFile)
	{
		cout << "Error in writing " << argv[3] << endl; 
		exit(-1);
	}

	for (int i=0; i<counter; i++)
		oFile << sample[i]; 

	oFile.close(); 
}
