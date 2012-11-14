#include <string>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <vector>
#include "CSampleIDWeight.h"

using namespace std; 

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		cout << argv[0] << " input_filename(binary) output_filename(text)"; 
		exit(-1); 
	} 
	CSampleIDWeight one_sample; 
	vector <CSampleIDWeight> sample;
	fstream iFile; 
	iFile.open(argv[1], ios::in|ios::binary); 
	if (!iFile)
	{
		cout << "Error in loading " << argv[1] << endl; 
		exit(-1); 
	}
	while(!read(iFile, &(one_sample)).eof())
		sample.push_back(one_sample); 

	iFile.close(); 

	ofstream oFile; 
	oFile.open(argv[2]); 
	if (!oFile)
	{
		cout << "Error in writing " << argv[2] << endl; 
		exit(-1);
	}

	for (int i=0; i<(int)(sample.size()); i++)
		oFile << sample[i]; 

	oFile.close(); 
}
