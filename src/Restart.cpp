
//------------------------------------------------------------------
//  **MEDYAN** - Simulation Package for the Mechanochemical
//               Dynamics of Active Networks, v4.0
//
//  Copyright (2015-2018)  Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the MEDYAN web page for more information:
//  http://www.medyan.org
//------------------------------------------------------------------

#include "Restart.h"

void Restart::readNetworkSetup() {
	//Flush the filestream pointer
	_inputFile.clear();
	//Go to the first entry in the file.
	_inputFile.seekg(0);

	string line;
	//get line
	while(getline(_inputFile, line)) {
		//Continue if commented
		if(line.find("#") != string::npos) { continue; }
		vector<string> lineVector = split<string>(line);
		if(line.size()>0) {
			if (lineVector[0] == "BEAD") {
				/*Get lines till the line is not empty*/
				getline(_inputFile, line);
				while (line.size() > 0) {
					cout << line << endl;
					//Split string based on space delimiter
					vector<string> lineVector = split<string>(line);
					_rBData.bidvec.push_back(atoi(((lineVector[0]).c_str())));
					//Copy coords & forces
					for (auto it = lineVector.begin() + 1;
					     it != lineVector.begin() + 5; it++)
						_rBData.coordvec.push_back(atof(((lineVector[0]).c_str())));
					for (auto it = lineVector.begin() + 5;
					     it != lineVector.begin() + 8; it++)
						_rBData.forceAuxvec.push_back(atof(((lineVector[0]).c_str())));
					getline(_inputFile, line);
				}
			}
		}
	}
}

void Restart::setupInitialNetwork() {
	for(unsigned int b=0;b<_rBData.bidvec.size();b++){
		auto bID = _rBData.bidvec[b];
		//Extract part of the vector.
		vector<floatingpoint> tempcoord(_rBData.coordvec.begin()+3*bID, _rBData.coordvec
		.begin()+3*bID+2);
		_subSystem->addTrackable<Bead>(tempcoord,nullptr,0);
		//Copy Forces
		for(unsigned int dim = 0; dim < 3; dim++)
			Bead::getDbData().forcesAux.data()[3*b+dim] = _rBData.coordvec.data()
					[3*bID+dim];
	}
	cout<<"Num beads "<<Bead::getBeads().size()<<endl;
}