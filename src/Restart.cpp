
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
			//Get lines till the line is not empty

				getline(_inputFile, line);
				while (line.size() > 0) {
					//Split string based on space delimiter
					vector<string> lineVector = split<string>(line);
					//Bead ID
					_rBData.bsidvec.push_back(atoi(((lineVector[0]).c_str())));
					//Filament ID
					_rBData.filidvec.push_back(atoi(lineVector[1].c_str()));
					//Filament pos
					_rBData.filpos.push_back(atoi(lineVector[2].c_str()));
					//Copy coords & forces
					for (auto it = lineVector.begin() + 3;
					     it != lineVector.begin() + 6; it++)
						_rBData.coordvec.push_back(atof((*it).c_str()));
					for (auto it = lineVector.begin() + 6;
					     it != lineVector.begin() + 9; it++)
						_rBData.forceAuxvec.push_back(atof((*it).c_str()));

					getline(_inputFile, line);
				}
			}
			else if (lineVector[0] == "CYLINDER") {
			//Get lines till the line is not empty

				getline(_inputFile, line);
				while (line.size() > 0) {
					restartCylData _rcyldata;

					//Split string based on space delimiter
					vector<string> lineVector = split<string>(line);
					//Cylinder stable index
					_rcyldata.cylsid = atoi(((lineVector[0]).c_str()));
					//Filament ID
					_rcyldata.filid = atoi((lineVector[1]).c_str());
					_rcyldata.filtype = atoi((lineVector[2]).c_str());
					_rcyldata.filpos = atoi((lineVector[3]).c_str());
					//Filament pos
					//Bead stable indices
					_rcyldata.beadsidpairvec[0] = atoi((lineVector[4]).c_str());
					_rcyldata.beadsidpairvec[1] = atoi((lineVector[5]).c_str());
					//minus end or plus end?
					_rcyldata.endstatusvec[0] = atoi((lineVector[6]).c_str());
					_rcyldata.endstatusvec[1] = atoi((lineVector[7]).c_str());
					//minus/plus end type
					_rcyldata.endtypevec[0] = atoi((lineVector[8]).c_str());
					_rcyldata.endtypevec[1] = atoi((lineVector[9]).c_str());
					//endmonomerpos
					_rcyldata.endmonomerpos[0] = atoi((lineVector[10]).c_str());
					_rcyldata.endmonomerpos[1] = atoi((lineVector[11]).c_str());
					//totalmonomers
					_rcyldata.totalmonomers = atoi((lineVector[12]).c_str());
					//eqlen
					_rcyldata.eqlen = atof((lineVector[13]).c_str());
					//append to the vector
					_rCDatavec.push_back(_rcyldata);

					//get next cylinder data
					getline(_inputFile, line);
				}
			}
			else if (lineVector[0] == "FILAMENT") {
			//Get lines till the line is not empty

				getline(_inputFile, line);
				while (line.size() > 0) {
					restartFilData _rfildata;

					//Split string based on space delimiter
					vector<string> lineVector = split<string>(line);
					//Filament ID
					_rfildata.filid = atoi((lineVector[0]).c_str());
					//Filament Type
					_rfildata.filType = atoi((lineVector[1]).c_str());
					//Cylinder id vec
					for (auto it = lineVector.begin() + 2; it != lineVector.end(); it++) {
						_rfildata.cylsidvec.push_back(atoi((*it).c_str()));
					}
					//append to the vector
					_rFDatavec.push_back(_rfildata);

					//get next filament data
					getline(_inputFile, line);
				}
			}

			else if (lineVector[0] == "LINKER") {
			//Get lines till the line is not empty

				getline(_inputFile, line);
				while (line.size() > 0) {
					cout << line << endl;
				}
			}

			else if (lineVector[0] == "MOTOR") {
				//Get lines till the line is not empty

				getline(_inputFile, line);
				while (line.size() > 0) {
					restartMotorData _rmdata;
					//Split string based on space delimiter
					vector<string> lineVector = split<string>(line);
					//Motor ID
					_rmdata.motorid = atoi((lineVector[0]).c_str());
					//Cyl ID1
					_rmdata.cylid1 = atoi((lineVector[1]).c_str());
					//Cyl ID2
					_rmdata.cylid2 = atoi((lineVector[2]).c_str());
					//pos1
					_rmdata.pos1 = atof((lineVector[3]).c_str());
					//pos2
					_rmdata.pos2 = atof((lineVector[4]).c_str());
					//pos2
					_rmdata.eqlen = atof((lineVector[5]).c_str());

					_rMDatavec.push_back(_rmdata);
					//get next filament data
					getline(_inputFile, line);
				}
			}
		}
	}
}

void Restart::setupInitialNetwork() {

	//Step 1. Create dummy filaments
	//Step 2. Create Beads
	//Step 3. Create cylinders & set plus/minus ends where necessary
	//Step 4. Associate cylinders with each filament by adding it to the cylinder vector
	// in the appropriate order starting from minus end all the way to plus end cylinder.

	map<int, Filament*> filamentmap;//makes it easy to access Filament pointer from
	// Filament ID. Filaments do not follow stable index protocol and hence need not have
	// continuous ID values.

	for(auto &fil : _rFDatavec) {
		//Create dummy filament
		fil.filamentpointer = _subSystem->addTrackable<Filament>(_subSystem, fil.filType);
		//override ID
		fil.filamentpointer->overrideId(fil.filid);
		//add to map
		filamentmap[fil.filid] = fil.filamentpointer;
	}
	cout<<"Num filaments "<<Filament::getFilaments().size()<<endl;

	for(unsigned int b=0;b<_rBData.bsidvec.size();b++){
		auto bID = _rBData.bsidvec[b];
		auto filptr = filamentmap[_rBData.filidvec[b]];
		//Extract part of the vector.
		vector<floatingpoint> tempcoord(_rBData.coordvec.begin()+3*bID, _rBData.coordvec
		.begin()+3*bID+3);
		//Initialize beads
		_subSystem->addTrackable<Bead>(tempcoord, filptr, _rBData.filpos[b]);
		//Copy Forces
		for(unsigned int dim = 0; dim < 3; dim++)
			Bead::getDbData().forcesAux.data()[3*b+dim] = _rBData.coordvec.data()[3*bID+dim];
	}
	cout<<"Num beads "<<Bead::getBeads().size()<<endl;

	for(auto &cyl : _rCDatavec){
		auto b1 = Bead::getBeads()[cyl.beadsidpairvec[0]];
		auto b2 = Bead::getBeads()[cyl.beadsidpairvec[1]];
		auto filptr = filamentmap[cyl.filid];
		auto _filType = cyl.filtype;
		//initialize cylinder
		Cylinder* c0 = _subSystem->addTrackable<Cylinder> (filptr, b1, b2, _filType,
		                                                   cyl.filpos, false, false,
		                                                   true, cyl.eqlen);
		cyl.cylinderpointer = c0;
		//set minusend or plusend
		if(cyl.endstatusvec[0])
			c0->setPlusEnd(true);
		else if(cyl.endstatusvec[1])
			c0->setMinusEnd(true);
	}
	cout<<"Num cylinders "<<Cylinder::getCylinders().size()<<endl;

	for(auto fil : _rFDatavec) {
		vector<Cylinder*> cylvector;
		for(auto cylsid:fil.cylsidvec){
			cylvector.push_back(_rCDatavec[cylsid].cylinderpointer);
		}
		fil.filamentpointer->initializerestart(cylvector, _rCDatavec);
	}


}