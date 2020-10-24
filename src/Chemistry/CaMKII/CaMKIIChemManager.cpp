
#include "CaMKIIChemManager.h"

#include <cmath>

#include "ChemManager.h"

#include "ChemCallbacks.h"
#include "CompartmentGrid.h"
#include "ReactionTemplate.h"
#include "BindingManager.h"

#include "SysParams.h"
#include "MathFunctions.h"

using namespace mathfunc;

void CaMKIIChemManager::setupBindingSites(ChemistryData &_chemData, int filType) {
	if (_chemData.CaMKII_BINDING_INDEX[filType] != "") {
		auto it = find(_chemData.speciesBound[filType].begin(),
					   _chemData.speciesBound[filType].end(),
					   _chemData.CaMKII_BINDING_INDEX[filType]);

		if (it == _chemData.speciesBound[filType].end()) {

			cout << "The camkiier binding site listed is not a valid bound species. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		} else {
			SysParams::CParams.camkiierBindingBoundIndex[filType] = it - _chemData.speciesBound[filType].begin();
		}
	}

	if (_chemData.CaMKII_BUNDLING_INDEX[filType] != "") {
		auto it = find(_chemData.speciesBound[filType].begin(),
					   _chemData.speciesBound[filType].end(),
					   _chemData.CaMKII_BUNDLING_INDEX[filType]);

		if (it == _chemData.speciesBound[filType].end()) {

			cout << "The camkiier bundling site listed is not a valid bound species. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		} else {
			SysParams::CParams.camkiierBundlingBoundIndex[filType] = it - _chemData.speciesBound[filType].begin();
		}
	}

}

void CaMKIIChemManager::setupBindingSitesInitCylinders(int filType) {
	//for initialization of cylinders, for CaMKII, based on Aravind suggestion

	if (SysParams::CParams.brancherBoundIndex[filType] != SysParams::CParams.camkiierBindingBoundIndex[filType] &&
		SysParams::CParams.linkerBoundIndex[filType] != SysParams::CParams.camkiierBindingBoundIndex[filType] &&
		SysParams::CParams.motorBoundIndex[filType] != SysParams::CParams.camkiierBindingBoundIndex[filType])
		SysParams::CParams.bindingIndices[filType].push_back(SysParams::CParams.camkiierBindingBoundIndex[filType]);

	if (SysParams::CParams.brancherBoundIndex[filType] != SysParams::CParams.camkiierBundlingBoundIndex[filType] &&
		SysParams::CParams.linkerBoundIndex[filType] != SysParams::CParams.camkiierBundlingBoundIndex[filType] &&
		SysParams::CParams.motorBoundIndex[filType] != SysParams::CParams.camkiierBundlingBoundIndex[filType] &&
		SysParams::CParams.camkiierBindingBoundIndex[filType] != SysParams::CParams.camkiierBundlingBoundIndex[filType])
		SysParams::CParams.bindingIndices[filType].push_back(SysParams::CParams.camkiierBundlingBoundIndex[filType]);

}

void CaMKIIChemManager::initCMonomerCaMKII(ChemistryData &_chemData, CMonomer *m, short filamentType, Compartment *c,
										   int &bIndex) {
	for (auto &ca : _chemData.speciesCaMKIIer[filamentType]) {
		SpeciesCaMKIIer *sca =
				c->addSpeciesCaMKIIer(SpeciesNamesDB::genUniqueFilName(ca));
		m->_speciesBound[bIndex] = sca;
		bIndex++;
	}
	for (auto &ca : _chemData.speciesCaMKIICylinder[filamentType]) {
		SpeciesCaMKIICylinder *sca =
				c->addSpeciesCaMKIICylinder(SpeciesNamesDB::genUniqueFilName(ca));
		sca->getRSpecies().setUpperLimitForN(CAMKII_MAX_COORD_NUMBER_UPPER_LIMIT);
		m->_speciesBound[bIndex] = sca;
		bIndex++;
	}
}

void CaMKIIChemManager::genFilBindingReactionsCaMKII(SubSystem *_subSystem, ChemistryData &_chemData, int filType,
													 CompartmentGrid *grid,
													 Compartment *C,
													 int &managerIndex, int &camkiiIndex) {

	double rMax, rMin;

	for (auto &r: _chemData.camkiibindingReactions[filType]) {

		vector<Species *> reactantSpecies;
		vector<Species *> productSpecies;

		vector<string> reactants = get<0>(r);
		vector<string> products = get<1>(r);

		//Checks on number of reactants, products
		if (reactants.size() != CAMKIIBINDINGREACTANTS ||
			products.size() != CAMKIIBINDINGPRODUCTS) {
			cout << "Invalid CaMKII Binding reaction. Number of reactants or products does not match. Exiting." << endl;
			exit(EXIT_FAILURE);
		}
		string camkiierName;

		//FIRST PRODUCT SPECIES MUST BE CAMKIIER
		short camkiierInt;

		auto product = products[0];
		if (product.find("CAMKIIER") != string::npos) {

			//look up species, make sure in list
			string name = product.substr(0, product.find(":"));
			auto it = find(_chemData.speciesCaMKIIer[filType].begin(), _chemData.speciesCaMKIIer[filType].end(), name);

			if (it != _chemData.speciesCaMKIIer[filType].end()) {

				camkiierName = name;

				//get position of iterator
				camkiierInt = distance(_chemData.speciesCaMKIIer[filType].begin(), it);
			} else {
				cout <<
					 "A CaMKIIer species that was included in a reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			cout <<
				 "First product listed in a CaMKII binding reaction must be camkiier. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}

		//FIRST REACTANT MUST BE BULK OR DIFFUSING
		auto reactant = reactants[0];

		if (reactant.find("BULK") != string::npos) {

			//Look up species, make sure in list
			string name = reactant.substr(0, reactant.find(":"));

			auto it = find_if(_chemData.speciesBulk.begin(), _chemData.speciesBulk.end(),
							  [name](tuple<string, int, double, double, string> element) {
								  return get<0>(element) == name ? true : false;
							  });

			if (it == _chemData.speciesBulk.end()) {
				cout <<
					 "A bulk species that was included in a CaMKII binding reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
			reactantSpecies.push_back(grid->findSpeciesBulkByName(name));
		} else if (reactant.find("DIFFUSING") != string::npos) {

			//Look up species, make sure in list
			string name = reactant.substr(0, reactant.find(":"));

			auto it = find_if(_chemData.speciesDiffusing.begin(), _chemData.speciesDiffusing.end(),
							  [name](tuple<string, int, double, double, double, string, int> element) {
								  return get<0>(element) == name ? true : false;
							  });
			if (it == _chemData.speciesDiffusing.end()) {
				cout <<
					 "A diffusing species that was included in a CaMKII binding reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
			reactantSpecies.push_back(C->findSpeciesByName(name));
		} else {
			cout <<
				 "First species listed in a CaMKII binding reaction must be either bulk or diffusing. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}

		//SECOND REACTANT SPECIES MUST BE BOUND
		reactant = reactants[1];
		if (reactant.find("BOUND") != string::npos) {

			//look up species, make sure in list
			string name = reactant.substr(0, reactant.find(":"));
			auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
			int position = 0;

			if (it != _chemData.speciesBound[filType].end()) {

				//get position of iterator
				position = distance(_chemData.speciesBound[filType].begin(), it);

				if (position != SysParams::CParams.camkiierBindingBoundIndex[filType]) {
					cout <<
						 "Second species listed in a CaMKII binding reaction must be the corresponding camkiier empty site. Exiting."
						 << endl;
					exit(EXIT_FAILURE);
				}

				//find the species single binding, push "BindingTemplate" + the name
				string bename = SpeciesNamesDB::genBindingName("BindingTemplate" + camkiierName, name);

				reactantSpecies.push_back(C->findSpeciesByName(bename));
			} else {
				cout <<
					 "A bound species that was included in a CaMKII binding reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			cout <<
				 "Second species listed in a CaMKII binding reaction must be bound. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}

		//Create reaction
		float onRate = get<2>(r);
		float offRate = get<3>(r);
		auto temp = SysParams::CaMKIIUnbindingBareRate;
		if (temp.size() > 0)
			temp[camkiierInt] = offRate;
		else
			temp.push_back(offRate);
		SysParams::CaMKIIUnbindingBareRate = temp;
		//get nucleation zone
		string nzstr = get<4>(r);
		NucleationZoneType nucleationZone;
		if (nzstr == "ALL")
			nucleationZone = NucleationZoneType::ALL;
		else if (nzstr == "BOUNDARY")
			nucleationZone = NucleationZoneType::BOUNDARY;
		else if (nzstr == "TOPBOUNDARY")
			nucleationZone = NucleationZoneType::TOPBOUNDARY;
		else {
			cout << "Nucleation zone type specified in a CaMKII binding reaction not valid. Exiting." << endl;
			exit(EXIT_FAILURE);
		}
		double nucleationDist = get<5>(r);

		ReactionBase *rxn = new Reaction<2, 0>(reactantSpecies, onRate);
		rxn->setReactionType(ReactionType::CAMKIIBINDING);

		C->addInternalReaction(rxn);

		//create manager
		CaMKIIBindingManager *bManager = new CaMKIIBindingManager(rxn, C, camkiierInt, camkiierName, filType,
																  nucleationZone, nucleationDist);
		C->addFilamentBindingManager(bManager);

		bManager->setMIndex(managerIndex++);

		//attach callback
		CaMKIIBindingCallback camkiicallback(bManager, onRate, offRate, _subSystem);
		ConnectionBlock rcb(rxn->connect(camkiicallback, false));
	}


	for (auto &r: _chemData.camkiibundlingReactions[filType]) {

		vector<Species *> reactantSpecies;
		vector<Species *> productSpecies;

		vector<string> reactants = get<0>(r);
		vector<string> products = get<1>(r);

		//Checks on number of reactants, products
		if (reactants.size() != CAMKIIBUNDLINGREACTANTS ||
			products.size() != CAMKIIBUNDLINGPRODUCTS) {
			cout << "Invalid CaMKII bundling reaction. Number of reactants or products does not match. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}
		string camkiierName;

		//FIRST PRODUCT SPECIES MUST BE CAMKIIER
		short camkiierProductInt, camkiierReactantInt;

		auto product = products[0];
		if (product.find("CAMKIIER") != string::npos) {

			//look up species, make sure in list
			string name = product.substr(0, product.find(":"));
			auto it = find(_chemData.speciesCaMKIIer[filType].begin(), _chemData.speciesCaMKIIer[filType].end(), name);

			if (it != _chemData.speciesCaMKIIer[filType].end()) {

				camkiierName = name;

				//get position of iterator
				camkiierProductInt = distance(_chemData.speciesCaMKIIer[filType].begin(), it);
			} else {
				cout <<
					 "A camkiier species that was included in a CaMKII bundling reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			cout <<
				 "Fourth species listed in a CaMKII bundling reaction must be camkiier. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}

		//FIRST REACTANT MUST BE CAMKIIER
		auto reactant = reactants[0];
		if (reactant.find("CAMKIIER") != string::npos) {

			//Look up species, make sure in list
			string name = reactant.substr(0, reactant.find(":"));
			auto it = find(_chemData.speciesCaMKIIer[filType].begin(), _chemData.speciesCaMKIIer[filType].end(), name);

			if (it != _chemData.speciesCaMKIIer[filType].end()) {

				camkiierName = name;

				//get position of iterator
				camkiierReactantInt = distance(_chemData.speciesCaMKIIer[filType].begin(), it);

				// The product and reactant should be the same.
				if (camkiierProductInt != camkiierReactantInt) {
					cout <<
						 "First species listed in a CaMKII bundling reaction must be the same camkiier as the product. Exiting."
						 << endl;
					exit(EXIT_FAILURE);
				}
			} else {
				cout <<
					 "A camkiier species that was included in a CaMKII bundling reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			cout <<
				 "Third species listed in a CaMKII bundling reaction must be bound. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}

		//SECOND REACTANT SPECIES MUST BE BOUND
		reactant = reactants[1];
		if (reactant.find("BOUND") != string::npos) {

			//look up species, make sure in list
			string name = reactant.substr(0, reactant.find(":"));
			auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(), name);
			int position = 0;

			if (it != _chemData.speciesBound[filType].end()) {

				//get position of iterator
				position = distance(_chemData.speciesBound[filType].begin(), it);

				if (position != SysParams::CParams.camkiierBundlingBoundIndex[filType]) {
					cout <<
						 "Second species listed in a CaMKII bundling reaction must be the corresponding camkiier empty site. Exiting."
						 << endl;
					exit(EXIT_FAILURE);
				}

				//find the species single binding, push
				string bename = SpeciesNamesDB::genBindingName(camkiierName, name);

				reactantSpecies.push_back(C->findSpeciesByName(bename));
			} else {
				cout <<
					 "A bound species that was included in a CaMKII bundling reaction was not initialized. Exiting."
					 << endl;
				exit(EXIT_FAILURE);
			}
		} else {
			cout <<
				 "Second species listed in a CaMKII bundling reaction must be bound. Exiting."
				 << endl;
			exit(EXIT_FAILURE);
		}
		//Create reaction
		float onRate = get<2>(r);
		float offRate = get<3>(r);
		auto temp = SysParams::CaMKIIUnbundlingBareRate;
		if (temp.size() > 0)
			temp[camkiierProductInt] = offRate;
		else
			temp.push_back(offRate);
		SysParams::CaMKIIUnbundlingBareRate = temp;
		rMin = get<4>(r);
		rMax = get<5>(r);
		int maxCoordination = get<6>(r);

		ReactionBase *rxn = new Reaction<1, 0>(reactantSpecies, onRate);
		rxn->setReactionType(ReactionType::CAMKIIBUNDLING);
		C->addInternalReaction(rxn);

		//create manager
		CaMKIIBundlingManager *bManager = new CaMKIIBundlingManager(rxn, C, camkiierProductInt, camkiierName,
																	filType, rMax, rMin, maxCoordination);
		C->addFilamentBindingManager(bManager);

		bManager->setNLIndex(camkiiIndex++);
		bManager->setMIndex(managerIndex++);

		//attach callback
		CaMKIIBundlingCallback camkiicallback(bManager, onRate, offRate, _subSystem);
		ConnectionBlock rcb(rxn->connect(camkiicallback, false));
	}
}

void CaMKIIChemManager::initializeManagers(SubSystem *_subSystem, Compartment *C0, int filType) {

	for (auto &manager : C0->getFilamentBindingManagers()) {

		CaMKIIBundlingManager *bManager; //added for CaMKII

		if ((bManager = dynamic_cast<CaMKIIBundlingManager *>(manager.get()))) {

			auto nl =
					new CylinderCylinderNL(bManager->getRMax() + SysParams::Geometry().cylinderSize[filType],
										   max(bManager->getRMin() - SysParams::Geometry().cylinderSize[filType], 0.0),
										   true);


			CaMKIIBundlingManager::_neighborLists.push_back(nl);
			_subSystem->addNeighborList(nl);
		}
	}

}

void CaMKIIChemManager::genSpeciesCaMKII(ChemistryData &_chemData, Compartment& protoCompartment) {

	for(int filType = 0; filType < SysParams::Chemistry().numFilaments; filType++) {
		for (auto &sb : _chemData.speciesCaMKIIer[filType]) {

			//look at camkii binding reaction that is associated with this species
			for (auto &rb : _chemData.camkiibindingReactions[filType]) {

				auto reactants = get<0>(rb);
				auto products = get<1>(rb);

				auto sb_bound = products[0].substr(0, products[0].find(":"));

				//basic check because we have not yet checked reactions
				if (reactants.size() != CAMKIIBINDINGREACTANTS ||
					products.size() != CAMKIIBINDINGPRODUCTS) {
					cout << "Invalid CaMKII binding reaction. Exiting." << " "
						 << reactants.size() << " " << CAMKIIBINDINGREACTANTS << " "
						 << products.size() << " " << CAMKIIBINDINGPRODUCTS << endl;
					exit(EXIT_FAILURE);
				}

				if (sb_bound == sb) {
					//look at bound species associated
					string bound = reactants[1].substr(0, reactants[1].find(":"));

					auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(),
								   bound);

					//quick check for validity
					if (it == _chemData.speciesBound[filType].end()) {
						cout <<
							 "A bound species that was included in a reaction was not initialized. Exiting."
							 << endl;
						exit(EXIT_FAILURE);
					}

					// add a single binding species with name sb + bound
					// added "BindingTemplate" to describe bundled CaMKII
					protoCompartment.addSpeciesSingleBinding(SpeciesNamesDB::genBindingName("BindingTemplate" + sb, bound));
				}
			}
			//look at camkii bundling reaction that is associated with this species
			for (auto &rb : _chemData.camkiibundlingReactions[filType]) {

				auto reactants = get<0>(rb);
				auto products = get<1>(rb);

				auto sb_bound = products[0].substr(0, products[0].find(":"));

				//basic check because we have not yet checked reactions
				if (reactants.size() != CAMKIIBUNDLINGREACTANTS ||
					products.size() != CAMKIIBUNDLINGPRODUCTS) {
					cout << "Invalid CaMKII bundling reaction. Exiting." << endl;
					exit(EXIT_FAILURE);
				}

				if (sb_bound == sb) {
					//look at bound species associated
					string bound = reactants[1].substr(0, reactants[1].find(":"));

					auto it = find(_chemData.speciesBound[filType].begin(), _chemData.speciesBound[filType].end(),
								   bound);

					//quick check for validity
					if (it == _chemData.speciesBound[filType].end()) {
						cout <<
							 "A bound species that was included in a reaction was not initialized. Exiting."
							 << endl;
						exit(EXIT_FAILURE);
					}

					//add a single binding species with name sb + bound
					protoCompartment.addSpeciesSingleBinding(SpeciesNamesDB::genBindingName(sb, bound));
				}
			}
		}

	}
}
