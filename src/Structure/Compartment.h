
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

#ifndef MEDYAN_Compartment_h
#define MEDYAN_Compartment_h

#include <vector>
#include <array>
#include <unordered_map>
#include <unordered_set>

#include "common.h"

#include "SpeciesContainer.h"
#include "ReactionContainer.h"
#include "BindingManager.h"
#include "HybridBindingSearchManager.h"
#include "Composite.h"
#include "ChemSim.h"
#ifdef SIMDBINDINGSEARCH
#include "dist_moduleV2/dist_driver.h"
#include "dist_moduleV2/dist_coords.h"
#endif
#include "MathFunctions.h"
#include "Structure/CellList.hpp"

//#include "BinGrid.h"
//#include "Bin.h"
//FORWARD DECLARATIONS
class BoundaryElement;
class Bead;
class Cylinder;

/// A container or holding Species and [Reactions](@ref Reactions).

/*! The Compartment class is a container for Species, internal [Reactions](@ref Reactions),
 *  and diffusion [Reactions](@ref Reactions) that can occur. A Compartment object keeps
 *  track of the above while also holding pointers to its neighbors in order to generate
 *  diffusion [Reactions](@ref Reactions) and control other interactions.
 *
 *  The Compartment also keeps Trackable elements in the SubSystem that are in its space, including
 *  [Beads](@ref Bead), [Cylinders](@ref Cylinder), and [BoundaryElements](@ref BoundaryElement).
 *
 *  Lastly, the Compartment holds a container of FilamentBindingManager for updating
 *  binding reactions local to this compartment space only.
 */

class Compartment : public Composite {

protected:
    ///CHEMICAL CONTAINERS
    SpeciesPtrContainerVector _species;  ///< Container with all species
                                         ///< in this compartment
    ReactionPtrContainerVector _internal_reactions; ///< Container with all internal
                                                    ///< reactions in compartment
    ReactionPtrContainerVector _diffusion_reactions; ///< Container with all diffusion
                                                     ///< reactions in compartment

    unordered_map<int,float> _diffusion_rates; ///< Diffusion rates of Species
                                               ///< in compartment



    /// All binding managers for this compartment
    vector<unique_ptr<FilamentBindingManager>> _bindingManagers;
    vector<unique_ptr<FilamentBindingManager>> _branchingManagers;

#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    //Each compartment has a single instance of hybridbindingmanagers
    HybridBindingSearchManager* _bindingsearchManagers = NULL;
#endif
    ///ELEMENT CONTAINERS
    unordered_set<BoundaryElement*> _boundaryElements; ///< Set of boundary element
                                                       ///< that are in this compartment

    unordered_set<Bead*> _beads; ///< Set of beads that are in this compartment

    vector<Compartment*> _neighbours; ///< Neighbors of the compartment (neighbors that
    unordered_map<Compartment*, size_t> _neighborIndex; ///< Spacial index of the neighbors of the same order as _neighbors
    ///< In 3D, the indices are in the order (x-, x+, y-, y+, z-, z+)

    // touch along faces
    vector<Compartment*> _enclosingneighbours; ///< Neighbors that envelop the compartment
    vector<Compartment*> _uniquepermuteneighbours; //Neighbors tht help to parse
    // through unique pairs of compartments to get all necessary binding sites.
    vector<short> _uniquepermuteneighboursstencil;
    vector<short> _enclosingneighboursstencil;///< enclosing compartments  have unique
    // position IDs between 0-26. This ID immediately helps us determine relative
    // position of the Enclosing Compartment with respect to the current compartment.


    ///OTHER COMPARTMENT PROPERTIES
    vector<floatingpoint> _coords;  ///< Coordinates of this compartment
    bool _activated = false; ///< The compartment is activated for diffusion

    floatingpoint _partialVolume = 1.0; ///< The volume fraction inside the
    // membrane/boundary
    ///< Might be changed to a list or a map when more membranes are involved
    array<floatingpoint, 6> _partialArea {{1.0, 1.0, 1.0, 1.0, 1.0, 1.0}}; ///< The area
    // inside the cell membrane
    ///<In the order of x, y and z, from smaller coordinate value neighbor to larger coordinate value
    ///< Might be changed to a list of arrays or a map of arrays when more membranes are involved

public:
    short _ID;
/*    vector<uint16_t> cindex_bs;
    vector<uint32_t> cID_bs;*/
#ifdef MOTORBIASCHECK
    size_t nummotorwalks=0;
#endif
    cell_list::CellListHeadUser< Cylinder, Compartment > cylinderCell; // Cell of cylinders

    /// Default constructor, only takes in number of dimensions
    Compartment() : _species(), _internal_reactions(),
    _diffusion_reactions(), _diffusion_rates(), _neighbours()  {
        }

    /// Constructor which clones another compartment
    Compartment(const Compartment &C) : _species(), _internal_reactions(),
                                        _diffusion_reactions(), _neighbours()
    {
        C.cloneSpecies(this);
        C.cloneReactions(this);
        _diffusion_rates = C._diffusion_rates;
        _activated = C._activated;

        // Should eventually clone beads, cylinders, boundary elements.... not clear yet
    }

    /// Assignment operator
    Compartment& operator=(const Compartment &other);

    /// Destructor
    /// @note noexcept is important here. Otherwise, gcc flags the constructor as
    /// potentially throwing, which in turn disables move operations by the STL
    /// containers. This behaviour is a gcc bug (as of gcc 4.703), and will presumbaly
    /// be fixed in the future.
    virtual ~Compartment() noexcept
    {
        clearNeighbouringReactions();
        clearReactions();
        clearSpecies();
        removeFromNeighboursList();

        // Should eventually delete beads, cylinders, boundary elements....not yet clear
    }

    /// get ID
    virtual int getId(){return _ID;}
    /// Applies SpeciesVisitor v to every Species* object directly owned by this node.
    /// This method needs to be overriden by descendent classes that contain Species.
    virtual bool apply_impl(SpeciesVisitor &v) override;

    /// Applies ReactionVisitor v to every ReactionBase* object directly owned by this
    /// node. This method needs to be overriden by descendent classes that contain
    /// ReactionBase.
    virtual bool apply_impl(ReactionVisitor &v) override;

    ///Set a compartment as active. Used at initialization.
    virtual void setAsActive() {_activated = true;}

    /// Activate a compartment. Has the following side effects:
    /// 1) Adds diffusion reactions between this compartment's diffusing
    ///    species and its neighboring active compartments
    virtual void activate(ChemSim* chem);

    /// Deactivate a compartment. Has the following sid effects:
    /// 0) Initially checks that all cylinders are removed
    ///    from this compartment. A compartment cannot be deactivated
    ///    unless this condition is already true.
    /// 1) Transfers copy numbers of all diffusing species from this
    ///    compartment and moves them to a neighboring active compartment.
    ///    If there are no neighboring active compartments, an error will result.
    /// 2) Removes all diffusion reactions involving diffusing species
    ///    in this compartment.
    virtual void deactivate(ChemSim* chem);

    ///Check if compartment is activated
    virtual bool isActivated() {return _activated;}

    ///Setter and getter for fs
    virtual void setCoordinates(vector<floatingpoint> coords) {_coords = coords;}
    virtual const vector<floatingpoint>& coordinates() {return _coords;}


    /// Transfer all species copy numbers from this compartment to neighboring
    /// active compartments. If no neighboring active compartments are present,
    /// throw an error.
    virtual void transferSpecies(int i);
    virtual void shareSpecies(int i);


    /// Removes all reactions from this compartment, diffusing and internal
    virtual void clearReactions() {

        _internal_reactions.reactions().clear();
        _diffusion_reactions.reactions().clear();
    }

    /// Clear all species from this compartment
    virtual void clearSpecies() {
        _species.species().clear();
    }

    /// Remove diffusion reactions between this compartment and its neighbors
    virtual void clearNeighbouringReactions() {

        for(auto &s : _species.species())
        {
            for(auto &n : _neighbours)
            {
                n->removeDiffusionReactions(s.get());
                n->removeInternalReactions(s.get());
            }
        }
    }
    /// Remove this compartment from the neighbor list
    virtual void removeFromNeighboursList() {

        for(auto &n : _neighbours)
        {
            n->removeNeighbour(this);
        }
    }

    /// Returns true if two Compartment objects are equal.
    /// Two Compartment objects are equal if each contains analogous Species and
    /// Reaction objects, in the same order
    friend bool operator==(const Compartment& a, const Compartment& b);

    /// Return true if two Compartment are not equal.
    /// @see operator ==(const Compartment& a, const Compartment& b) above
    friend bool operator !=(const Compartment& a, const Compartment& b)
    {
        return !(a==b);
    }

    /// This is a species container
    virtual bool isSpeciesContainer() const override {return true;}
    /// This is a reaction container
    virtual bool isReactionsContainer() const override {return true;}
    /// Returns compartment name
    virtual string getFullName() const override {return "Compartment";};
    /// Returns the number of species in this compartment
    size_t numberOfSpecies() const override {return _species.species().size();}
    /// Returns the number of internal reactions in this compartment
    size_t numberOfInternalReactions() const {
        return _internal_reactions.reactions().size();
    }
    /// Returns the total number of reactions in this compartment, diffusing and
    /// internal
    size_t numberOfReactions() const override {
        return _internal_reactions.reactions().size() +
               _diffusion_reactions.reactions().size();
    }

    //@{
    /// Species finder functions
    Species* findSpeciesByName(const string &name) {
        return _species.findSpeciesByName(name);
    }
    Species* findSpeciesByIndex (size_t index) {
        return _species.findSpeciesByIndex(index);
    }
    Species* findSpeciesByMolecule (int molecule) {
        return _species.findSpeciesByMolecule(molecule);
    }
    Species* findSimilarSpecies (const Species &s) {
        return _species.findSimilarSpecies(s);
    }
    //@}

    ///Remove species from this compartment
    size_t removeSpecies(Species* species) {return _species.removeSpecies(species);}

    /// Finds a similar internal reaction, see ReactionBase function
    ReactionBase* findSimilarInternalReaction (const ReactionBase &r) {
        return _internal_reactions.findSimilarReaction(r);
    }

    /// Finds a similar diffusion reaction, see ReactionBase function
    ReactionBase* findSimilarDiffusionReaction (const ReactionBase &r) {
        return _diffusion_reactions.findSimilarReaction(r);
    }

    /// Remove all diffusion reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeDiffusionReactions (Species* s) {
        _diffusion_reactions.removeReactions(s);
    }

    /// Remove all internal reactions that have a given species
    /// @param s - species whose reactions should be removed
    virtual void removeInternalReactions (Species* s) {
        _internal_reactions.removeReactions(s);
    }

    /// Remove a diffusion reaction
    virtual void removeDiffusionReaction(ReactionBase *r) {
        _diffusion_reactions.removeReaction(r);
    }

    /// Remove an internal reaction
    virtual void removeInternalReaction(ReactionBase *r) {
        _internal_reactions.removeReaction(r);
    }

    /// Add a unique species pointer to this compartment
    Species* addSpeciesUnique (unique_ptr<Species> &&species, float diff_rate = -1.0) {
        Species *sp = _species.addSpeciesUnique(move(species));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=diff_rate;
        return sp;
    }

    /// Add a unique internal reaction pointer to this compartment
    ReactionBase* addInternalReactionUnique (unique_ptr<ReactionBase> &&reaction) {
        ReactionBase *r = _internal_reactions.addReactionUnique(move(reaction));
        r->setParent(this);
        return r;
    }

    /// Add a unique diffusing reaction pointer to this compartment
    ReactionBase* addDiffusionReactionUnique (unique_ptr<ReactionBase> &&reaction) {
        ReactionBase *r = _diffusion_reactions.addReactionUnique(move(reaction));
        r->setParent(this);
        return r;
    }

    /// Add an internal reaction pointer to this compartment. Make unique
    ReactionBase* addInternalReaction (ReactionBase* r) {
        _internal_reactions.addReactionUnique(unique_ptr<ReactionBase>(r));
        r->setParent(this);
        return r;
    }

    /// Add a diffusion reaciton pointer to this compartment. Make unique
    ReactionBase* addDiffusionReaction (ReactionBase* r) {
        _diffusion_reactions.addReactionUnique(unique_ptr<ReactionBase>(r));
        r->setParent(this);
        return r;
    }

    /// Add a diffusing species to this compartment
    /// @param args - any number of SpeciesDiffusing objects
    template<typename ...Args>
    SpeciesDiffusing* addSpeciesDiffusing(Args&& ...args) {
        SpeciesDiffusing *sp =
        (SpeciesDiffusing*)(_species.addSpecies<SpeciesDiffusing>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a filament species to this compartment
    /// @param args - any number of SpeciesFilament objects
    template<typename ...Args>
    SpeciesFilament* addSpeciesFilament(Args&& ...args) {
        SpeciesFilament *sp =
        (SpeciesFilament*)(_species.addSpecies<SpeciesFilament>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a plus end species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesPlusEnd* addSpeciesPlusEnd(Args&& ...args) {
        SpeciesPlusEnd *sp =
        (SpeciesPlusEnd*)(_species.addSpecies<SpeciesPlusEnd>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a minus end species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesMinusEnd* addSpeciesMinusEnd(Args&& ...args) {
        SpeciesMinusEnd *sp =
        (SpeciesMinusEnd*)(_species.addSpecies<SpeciesMinusEnd>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a bound species to this compartment
    /// @param args - any number of SpeciesBound objects
    template<typename ...Args>
    SpeciesBound* addSpeciesBound(Args&& ...args) {
        SpeciesBound *sp =
        (SpeciesBound*)(_species.addSpecies<SpeciesBound>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a linker species to this compartment
    /// @param args - any number of SpeciesLinker objects
    template<typename ...Args>
    SpeciesLinker* addSpeciesLinker(Args&& ...args) {
        SpeciesLinker *sp =
        (SpeciesLinker*)(_species.addSpecies<SpeciesLinker>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a motor species to this compartment
    /// @param args - any number of SpeciesMotor objects
    template<typename ...Args>
    SpeciesMotor* addSpeciesMotor(Args&& ...args) {
        SpeciesMotor *sp =
        (SpeciesMotor*)(_species.addSpecies<SpeciesMotor>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a brancher species to this compartment
    /// @param args - any number of SpeciesBrancher objects
    template<typename ...Args>
    SpeciesBrancher* addSpeciesBrancher(Args&& ...args) {
        SpeciesBrancher *sp =
        (SpeciesBrancher*)(_species.addSpecies<SpeciesBrancher>(forward<Args>(args)...));
        sp->setParent(this);
        _diffusion_rates[sp->getMolecule()]=-1.0;
        return sp;
    }

    /// Add a single binding species to this compartment
    /// @param args - any number of SpeciesSingleBinding objects
    template<typename ...Args>
    SpeciesSingleBinding* addSpeciesSingleBinding(Args&& ...args) {
        SpeciesSingleBinding *sb =
        (SpeciesSingleBinding*)(_species.addSpecies<SpeciesSingleBinding>(forward<Args>(args)...));
        sb->setParent(this);
        _diffusion_rates[sb->getMolecule()]=-1.0;
        return sb;
    }

    /// Add a pair binding species to this compartment
    /// @param args - any number of SpeciesPairBinding objects
    template<typename ...Args>
    SpeciesPairBinding* addSpeciesPairBinding(Args&& ...args) {
        SpeciesPairBinding *sb =
        (SpeciesPairBinding*)(_species.addSpecies<SpeciesPairBinding>(forward<Args>(args)...));
        sb->setParent(this);
        _diffusion_rates[sb->getMolecule()]=-1.0;
        return sb;
    }

    /// Add an internal reaction to this compartment
    template<unsigned short M, unsigned short N, typename ...Args>
    ReactionBase* addInternalReaction (Args&& ...args) {
        ReactionBase *r = _internal_reactions.addReaction<M,N>(forward<Args>(args)...);
        r->setParent(this);
        return r;
    }

    /// Add an internal reaction to this compartment
    /// @param species, rate - specifying the species and rate that should be assigned
    template<template <unsigned short M, unsigned short N> class RXN, unsigned short M, unsigned short N>
    ReactionBase* addInternal(initializer_list<Species*> species, float rate) {
        ReactionBase *r = _internal_reactions.add<RXN,M,N>(species,rate);
        r->setParent(this);
        return r;
    }

    /// Add a diffusion reaction to this compartment
    template<typename ...Args>
    ReactionBase* addDiffusionReaction (Args&& ...args) {
        ReactionBase *r = _diffusion_reactions.addReaction<1,1>(forward<Args>(args)...);
        r->setParent(this);
        return r;
    }

    ///Add branching manager to this compartment. Used only in SIMD case
    void addBranchingBindingManager(FilamentBindingManager* m){
    	_branchingManagers.emplace_back(m);
    }

    /// Add a binding manager to this compartment
    void addFilamentBindingManager(FilamentBindingManager* m) {
        _bindingManagers.emplace_back(m);
    }
#if defined(HYBRID_NLSTENCILLIST) || defined(SIMDBINDINGSEARCH)
    void addHybridBindingSearchManager(HybridBindingSearchManager* bsm){
        if(_bindingsearchManagers == NULL)
            _bindingsearchManagers = bsm;
        else{
            cout<<"Hybrid Binding Search Manager exists. Compartment cannot have multiple"
                    " hybrid binding search managers. Exiting.."<<endl;
            exit(EXIT_FAILURE);
        }
    }

    HybridBindingSearchManager* getHybridBindingSearchManager(){
        return _bindingsearchManagers;
    }
#endif
#ifdef SIMDBINDINGSEARCH
    dist::Coords bscoords;
    vector<dist::Coords> bscoords_section;
	vector<dist::Coords> bscoords_section_linker;
	vector<dist::Coords> bscoords_section_motor;

    vector<int> Cyldcindexvec;
    vector<int> CylcIDvec;

    template<bool LinkerorMotor>
    dist::Coords& getSIMDcoordsV3(short i, short filamentType){
        if(LinkerorMotor)
            return bscoords_section_linker[filamentType*27 + i];
        else
            return bscoords_section_motor[filamentType*27 + i];
    }
    /*Each compartment is partitioned into 27 sub-volumes. Binding sites are allocated
    to relevant sub-volumes. It is worth noting that the 27 volumes  can be overlapping.
    Binding distance determines if the volumes are overlapping or not.*/
    vector<floatingpoint> partitionedcoordx[27], partitionedcoordy[27], partitionedcoordz[27];
    vector<uint32_t>  cindex_bs_section[27];
    vector<uint32_t> finfo_bs_section[27];

    void SIMDcoordinates_section();
    void SIMDcoordinates4linkersearch_section(bool isvectorizedgather);
    void SIMDcoordinates4motorsearch_section(bool isvectorizedgather);
    void getpartition3Dindex(int (&indices)[3], vector<floatingpoint> coord);
    template<bool rMaxvsCmpsize>
    void getpartitionindex(int (&indices)[3], vector<floatingpoint> coord,
                                floatingpoint (&cmpcornercoords)[6]);

	void addcoord(vector<floatingpoint> coord, uint32_t index, uint32_t cylfinfo, short i);
    bool checkoccupancy(Cylinder* cyl, short it, short _filamentType, short bstatepos);
    bool checkoccupancy(vector<vector<bool>>& boundstate, short bstatepos, int pos);
    void addcoordtopartitons(int (&pindices)[3], vector<floatingpoint> coord, uint32_t
    index, uint32_t cylfinfo);
    void addcoordtopartitons_smallrmax(int (&pindices)[3], vector<floatingpoint> coord,
                                  uint16_t index, uint32_t cylfinfo);
    template<bool rMaxvsCmpsize>
    void addcoordtorMaxbasedpartitons(int (&pindices)[3], vector<floatingpoint> coord,
                                       uint32_t index, uint32_t cylfinfo);

    void deallocateSIMDcoordinates();
#endif
    /// Get binding managers for this compartment
    vector<unique_ptr<FilamentBindingManager>>& getFilamentBindingManagers() {
        return _bindingManagers;
    }

    //Get BranchingBindingManager
    vector<unique_ptr<FilamentBindingManager>>& getBranchingManagers() {
        return _branchingManagers;
    }

    /// Get a specific motor binding manager from this compartment
    MotorBindingManager* getMotorBindingManager(int motorType) {

        MotorBindingManager* mp;

        for(auto it = _bindingManagers.begin(); it != _bindingManagers.end(); it++) {

            //find the correct manager and type
            if((mp = dynamic_cast<MotorBindingManager*>((*it).get())) && (*it)->getBoundInt() == motorType)
                return mp;
        }

        return nullptr;
    }

    ///Add a boundary element to this compartment
    void addBoundaryElement(BoundaryElement* be) {_boundaryElements.insert(be);}

    ///Remove a boundary element from this compartment
    ///@note does nothing if boundary element is not in compartment
    void removeBoundaryElement(BoundaryElement* be) {
        auto it = _boundaryElements.find(be);
        if(it != _boundaryElements.end()) _boundaryElements.erase(it);
    }
    ///Check if boundary element is in this container
    bool hasBoundaryElement(BoundaryElement* be) {
        auto it = _boundaryElements.find(be);
        return (it != _boundaryElements.end());
    }

    ///get the boundary elements in this compartment
   unordered_set<BoundaryElement*>& getBoundaryElements() {return _boundaryElements;}

    ///get the cylinders in this compartment
    auto getCylinders() const { return cylinderCell.manager->getElementPtrs(cylinderCell); }

    /// Get the diffusion rate of a species
    /// @param - species_name, a string
    float getDiffusionRate(string species_name) {
        int molecule = SpeciesNamesDB::stringToInt(species_name);
        return _diffusion_rates[molecule];
    }

    /// Set the diffusion rate of a species in the compartment
    void setDiffusionRate(Species *sp, float diff_rate) {
        int molecule = sp->getMolecule();
        _diffusion_rates[molecule]=diff_rate;
    }

    /// Set the diffusion rate of a species in the compartment
    /// @param - molecule, an integer representing the species in the speciesDB
    void setDiffusionRate(int molecule, float diff_rate) {
        _diffusion_rates[molecule]=diff_rate;
    }

    /// Set the diffusion rate of a species in the compartment
    /// @param - species_name, a string
    void setDiffusionRate(string species_name, float diff_rate) {
        int molecule = SpeciesNamesDB::stringToInt(species_name);
        _diffusion_rates[molecule]=diff_rate;
    }

    /// Add a neighboring compartment to this compartments list of neighbors
    void addNeighbour(Compartment *comp, size_t spatialIndex) {
        auto nit = find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit==_neighbours.end()) {
            _neighbours.push_back(comp);
            _neighborIndex[comp] = spatialIndex;
        }
        else
            throw runtime_error(
                                "Compartment::addNeighbour(): Compartment is already a neighbour");
    }

    /// Remove a neighboring compartment
    void removeNeighbour(Compartment *comp) {
        auto nit = find(_neighbours.begin(),_neighbours.end(), comp);
        if(nit!=_neighbours.end())
            _neighbours.erase(nit);
        _neighborIndex.erase(comp);
    }

    /// Add a neighboring compartment to this compartments list of neighbors
    void addenclosingNeighbour(Compartment *comp, int stencilpos) {
        auto nit = find(_enclosingneighbours.begin(),_enclosingneighbours.end(), comp);
        if(nit==_enclosingneighbours.end()) {
            _enclosingneighbours.push_back(comp);
            _enclosingneighboursstencil.push_back(stencilpos);
        }
        else
            throw runtime_error(
                    "Compartment::addenclosingNeighbour(): Compartment is already a "
                    "neighbour");
    }


    void adduniquepermuteNeighbour(Compartment *comp, int stencilpos) {
        auto nit = find(_uniquepermuteneighbours.begin(),_uniquepermuteneighbours.end(), comp);
        if(nit==_uniquepermuteneighbours.end()) {
            _uniquepermuteneighbours.push_back(comp);
            _uniquepermuteneighboursstencil.push_back(stencilpos);
        }
        else
            throw runtime_error(
                    "Compartment::addenuniquepermuteNeighbour(): Compartment is already a "
                    "neighbour");
    }

    /// Remove a neighboring compartment
    void removeenclosingNeighbour(Compartment *comp) {
        for(int i = 0; i < _enclosingneighbours.size(); i++){
            if(comp == _enclosingneighbours[i]){
                _enclosingneighbours.erase(_enclosingneighbours.begin() + i);
                _enclosingneighboursstencil.erase(_enclosingneighboursstencil.begin() + i);
                break;
            }
        }
    }

    vector<Compartment*> getenclosingNeighbours(){
        return _enclosingneighbours;
    }

    void removeuniquepermuteNeighbour(Compartment *comp) {
        for(int i = 0; i < _uniquepermuteneighbours.size(); i++){
            if(comp == _uniquepermuteneighbours[i]){
                _uniquepermuteneighbours.erase(_uniquepermuteneighbours.begin() + i);
                _uniquepermuteneighboursstencil.erase(_uniquepermuteneighboursstencil.begin() + i);
                break;
            }
        }
    }
    vector<short> getuniquepermuteneighborsstencil(){
        return _uniquepermuteneighboursstencil;
    }
    vector<Compartment*> getuniquepermuteNeighbours(){
        return _uniquepermuteneighbours;
    }

    /// Clone the species values of another compartment into this one
    void cloneSpecies (Compartment *target) const {
        assert(target->numberOfSpecies()==0);
        for(auto &s : _species.species()){
            Species* sClone = s->clone();
            target->addSpeciesUnique(unique_ptr<Species>(sClone));
        }
    }

    /// Clone the reaction values of another compartment into this one
    void cloneReactions (Compartment *target) const {
        assert(target->numberOfReactions()==0);
        for(auto &r : _internal_reactions.reactions()){

            auto rClone = r->clone(target->_species);
            rClone->setVolumeFrac(target->getVolumeFrac());
            target->addInternalReaction(rClone);
        }
    }

    /// Clone both the species and compartments into this compartment
    void cloneSpeciesReactions(Compartment* C) {
        if(_activated) this->cloneSpecies(C);
        this->cloneReactions(C);
        C->_diffusion_rates = this->_diffusion_rates;

    }

    /// Clone a compartment
    /// @note - this does not clone the neighbors, just reactions and species
    virtual Compartment* clone() {
        Compartment *C = new Compartment(*this);
        return C;
    }

    /// Generate diffusion reactions between this compartment and another
    ///@return - a vector of reactionbases that was just added
    vector<ReactionBase*> generateDiffusionReactions(Compartment* C);

    /// Generate all diffusion reactions for this compartment and its neighbors
    ///@return - a vector of reactionbases that was just added
    vector<ReactionBase*> generateAllDiffusionReactions();

    /// Generates all diffusion reactions between this compartment and its neighbors
    /// in addition to generating reverse reactions
    ///@return - a vector of reactionbases that was just added
    vector<ReactionBase*> generateAllpairsDiffusionReactions();

    /// Remove diffusion reactions between this compartment and another
    ///@return - a vector of reactionbases that was just removed
    void removeDiffusionReactions(ChemSim* chem, Compartment* C);


    /// Remove all diffusion reactions for this compartment and its neighbors
    ///@return - a vector of reactionbases that was just removed
    void removeAllDiffusionReactions(ChemSim* chem);


    /// Gives the number of neighbors to this compartment
    size_t numberOfNeighbours() const {return _neighbours.size();}

    /// Gives the number of enclosing neighbors to this compartment
    size_t numberOfenclosingNeighbours() const {return _enclosingneighbours.size();}


    ///Get the species container vector
    SpeciesPtrContainerVector& getSpeciesContainer() {return _species;}
    const SpeciesPtrContainerVector& getSpeciesContainer() const {return _species;}

    ///Get the internal reaction container vector
    ReactionPtrContainerVector& getInternalReactionContainer() {return _internal_reactions;}
    const ReactionPtrContainerVector& getInternalReactionContainer() const {return _internal_reactions;}

    ///Get the diffusion reaction container vector
    ReactionPtrContainerVector& getDiffusionReactionContainer() {return _diffusion_reactions;}
    const ReactionPtrContainerVector& getDiffusionReactionContainer() const {return _diffusion_reactions;}

    /// Get the vector list of neighbors to this compartment
    vector<Compartment*>& getNeighbours() {return _neighbours;}

    const vector<Compartment*>& getNeighbours() const {return _neighbours;}

    /// Print the species in this compartment
    void printSpecies() {_species.printSpecies();}
    /// Print the reactions in this compartment
    void printReactions() {
        _internal_reactions.printReactions();
        _diffusion_reactions.printReactions();
    }

    /// Check if all species are unique pointers
    bool areAllSpeciesUnique () {return _species.areAllSpeciesUnique();}

    /// Adds the reactions of this compartment to the ChemSim object
    /// @param - chem, a ChemSim object that runs the reaction-diffusion algorithm
    virtual void addChemSimReactions(ChemSim* chem) {
        for(auto &r : _internal_reactions.reactions()) chem->addReaction(r.get());
        for(auto &r : _diffusion_reactions.reactions()) chem->addReaction(r.get());
    }

    /// Print properties of this compartment
    virtual void printSelf() override {
        cout << this->getFullName() << "\n"
        << "Number of neighbors: " << numberOfNeighbours() << endl;
        printSpecies();
        cout << "Reactions:" << endl;
        printReactions();
    }

    //GetType implementation just returns zero (no Compartment types yet)
    virtual int getType() override {return 0;}

    // Helper function for getting the result of geometry from a approximately planar slice
    void getSlicedVolumeArea();
    // Helper function that does not scale rates
    void getNonSlicedVolumeArea();

    // Properties (public variables and getters and setters for private variables)
    bool boundaryInteresting = false; // A marker indicating this compartment is near a certain boundary

    //_partialVolume is the volume fraction
    floatingpoint getPartialVolume()const { return _partialVolume; }
    void setPartialVolume(floatingpoint partialVolume) { _partialVolume = partialVolume; }
    floatingpoint getVolumeFrac()const {
        return _partialVolume;
    }
    const array<floatingpoint, 6>& getPartialArea()const { return _partialArea; }
    void setPartialArea(const array<floatingpoint, 6>& partialArea) { _partialArea =
    partialArea; }

};
#ifdef SIMDBINDINGSEARCH
template<>
void Compartment::getpartitionindex<true>(int (&indices)[3], vector<floatingpoint> coord,
                             floatingpoint (&cmpcornercoords)[6]);
template<>
void Compartment::getpartitionindex<false>(int (&indices)[3], vector<floatingpoint> coord,
                              floatingpoint (&cmpcornercoords)[6]);
template<>
void Compartment::addcoordtorMaxbasedpartitons<true>(int (&pindices)[3], vector<floatingpoint>
        coord, uint32_t index, uint32_t cylfinfo);
template<>
void Compartment::addcoordtorMaxbasedpartitons<false>(int (&pindices)[3], vector<floatingpoint>
        coord, uint32_t index, uint32_t cylfinfo);
#endif
#endif
