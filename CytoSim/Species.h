//
//  Species.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/20/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

/**
 * @file Species.h
 * @brief this header file will contain defininition of the Species hieararchy and associated DB and helper classes.
 * @author Garegin Papoian *
 * @date 5/12/2012 
 */

#ifndef CytoSim_Experimenting_Species_h
#define CytoSim_Experimenting_Species_h

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>
#include <cassert>
#include <stdexcept>

#include <type_traits>

#include "common.h"
#include "utility.h"
#include "RSpecies.h"

class System;

namespace chem {
    
    class Composite;
    
    /// SpeciesNamesDB class is used to associate unique integers with character based names of Species
    /*! Often Species of the same type, let's say "Arp2/3" can be found in different forms, for example
     *  in cytosol vs bound to a filament. The corresponding specific Species will be distinct, however,
     *  they will share the same name, as returned by Species::getName(). This, in turn, helps to match 
     *  "Arp2/3" bound to filaments with diffusing "Arp2/3". SpeciesNamesDB associates a unique integer 
     *  with Species of the same type (e.g. "Arp2/3", regardless whether it is bound or not), and provides 
     *  conversion functions fron integer to 
     *  std::string (SpeciesNamesDB::intToString(int i)) and vice versa (SpeciesNamesDB::stringToInt (std::string name)). 
     *  class SpeciesNamesDB is a singleton, and should be used by first calling SpeciesNamesDB::Instance() method.
     *
     *  @code
     *  int y = SpeciesNamesDB::Instance()->stringToInt("Arp2/3"); // let's say y=2
     *  std::string x = SpeciesNamesDB::Instance()->intToString(2); // then x should be "Arp2/3"
     *  @endcode
     *
     * SpeciesNamesDB also has a function to generate a unique name given a string seed. This is particularly
     * useful for when a filament species needs a unique name, say with the seed "Actin"
     *
     */  
    class SpeciesNamesDB {
    private: 
        static SpeciesNamesDB* _instance; ///< the singleton instance
        SpeciesNamesDB() {}
    private: 
        std::unordered_map<std::string,int> _map_string_int;
        std::vector<std::string> _vec_int_string;
        unsigned long _num = 0; ///<used to generate unique names
    public:
        /// returns the unique instance of the singleton, which can be used to access the names DB
        static SpeciesNamesDB* Instance();
        
        /// Given an integer "i", returns the std::string associated with that integer. Throws an out of range exception.  
        std::string intToString(unsigned int i) const {
            if (i>=_vec_int_string.size())
                throw std::out_of_range("SpeciesNamesDB::intToString(int i) index error:[" + std::to_string(i) +"], while the vector size is " + std::to_string(_vec_int_string.size()));
            return _vec_int_string[i];
        }
        
        /// Given a std::string "name", returns the unique integer associated with that string. If the 
        /// string does not exist yet, it is created.
        int stringToInt (std::string name) {
            auto mit = _map_string_int.find(name);
            if(mit == _map_string_int.end()){
                _vec_int_string.push_back(name);
                _map_string_int[name]=_vec_int_string.size()-1;
                return _vec_int_string.size()-1;
            }
            else
                return mit->second;
        }

        ///Generate a unique name based on a seed name (just adds integer value to end of string)
        std::string generateUniqueName(std::string name) {
            std::string uniqueName = name + std::to_string(_num);
            if(_map_string_int.find(uniqueName) != _map_string_int.end())
                return uniqueName;
            else
                return name + std::to_string(++_num);
        }

        /// Clear the contents of the database
        void clear() {
            _map_string_int.clear();
            _vec_int_string.clear();
        }
    };
        
    /// Species class represents chemical molecules, tracks their copy number and can be used in [Reactions](@ref Reaction).
    /*! This class represents chemical species, such as G-Actin. As a class, it provides Species' name, the current
     *  copy number of molecules, Species type (e.g. SType::Bulk), and a few other characteristics of a Species. 
     *  This class is synergetic with the Reaction, since Species can be added to [Reactions](@ref Reaction). 
     *
     *  As a type, Species is composed of the following primary fields: type, name, copy number. A subclass may added 
     *  additional primary fields. Two Species are mathematically equal, if their primary fields are equal. This 
     *  means that applying the copy constructor will guarantee that primary fields will be equal:
     *  @code
     *  SpeciesBulk A{"Arp2/3",25};
     *  SpeciesBulk B(A);
     *  assert(A==B); // A must be equal to B
     *  @endcode
     *
     *  @note Each Species owns an unshared RSpecies object. However, the RSpecies field is not among the primary fields
     *  defining the Species identity (hence, the equality operator). In particular, upon copying, the source and target 
     *  Species will have distinct RSpecies fields. However, in a fully constructed program, the source and target species 
     *  (after applying the copy constructor) should eventually be involved in the same set of equivalent reactions. This
     *  means that their respective reactions will return true when the equlaity operator is applied.
     *
     *  @note The Species class allows callbacks (see makeSignaling and related methods). 
    */
    class Species {
    private: //Variables
        int _molecule; ///< unique id identifying the molecule (e.g. the integer id corresponding to "Arp2/3") 
        RSpecies* _rspecies; ///< pointer to RSpecies; Species is responsible for creating and destroying RSpecies
        Composite *_parent; ///< pointer to the "container" object holding this Species (could be a nullptr)
    public:
        /// Default Constructor; Should not be used by the end users - only internally (although it is not marked as private)
        Species()  : _parent(nullptr) {
            _molecule=SpeciesNamesDB::Instance()->stringToInt("");
            _rspecies = new RSpecies(*this);
//            std::cout << "Species(): Default ctor called, creating ptr=" << this << std::endl;
        }
        
        /// The constructor for this base class of Species should not be called directly - only by the concrete subclasses
        /// @param name - a string for the Species name associated with this Species. For example, "G-Actin" or "Arp2/3"
        /// @param type_enum - SType enum, such as SType::Diffusing
        /// @param n - copy number
        /// @param ulim - upper limit for this species' copy number
        Species (const std::string &name, species_copy_t n=0, species_copy_t ulim=max_ulim)  : _parent(nullptr)
        {
            _molecule=SpeciesNamesDB::Instance()->stringToInt(name);
            _rspecies = new RSpecies(*this, n, ulim);
//            std::cout << "Species (const std::string &name, species_copy_t n=0): Main ctor called, creating ptr=" << this << std::endl;
            
        }
        
        /// Copy constructor
        /// @note The associated RSpecies subpart is not copied, but a new one is created. This means that 
        /// the copied destination Species won't be included in any Reaction interactions of the original 
        /// source Species. The species copy numbered is copied to the target.
        /// The A Species parent attriute is not copied, but set to nullptr. 
        Species (const Species &rhs)  : _molecule(rhs._molecule), _parent(nullptr) {
#ifdef TRACK_UPPER_COPY_N
            _rspecies = new RSpecies(*this, rhs.getN(), rhs.getUpperLimitForN());
#else
            _rspecies = new RSpecies(*this, rhs.getN());
#endif
//            std::cout << "Species(const Species &rhs): copy constructor called, old ptr=" << &rhs << ", new ptr=" << this << std::endl;
        }
        
        /// Move constructor - makes it possible to easily add Species to STL containers, such as vector<Species>
        /// One has to be careful with vector<Species> as opposed to vector<Species*>. The latter is "safe" in 
        /// terms of being internally moved around by vector.resize, etc., but vector<Species> would copy Species 
        /// if a move constructor is not available. This will then destruct the original Species (with its 
        /// associated RSpecies), hence, the Reaction network of the Species will be lost. Since this behavior is 
        /// most likely not desired or unacceptable, the internal vector operations should only "move" 
        /// the Species around, without copying. Moving transfers the RSpecies pointer from source to target, 
        /// stealing resources from the source, leaving it for destruction. The Species parent attriute is moved it. 
        Species (Species &&rhs) noexcept
        : _molecule(rhs._molecule), _rspecies(rhs._rspecies), _parent(rhs._parent) {
            rhs._rspecies = nullptr;
//            std::cout << "Species(Species &&rhs): move constructor called, old ptr=" << &rhs << ", new ptr=" << this << std::endl;
        }
        
        /// Assignment operator
        /// An assignment A = B copies the name of B to A. It destroys the Reaction interactions of A, and resents 
        /// them to a blank value (i.e. A won't be involced in any Reactions). B's Reaction interactions are not copied.
        /// The copy number of B is copied to A. The A Species parent attriute is not copied, but set to nullptr. 
        Species& operator=(const Species& rhs)  {
//            std::cout << "Species& operator=(const Species& rhs):" << this << ", " << &rhs << std::endl; 
            _molecule = rhs._molecule;
#ifdef TRACK_UPPER_COPY_N
            _rspecies = new RSpecies(*this, rhs.getN(), rhs.getUpperLimitForN());
#else
            _rspecies = new RSpecies(*this, rhs.getN());
#endif
            _parent = nullptr;
            return *this;
        }
        
        /// Move assignment is needed for the same reasons as move constructor.
        /// @see Species (Species &&rhs)
        Species& operator=(Species&& rhs)  {
//            std::cout << "Species& operator=(Species&& rhs):" << this << ", " << &rhs << std::endl; 
            _molecule = rhs._molecule;
            _rspecies = rhs._rspecies;
            rhs._rspecies = nullptr;
            _parent=rhs._parent;
            return *this;
        }
        
        virtual Species* clone() {
            return new Species(*this);
        }
        
        Composite* getParent() {return _parent;}
        
        void setParent (Composite *other) {_parent=other;}
        
        bool hasParent() const {return _parent!=nullptr? true : false;}
        
        Composite* getRoot();
        
        /// Return a reference to RSpecies. Notice that value copying won't be allowed 
        /// because RSpecies is not copyable.
        RSpecies& getRSpecies () {return *_rspecies;}
        
        /// Return a constant reference to RSpecies. 
        const RSpecies& getRSpecies () const {return *_rspecies;}
        
        /// Sets the copy number for this Species. 
        /// @param n should be a non-negative number, but no checking is done in run time
        /// @note The operation does not emit any signals about the copy number change.
        void setN(species_copy_t n) {_rspecies->setN(n);}
        
        /// Return the current copy number of this Species
        species_copy_t getN () const {return _rspecies->getN();}
        
#ifdef TRACK_UPPER_COPY_N
        /// Return the upper limit for the copy number of this Species
        species_copy_t getUpperLimitForN() const {return _rspecies->getUpperLimitForN();}
#endif
        
        /// Return this Species' name
        std::string getName() const {return SpeciesNamesDB::Instance()->intToString(_molecule);}
        
        /// Return the molecule index associated with this Species' (as int)
        int getMolecule() const {return _molecule;}
        
#ifdef RSPECIES_SIGNALING
        /// Return true if this Species emits signals on copy number change
        bool isSignaling () const {return _rspecies->isSignaling();}
        
        /// Set the signaling behavior of this Species
        /// Gillespie-like simulation algorithm)
        void startSignaling () {_rspecies->startSignaling();}
        
        /// Destroy the signal associated with this Species
        /// @note To start signaling again, makeSignaling(...) needs to be called
        void stopSignaling () {_rspecies->stopSignaling();}
        
        /// Connect the callback, rspecies_callback to a signal corresponding to RSpecies *s.
        /// @param std::function<void (RSpecies *, int)> const &RSpecies_callback - a function object to be called (a slot)
        /// @param int priority - lower priority slots will be called first. Default is 5 Do not use priorities 1 and 2 
        ///                       unless absolutely essential.
        /// @return a connection object which can be used to later disconnect this particular slot or temporarily block it
        boost::signals2::connection connect(std::function<void (RSpecies *, int)> const &RSpecies_callback, int priority=5);
#endif
        
        /// Returns true if two Species objects are equal.
        /// This function would accept derived class of Species, such as SpeciesBulk
        /// Two Species are equal if their SType(s) are equal (i.e. are of the same class),
        /// their names are equal.
        friend bool operator ==(const Species& a, const Species& b)
        {
            if (a.getMolecule() != b.getMolecule() or typeid(a) != typeid(b))
                return false;
            return true;
        }
        
        /// Return true if two Species are not equal.         
        /// @see operator ==(const Species& a, const Species& b) above
        friend bool operator !=(const Species& a, const Species& b){
            return !(a==b);
        }
        
        // virtual methods
        
        /// Virtual destructor
        /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
        /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug 
        /// (as of gcc 4.703), and will presumbaly be fixed in the future.
        virtual ~Species () noexcept
        {
            if(_rspecies!=nullptr)
                delete _rspecies;
        };
                
        /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Bulk}"
        virtual std::string getFullName() const {return getName() + "{None}";}
        
        virtual size_t countSpecies() const {return 1;}
    };
    
    /// SpeciesBulk should be used for Species without spatial information (i.e. well-mixed in the container)
    class SpeciesBulk : public Species {
    public:
        /// Default constructor
        SpeciesBulk()  : Species() {}
        
        /// The main constructor 
        /// @param name - Example, "G-Actin" or "Arp2/3"
        /// @param n - copy number
        SpeciesBulk (const std::string &name, species_copy_t n=0, species_copy_t ulim=max_ulim)  :  Species(name, n, ulim) {};
        
        /// Copy constructor
        SpeciesBulk (const SpeciesBulk &rhs)  : Species(rhs) {}
        
        /// Move constructor
        SpeciesBulk (SpeciesBulk &&rhs) noexcept : Species(std::move(rhs)) {
        }
        
        /// Regular Assignment
        SpeciesBulk& operator=(const SpeciesBulk& rhs)  {
            Species::operator=(rhs);
            return *this;
        }
        
        /// Move assignment
        SpeciesBulk& operator=(SpeciesBulk&& rhs) 
        {
            Species::operator=(std::move(rhs));
            return *this;
        }
        
        virtual SpeciesBulk* clone() {
            return new SpeciesBulk(*this);
        }
        
        /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Bulk}"
        virtual std::string getFullName() const {return getName() + "{Bulk}";}
        
        /// Default destructor
        ~SpeciesBulk () noexcept {};
    };
    
    /// SpeciesDiffusing should be used for Species which can move spatially from one compartment to the neighboring one (i.e. they are the stochastic analogue of determenistic reaction-diffusion processes)
    class SpeciesDiffusing : public Species {
    public:
        /// Default constructor
        SpeciesDiffusing()  : Species() {}
        
        /// The main constructor 
        /// @param name - Example, "G-Actin" or "Arp2/3"
        /// @param n - copy number
        SpeciesDiffusing (const std::string &name, species_copy_t n=0, species_copy_t ulim=max_ulim)  :  Species(name, n, ulim) {};
        
        /// Copy constructor
        SpeciesDiffusing (const SpeciesDiffusing &rhs)  : Species(rhs) {}
        
        /// Move constructor
        SpeciesDiffusing (SpeciesDiffusing &&rhs) noexcept : Species(std::move(rhs)) {
        }
        
        /// Regular Assignment
        SpeciesDiffusing& operator=(const SpeciesDiffusing& rhs)  {
            Species::operator=(rhs);
            return *this;
        }
        
        /// Move assignment
        SpeciesDiffusing& operator=(SpeciesDiffusing&& rhs) 
        {
            Species::operator=(std::move(rhs));
            return *this;
        }
        
        virtual SpeciesDiffusing* clone() {
            return new SpeciesDiffusing(*this);
        }
        
        /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Diffusing}"
        virtual std::string getFullName() const {return getName() + "{Diffusing}";}
        
        /// Default destructor
        ~SpeciesDiffusing () noexcept {};
    };
    
        /// SpeciesFilament should be used for species that can be in a filament. These species can not move
    /// cross-compartment.
    class SpeciesFilament : public Species {
    public:
        /// Default constructor
        SpeciesFilament()  : Species() {}
        
        /// The main constructor
        /// @param name - Example, "G-Actin" or "Arp2/3"
        /// @param n - copy number
        SpeciesFilament (const std::string &name, species_copy_t n=0, species_copy_t ulim=1)
            :  Species(name, n, ulim) {};
        
        /// Copy constructor
        SpeciesFilament (const SpeciesFilament &rhs)  : Species(rhs) {}
        
        /// Move constructor
        SpeciesFilament (SpeciesFilament &&rhs) noexcept : Species(std::move(rhs)) {
        }
        
        /// Regular Assignment
        SpeciesFilament& operator=(const SpeciesFilament& rhs)  {
            Species::operator=(rhs);
            return *this;
        }
        
        /// Move assignment
        SpeciesFilament& operator=(SpeciesFilament&& rhs)
        {
            Species::operator=(std::move(rhs));
            return *this;
        }
        
        virtual SpeciesFilament* clone() {
            return new SpeciesFilament(*this);
        }
        
        /// Return the full name of this Species in a std::string format (e.g. "Actin{Filament}"
        virtual std::string getFullName() const {return getName() + "{Filament}";}
        
        /// Default destructor
        ~SpeciesFilament () noexcept {};
    };
    
    /// SpeciesBound should be used for species that can be bound to a filament.
    /// These species can not move cross-compartment.
    class SpeciesBound : public Species {
    public:
        /// Default constructor
        SpeciesBound()  : Species() {}
        
        /// The main constructor
        /// @param name - Example, "G-Actin" or "Arp2/3"
        /// @param n - copy number
        SpeciesBound (const std::string &name, species_copy_t n=0, species_copy_t ulim=1)
        :  Species(name, n, ulim) {};
        
        /// Copy constructor
        SpeciesBound (const SpeciesBound &rhs)  : Species(rhs) {}
        
        /// Move constructor
        SpeciesBound (SpeciesBound &&rhs) noexcept : Species(std::move(rhs)) {
        }
        
        /// Regular Assignment
        SpeciesBound& operator=(const SpeciesBound& rhs)  {
            Species::operator=(rhs);
            return *this;
        }
        
        /// Move assignment
        SpeciesBound& operator=(SpeciesBound&& rhs)
        {
            Species::operator=(std::move(rhs));
            return *this;
        }
        
        virtual SpeciesBound* clone() {
            return new SpeciesBound(*this);
        }
        
        /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Bound}"
        virtual std::string getFullName() const {return getName() + "{Bound}";}
        
        /// Default destructor
        ~SpeciesBound () noexcept {};
    };
  
    
} // end of chem namespace 

/// Print self into an iostream
std::ostream& operator<<(std::ostream& os, const chem::Species& s);

#endif
