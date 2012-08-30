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
#include "Component.h"
#include "RSpecies.h"

class System;

namespace chem {
    
    
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
     */  
    class SpeciesNamesDB {
    private: 
        static SpeciesNamesDB* _instance; ///< the singleton instance
        SpeciesNamesDB() {}
    private: 
        std::unordered_map<std::string,int> _map_string_int;
        std::vector<std::string> _vec_int_string;
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
    class Species : public Component {
    private: //Variables
        int _molecule; ///< unique id identifying the molecule (e.g. the integer id corresponding to "Arp2/3") 
        std::unique_ptr<RSpecies> _rspecies; ///< pointer to RSpecies; Species is responsible for creating and destroying RSpecies
    public:
        /// Default Constructor; Should not be used by the end users - only internally (although it is not marked as private)
        Species()  : Component() {
            _molecule=SpeciesNamesDB::Instance()->stringToInt("");
            _rspecies = std::unique_ptr<RSpecies>(new RSpecies(*this));
//            std::cout << "Species(): Default ctor called, creating ptr=" << this << std::endl;
        }
        
        /// The constructor for this base class of Species should not be called directly - only by the concrete subclasses
        /// @param name - a string for the Species name associated with this Species. For example, "G-Actin" or "Arp2/3"
        /// @param type_enum - SType enum, such as SType::Diffusing
        /// @param n - copy number
        /// @param ulim - upper limit for this species' copy number
        Species (const std::string &name, species_copy_t n=0, species_copy_t ulim=max_ulim)  : Component()
        {
            _molecule=SpeciesNamesDB::Instance()->stringToInt(name);
            _rspecies = std::unique_ptr<RSpecies>(new RSpecies(*this, n, ulim));
//            std::cout << "Species (const std::string &name, species_copy_t n=0): Main ctor called, creating ptr=" << this << std::endl;
            
        }
        
        /// Copy constructor
        /// @note The associated RSpecies subpart is not copied, but a new one is created. This means that 
        /// the copied destination Species won't be included in any Reaction interactions of the original 
        /// source Species. The species copy numbered is copied to the target.
        Species (const Species &rhs)  : Component(), _molecule(rhs._molecule) {
            _rspecies = std::unique_ptr<RSpecies>(new RSpecies(*this, rhs.getN(), rhs.getUpperLimitForN()));
//            std::cout << "Species(const Species &rhs): copy constructor called, old ptr=" << &rhs << ", new ptr=" << this << std::endl; 
        }
        
        /// Move constructor - makes it possible to easily add Species to STL containers, such as vector<Species>
        /// One has to be careful with vector<Species> as opposed to vector<Species*>. The latter is "safe" in 
        /// terms of being internally moved around by vector.resize, etc., but vector<Species> would copy Species 
        /// if a move constructor is not available. This will then destruct the original Species (with its 
        /// associated RSpecies), hence, the Reaction network of the Species will be lost. Since this behavior is 
        /// most likely not desired or unacceptable, the internal vector operations should only "move" 
        /// the Species around, without copying. Moving transfers the RSpecies pointer from source to target, 
        /// stealing resources from the source, leaving it for destruction.
        Species (Species &&rhs) noexcept
        : Component(), _molecule(rhs._molecule), _rspecies(std::move(rhs._rspecies)) {
//            std::cout << "Species(Species &&rhs): move constructor called, old ptr=" << &rhs << ", new ptr=" << this << std::endl;
        }
        
        /// Assignment operator
        /// An assignment A = B copies the name of B to A. It destroys the Reaction interactions of A, and resents 
        /// them to a blank value (i.e. A won't be involced in any Reactions). B's Reaction interactions are not copied.
        /// The copy number of B is copied to A.
        Species& operator=(const Species& rhs)  {
//            std::cout << "Species& operator=(const Species& rhs):" << this << ", " << &rhs << std::endl; 
            _molecule = rhs._molecule;
            _rspecies = std::unique_ptr<RSpecies>(new RSpecies(*this, rhs.getN(), rhs.getUpperLimitForN()));
            return *this;
        }
        
        /// Move assignment is needed for the same reasons as move constructor.
        /// @see Species (Species &&rhs)
        Species& operator=(Species&& rhs)  {
//            std::cout << "Species& operator=(Species&& rhs):" << this << ", " << &rhs << std::endl; 
            _molecule = rhs._molecule;
            _rspecies = std::move(rhs._rspecies);
            return *this;
        }
        
        /// Return a reference to RSpecies. Notice that value copying won't be allowed 
        /// because RSpecies is not copyable.
        RSpecies& getRSpecies () {return (*_rspecies.get());}
        
        /// Return a constant reference to RSpecies. 
        const RSpecies& getRSpecies () const {return (*_rspecies.get());}
        
        /// Sets the copy number for this Species. 
        /// @param n should be a non-negative number, but no checking is done in run time
        /// @note The operation does not emit any signals about the copy number change.
        void setN(species_copy_t n) {_rspecies->setN(n);}
        
        /// Return the current copy number of this Species
        species_copy_t getN () const {return _rspecies->getN();}
        
        /// Return the upper limit for the copy number of this Species
        species_copy_t getUpperLimitForN() const {return _rspecies->getUpperLimitForN();}
        
        /// Return this Species' name
        std::string getName() const {return SpeciesNamesDB::Instance()->intToString(_molecule);}
        
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
        
        /// Returns true if two Species objects are equal.
        /// This function would accept derived class of Species, such as SpeciesBulk
        /// Two Species are equal if their SType(s) are equal (i.e. are of the same class),
        /// their names are equal and their copy numbers are equal
        /// Note that specialized versions can be defined in derived classes 
        /// by re-implementing: @see virtual bool is_equal(const Species& b) const 
        friend bool operator ==(const Species& a, const Species& b)
        {
            if (a.getN() != b.getN() or a.getName() != b.getName() or typeid(a) != typeid(b))
                return false;
            return a.is_equal(b); // Invoke virtual is_equal via derived subclass of a (which should be the same as b)
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
        virtual ~Species () noexcept {};
                
        /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Bulk}"
        virtual std::string getFullName() const {return getName() + "{None}";}
        
        /// Default implementation returns true; subclasses should specialize it, if needed
        virtual bool is_equal(const Species& b) const {
            return true;
        }
        
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
        
        /// Return the full name of this Species in a std::string format (e.g. "Arp2/3{Diffusing}"
        virtual std::string getFullName() const {return getName() + "{Diffusing}";}
        
        /// Default destructor
        ~SpeciesDiffusing () noexcept {};
    };
} // end of chem namespace 

/// Print self into an iostream
std::ostream& operator<<(std::ostream& os, const chem::Species& s);

#endif
