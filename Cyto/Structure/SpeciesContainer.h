//
//  SpeciesContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 8/30/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef __CytoSim__SpeciesContainer__
#define __CytoSim__SpeciesContainer__

#include <iostream>

#include "common.h"
#include "Species.h"

/// An abstract interface for a container of pointers to Species objects
class SpeciesPtrContainerIFace {
public:
    /// Add species to the container. The memory of species is owned by the container
    virtual Species* addSpecies(Species *species) = 0;

    /// Add species to the container. The memory of species is owned by the container
    virtual Species* addSpeciesUnique(unique_ptr<Species> &&species) = 0;

//        virtual Species* addSpecies(const string &name, species_copy_t copy) = 0;

    /// Remove all species matching "name" from the container. The memories are freed.
    virtual size_t removeSpecies(const string &name) = 0;

    /// Remove species from the container. The memory is freed.
    virtual size_t removeSpecies(Species* species) = 0;
    
    /// Empty the container. All memories are freed.
    virtual void clear() = 0;
    
    /// Return a pointer to Species which has a name matching the argument.
    /// @note The first match is returned.
    virtual Species* findSpeciesByName(const string &name) const = 0;
    
    /// Return a pointer to Species having specified index in the container 
    virtual Species* findSpeciesByIndex (size_t index) const = 0;
    
    /// Return a pointer to Species matching the molecule field of Species
    /// @note The first match is returned.
    virtual Species* findSpeciesByMolecule (int molecule) const = 0;
    
    /// Return a pointer to Species which satisfies the equality operator with s
    /// @note The first match is returned.
    virtual Species* findSimilarSpecies (const Species &s) const = 0;
    
    /// Return true, if all contained Species are different from each other as determined by Species.getMolecule() method.
    virtual bool areAllSpeciesUnique () const = 0;
    
    /// Return the number of Species in the container
    virtual size_t size() const = 0;
    
    /// Print all Species contained by the container
    virtual void printSpecies() const = 0;
};




/// A concrete class implementing the SpeciesPtrContainerIFace, using vector<unique_ptr<Species>> as the container implementation
class SpeciesPtrContainerVector : public  SpeciesPtrContainerIFace {
protected:
    vector<unique_ptr<Species>> _species;///< Species pointers container
public:
    /// Default Constructir
    SpeciesPtrContainerVector() = default;
    
    /// Copying is not allowed
    SpeciesPtrContainerVector(const SpeciesPtrContainerVector &) = delete;
    
    /// The equality operator is not allowed
    SpeciesPtrContainerVector& operator=(SpeciesPtrContainerVector &) = delete;  // no assignment
    
//        friend void swap(SpeciesPtrContainerVector& first, SpeciesPtrContainerVector& second) // nothrow
//        {
//            // enable ADL (not necessary in our case, but good practice)
//            using swap;
//            swap(first._species, second._species);
//        }
    
    
    /// Empty the container. All memories are freed.
    virtual void clear() override {_species.clear();}
    
    /// Add species to the container. The memory of species is owned by the container
    virtual Species* addSpecies(Species *species) {
        _species.emplace_back(unique_ptr<Species>(species));
        return _species.back().get();
    }
    
    /// Add species to the container. The memory of species is owned by the container.
    /// @param species is a unique_ptr<Species> object, which needs to be a rvalue (e.g.
    /// move(...) was used to make it so.
    virtual Species* addSpeciesUnique (unique_ptr<Species> &&species) override {
        _species.emplace_back(move(species));
        return _species.back().get();
    }
    
    /// Add species of type T to the container, forwarding Args to the corresponding Species constructor 
    template<typename T, typename ...Args>
    Species* addSpecies( Args&& ...args )
    {
        _species.emplace_back(unique_ptr<Species>( new T( forward<Args>(args)...) ));
        //        _species.emplace_back(make_unique(Args...));
        return _species.back().get();
    }
    
    /// Remove all species matching "name" from the container. The memories are freed.
    virtual size_t removeSpecies(const string &name) override {
        size_t counter=0;
        while(true){
            auto child_iter = find_if(_species.begin(),_species.end(),
                                           [&name](const unique_ptr<Species> &element)
                                           {
                                               return element->getName()==name ? true : false;
                                           });
            if(child_iter!=_species.end()){
                _species.erase(child_iter);
                ++counter;
            }
            else
                break;
        }
        return counter;
    }
    
    /// Remove all species from the container. The memory is freed.
    virtual size_t removeSpecies (Species* species) override {
        size_t counter=0;
        while(true) {
            auto child_iter = find_if(_species.begin(),_species.end(),
                                           [&species](const unique_ptr<Species> &element)
                                           {
                                               return element.get()==species ? true : false;
                                           });
            if(child_iter!=_species.end()){
                _species.erase(child_iter);
                ++counter;
            }
            else
                break;
        }
        return counter;
    }
    
    /// Return a pointer to Species which has a name matching the argument. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual Species* findSpeciesByName(const string &name) const override {
        auto child_iter = find_if(_species.begin(),_species.end(),
                                       [&name](const unique_ptr<Species> &element)
                                       {
                                           return element->getName()==name ? true : false;
                                       });
        if(child_iter!=_species.end())
            return child_iter->get();
        else
            return nullptr;
    }
    
    /// Return a pointer to Species matching the molecule field of Species. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual Species* findSpeciesByMolecule (int molecule) const override {
        auto child_iter = find_if(_species.begin(),_species.end(),
                                       [molecule](const unique_ptr<Species> &element)
                                       {
                                           return element->getMolecule()==molecule ? true : false;
                                       });
        if(child_iter!=_species.end())
            return child_iter->get();
        else
            return nullptr;
    }
    
    /// Return a pointer to Species having specified index in the container.
    /// @note Array bounds are not checked, so index must be valid.
    virtual Species* findSpeciesByIndex (size_t index) const  override {
        return _species[index].get();
    }
    
    /// Return a pointer to Species which satisfies the equality operator with s. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual Species* findSimilarSpecies (const Species &s) const  override {
        auto it = find_if(_species.begin(),_species.end(),
                               [&s](const unique_ptr<Species> &element)
                               {return s==(*element);});
        if(it!=_species.end())
            return it->get();
        return nullptr;
    }
    
    
    /// Return a reference to the underlying vector<unique_ptr<Species>> container.
    vector<unique_ptr<Species>>& species() {return _species;}
    
    /// Return a const reference to the underlying vector<unique_ptr<Species>> container.
    const vector<unique_ptr<Species>>& species() const {return _species;}
    
    /// Returns true, if all contained Species are different from each other as determined by Species.getMolecule() method.
    virtual bool areAllSpeciesUnique () const override {
        vector<int> molecs;
        transform(_species.cbegin(),_species.cend(), back_inserter(molecs),
                       [](const unique_ptr<Species> &us)
                            {return us->getMolecule();});
        sort(molecs.begin(), molecs.end());
        auto vit = adjacent_find(molecs.begin(), molecs.end());
        if(vit==molecs.end())
            return true;
        else
            return false;
    }
    
    /// Return the number of Species in the container
    virtual size_t size() const override {return _species.size();}
    
    /// Print all Species contained by the container
    virtual void printSpecies() const override {
        for(auto &s : _species)
            cout << (*s.get()) << endl;
    }

};

/// An abstract interface for a container of Species objects
class SpeciesContainerIFace {
public:
    /// Remove all Species with the specified index from the container
    virtual size_t removeSpecies(size_t index) = 0;
    
    /// Remove all Species with the specified name from the container
    virtual size_t removeSpecies(const string &name) = 0;
    
    /// Return a reference to the first Species having the specified name.
    /// @note Throw an exception if the Species is not found.
    virtual Species& findSpecies(const string &name) = 0;
    
    /// Find the index of Species with the specified name in the container.
    /// @note The first match is returned. An exception is thrown if the Species is not found.
    virtual size_t findSpeciesIndex(const string &name) const = 0;
    
    /// Return a reference to Species at the position index in the container
    virtual Species& findSpecies (size_t index) = 0;
    
    /// Return a reference to Species which satisfies the equality operator with s. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual Species& findSimilarSpecies (const Species &s) = 0;
    
    /// Return a reference to Species matching the molecule field of Species. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual Species& findSpeciesByMolecule (int molecule) = 0;

    /// Return true, if all contained Species are different from each other as determined by Species.getMolecule() method.
    virtual bool areAllSpeciesUnique () const = 0;
    
    /// Return the number of Species in the container
    virtual size_t size() const = 0;
    
    /// Print all Species in the container
    virtual void printSpecies() const {}
};




/// A concrete class implementing the SpeciesContainerIFace, using vector<SpeciesSpecific> as the container implementation.

/*! Because in the class the Species are held by value, and not as pointers, the container must be homogeneous,
 *  i.e. consist of Species of the same derived type. This is the SpeciesSpecific template paramter, which 
 *  for example can be SpeciesDiffuse.
 *
 *  When adding Species to this class, the return value of the addSpecies() method could be remembered, so 
 *  to access that Species in the future by index.
 */
template <class SpeciesSpecific>
class SpeciesContainerVector : public  SpeciesContainerIFace {
protected:
    vector<SpeciesSpecific> _species; ///< The container of Species of type SpeciesSpecific
public:
    /// Add species of type SpeciesSpecific to the container, forwarding Args to the corresponding Species constructor.
    /// @return the index of the added Species in the container
    template<typename ...Args>
    size_t addSpecies(Args&& ...args){
//            cout << "SpeciesContainerVector::addSpecies()..." << endl;
        _species.emplace_back(forward<Args>(args)...);
        return _species.size()-1;
    }
    
    /// Remove all Species with the specified index from the container
    virtual size_t removeSpecies(size_t index) override {
        assert(index<size() && "SpeciesContainerVector::removeSpecies(): Invalid index.");
        _species.erase(_species.begin()+index);
        return 1;
    }
    
    /// Remove all Species with the specified name from the container
    virtual size_t removeSpecies(const string &name) override {
        size_t counter = 0;
        while(true) {
            auto child_iter = find_if(_species.begin(),_species.end(),
                                           [&name](const Species &element)
                                           {
                                               return element.getName()==name ? true : false;
                                           });
            if(child_iter!=_species.end()){
                _species.erase(child_iter);
                ++counter;
            }
            else
                break;
        }
        return counter;
    }
    
    /// Return a reference to the first Species having the specified name.
    /// @note Throw an exception if the Species is not found.
    virtual SpeciesSpecific& findSpecies(const string &name) override {
        auto child_iter = find_if(_species.begin(),_species.end(),
                                       [&name](const Species &element)
                                       {
                                           return element.getName()==name ? true : false;
                                       });
        if(child_iter!=_species.end())
            return (*child_iter);
        else
            throw out_of_range("Species::findSpecies(): The name was not found");
    }
    
    /// Find the index of Species with the specified name in the container.
    /// @note The first match is returned. An exception is thrown if the Species is not found.
    virtual size_t findSpeciesIndex(const string &name) const override {
        size_t index = 0;
        for(auto &s : _species){
            if(s.getName()==name)
                return index;
            else
                ++index;
        }
        throw out_of_range("Species::findSpecies(): The name was not found");
    }

    
    /// Return a reference to Species at the position index in the container.
    /// @note No bound checking is done, hence, index must be valid.
    virtual SpeciesSpecific& findSpecies (size_t index) override {
        return _species[index];
    }
    
    /// Returns true, if all contained Species are different from each other as determined by Species.getMolecule() method.
    virtual bool areAllSpeciesUnique () const override {
        vector<int> molecs;
        transform(_species.cbegin(),_species.cend(), back_inserter(molecs),
                       [](const Species &s)
                       {return s.getMolecule();});
        sort(molecs.begin(), molecs.end());
        auto vit = adjacent_find(molecs.begin(), molecs.end());
        if(vit==molecs.end())
            return true;
        else
            return false;
    }
    
    /// Return a reference to Species which satisfies the equality operator with s. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual SpeciesSpecific& findSimilarSpecies (const Species &s) override {
        auto it = find_if(_species.begin(),_species.end(),
                               [&s](const Species &element)
                               {return s==element;});
        if(it!=_species.end())
            return *it;
        else
            throw out_of_range("SpeciesContainerVector::findSimilarSpecies(): The name was not found");
    }
    
    /// Return a reference to Species matching the molecule field of Species. Otherwise, return a nullptr.
    /// @note The first match is returned.
    virtual SpeciesSpecific& findSpeciesByMolecule (int molecule) override {
        auto child_iter = find_if(_species.begin(),_species.end(),
                                       [molecule](const Species &element)
                                       {
                                           return element.getMolecule()==molecule ? true : false;
                                       });
        if(child_iter!=_species.end())
            return *child_iter;
        else
            throw out_of_range("SpeciesContainerVector::findSpeciesByMolecule(): The molecule was not found");
    }
    
    /// Return the number of Species in the container
    virtual size_t size() const override {return _species.size();}
    
    /// Print all Species in the container
    virtual void printSpecies() const override {
        for(auto &s : _species)
            cout << s << endl;
    }
    
    /// Return a reference to the underlying vector<unique_ptr<Species>> container.
    vector<SpeciesSpecific>& species() {return _species;}
    
    /// Return a const reference to the underlying vector<unique_ptr<Species>> container.
    const vector<SpeciesSpecific>& species() const {return _species;}

};

#endif /* defined(__CytoSim__SpeciesContainer__) */
