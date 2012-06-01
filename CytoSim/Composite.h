//
//  Composite.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/29/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Composite_h
#define CytoSim_Composite_h

#include <string>
#include <vector>

#include "Component.h"
#include "Species.h"

namespace chem {
    
class Composite : public Component {
private:
    std::vector<std::unique_ptr<Species>> _species;
    std::vector<std::unique_ptr<Composite>> _children;
    
public: //should be turned into protected
    void addSpeciesUniq(std::unique_ptr<Species> &&child_species) {
        _species.push_back(std::move(child_species));
        _species.back()->setParent(this);
    }
    
    template<typename T, typename ...Args>
    void addSpecies( Args&& ...args )
    {
        _species.push_back(std::unique_ptr<T>( new T( std::forward<Args>(args)...) ));
        _species.back()->setParent(this);
        //        _species.emplace_back(make_unique(Args...));
    }
    
public:
    Composite() :  Component() {}
    
    virtual std::string getFullName() const {return "Composite";}; 
    
    virtual ~Composite() noexcept {}
    
    virtual void addChild (std::unique_ptr<Composite> &&child) {
        _children.push_back(std::move(child));
        _children.back()->setParent(this);
    }
    
    virtual std::vector<std::unique_ptr<Composite>>& children () {return _children;}

    virtual const std::vector<std::unique_ptr<Composite>>& children () const {return _children;}

    virtual Composite* children (size_t i) {return _children[i].get();}
    
    virtual Species* species(size_t i) {return _species[i].get();}
        
    virtual std::vector<std::unique_ptr<Species>>& species() {
        return _species;
    }
    
    virtual const std::vector<std::unique_ptr<Species>> & species() const {
        return _species;
    }
    
    virtual size_t countSpecies() const {
        size_t res = species().size();
        for(auto &c : _children)
            res+=c->countSpecies();
        return res;
    }

};

} // end of chem
#endif
