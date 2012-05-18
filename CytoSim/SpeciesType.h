//
//  SpeciesType.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/17/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_SpeciesType_h
#define CytoSim_SpeciesType_h

#include <boost/flyweight.hpp>

enum class SType : unsigned char {
    Bulk = 0, ///< Species that have no spatial association (i.e. are "well-mixed") 
    Diffusing = 1, ///< Species that diffuse between cytosolic compartments 
    Membrane = 2, ///< Species that diffuse within a membrane 
    Filament=3, ///< Species that comprise filaments (such as F-Actin)
    Walking=4, ///< Species that can walk ("convectively") on filaments (like Myosin X)
    Motors=5 ///< Species that bind to filaments and generate forces (like Myosin II)
    };
    
    
    /// SpeciesType class represents the type of chemical Species (for example, a Species with name "Arp2/3" and SType Diffusing)
    
    /*! This class describes the type of Species. The description is based on a name (std::string) and SType (Bulk, Diffusing, etc.). 
     *  @note In the future may implement the Logging and Persistence interfaces. 
     *  @note SpeciesType is used by Species with the flyweight pattern, since in a large system millions of unique Species may exist,
     *        but only a handful of SpeciesType objects. For example: one ("F-Actin",Filament) may correspond to millions of individual 
     *        Species.
     */
    class SpeciesType {
    private:
        std::string _name; ///< the descriptive name associated with this Species 
        SType _type; ///< the type of species, such as Bulk, Diffusing, etc.
        static std::vector<std::string> _vec_type_name; ///< this variable is used to help translate strings to enum values for the SType
    public:
        ///Given a species name as a string and its SType, constructs the SpeciesType
        SpeciesType(const std::string &name, SType type) : _name(name), _type(type) {}
        
        ///Given a species name as a string and its type as a string, constructs the SpeciesType
        SpeciesType(const SpeciesType &st) : _name(st._name), _type(st._type) {}
        
        ///The move constructor. May be needed for boost::flyweight?
        SpeciesType(SpeciesType &&st) : _name(std::move(st._name)), _type(st._type) {}

        ///Swap non-member function is needed for the implicit assignment operator
        void swap(SpeciesType &other)
        {   
            std::cout << "SpeciesType::swap(...) was called" << std::endl; 
            _name.swap(other._name);
            std::swap(_type,other._type);
        }
        
        ///The lvalue assignment operator
        SpeciesType& operator=(const SpeciesType &other) 
        {
            _name = other._name;
            _type = other._type;
            return *this;
        } 
        
        ///Move assignment operator. 
        SpeciesType& operator=(SpeciesType&& st){
            _name=std::move(st._name);
            _type=st._type;
            return *this;
        }
                
        ///Equality operator use both species name and type to compare to SpeciesType objects
        bool operator==(const SpeciesType& species_type) const 
        {
            return species_type.getName()==_name && species_type.getType()==_type;
        }
        
        ///Inequality operator use both species name and type to compare to SpeciesType objects
        bool operator!=(const SpeciesType& species_type) const 
        {
            return !(species_type==(*this));
        }
        
        ///Returns the name associated with this SpeciesType as a string
        std::string getName() const {return _name;}
        
        
        ///Returns the SType associated with this SpeciesType 
        SType getType() const {return _type;}
        
        ///Returns a string representing this SpeciesType by concatanating its name and type 
        std::string getTypeAsString () const {return _vec_type_name[static_cast<int>(_type)];}
        
        ///Return if true if this SpeciesType has name and SType given as parameters to this function 
        bool is_of_type(const std::string &name, SType type) const {return name==_name && type==_type;}
        
        ///boost::flyweight needs a hashing function, defined here.
        friend std::size_t hash_value(SpeciesType const& species_type){
            std::size_t seed = 0;
            int type=static_cast<int>(species_type.getType());
            boost::hash_combine(seed,species_type.getName());
            boost::hash_combine(seed,type);
            //        std::cout << species_type.getName()+"[" + species_type.getTypeAsString() << "] hash_value called...\n";
            return seed;
        }
        
        friend std::ostream& operator<<(std::ostream& os, const SpeciesType& st){
            os << st.getName() << "[" << st.getTypeAsString() << "]";
            return os;
        }
        
//        friend istream& operator>>(istream& is, SpeciesType& st) 
//        { 
//            string temp; 
//            is >> temp; 
//            st.set(temp); 
//            return is; 
//        }
    };
#endif
