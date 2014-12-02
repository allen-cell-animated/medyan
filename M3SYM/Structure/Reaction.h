
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for installation and documentation:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Reaction_h
#define M3SYM_Reaction_h

#include "common.h"

#include "ReactionBase.h"
#include "ChemSim.h"

/// Represents a concrete chemical reaction, such as A + B -> C, where M is the number of reactants and N is the number of products.

/*! Reaction<M,N> encodes a chemical reaction between M reactants and N products. It follows the ReactionBase interface, where 
 *  many methods are defined. Most of the methods defined in Reaction<M,N> are specific implementations of virtual functions 
 *  declared in ReactionBase. Hence, to a very large degree, Reaction<M,N> is an implementation class, while ReactionBase 
 *  provides the public interface for reaction objects. The documentation of the latter should be mainly consulted when 
 *  working with reactions.
 */
template <unsigned short M, unsigned short N>
    class Reaction : public ReactionBase {
    private:
        array<RSpecies*, M+N> _rspecies;///< An array of RSpecies objects (reactants followed by products)
    public:
        /// The main constructor:
        /// @param species - are reactants and products put together into a single list (starting from reactants)
        /// @param rate - the rate constant for this ReactionBase
        Reaction(initializer_list<Species*> species, float rate = 0.0, bool isProtoCompartment = false) : ReactionBase(rate, isProtoCompartment)
        {
            initializeSpecies(species);
        }
        
        /// The main constructor:
        /// @param it_begin - an iterator to the beginning of an RSpecies* container
        /// @param it_end - an iterator to the end of an RSpecies* container
        /// @param rate - the rate constant for this ReactionBase
        template <typename InputContainer>
        Reaction(const InputContainer &species, float rate = 0.0, bool isProtoCompartment = false) : ReactionBase(rate, isProtoCompartment)
        {
            //            cout << "Reaction<M,N>(const vector<Species*> &species, float rate) called..." << endl;
            initializeSpecies(species);
        }
        
        /// no copying (including all derived classes)
        Reaction (const Reaction &rb) = delete;
        
        /// no assignment (including all derived classes)
        Reaction& operator=(Reaction &rb) = delete;

#ifdef BOOST_MEM_POOL
        /// Advanced memory management
        void* operator new(size_t size);
        
        void operator delete(void* ptr) noexcept;
#endif
        /// Destructor
        /// Tell Rspecies to remove this Reaction from its internal lists of reactions
        /// @note noexcept is important here. Otherwise, gcc flags the constructor as potentially throwing,
        /// which in turn disables move operations by the STL containers. This behaviour is a gcc bug
        /// (as of gcc 4.703), and will presumbaly be fixed in the future.
        virtual ~Reaction() noexcept
        {
            for(auto i=0U; i<M; ++i) _rspecies[i]->removeAsReactant(this);
            for(auto i=M; i<(M+N); ++i) _rspecies[i]->removeAsProduct(this);
        }
        
        /// Returns a pointer to the first element of array<RSpecies*, M+N>
        /// The type of the pointer is RSpecies**. In conjunction with getM() and getN(),
        /// this pointer can be used to iterate over RSpecies associated with this reaction.
        inline virtual RSpecies** rspecies() override {return &_rspecies[0];}
        
        /// Return a list of reactions which rates would be affected if this
        /// reaction were to be executed.
        virtual vector<ReactionBase*> getAffectedReactions() override
        {
            unordered_set<ReactionBase*> rxns;
            for(auto s : _rspecies){
                
                for(auto it = s->beginReactantReactions(); it != s->endReactantReactions(); it++) {
                    ReactionBase* r = (*it);
                    rxns.insert(r);
                }
            }
            rxns.erase(this);
            return vector<ReactionBase*>(rxns.begin(),rxns.end());
        }
        
    protected:
        /// An implementation method used by the constructor.
        template <typename InputContainer>
        void initializeSpecies(const InputContainer &species)
        {
            assert(species.size()==(M+N) && "Reaction<M,N> Ctor: The species number does not match the template M+N");
            transform(species.begin(),species.end(),_rspecies.begin(),
                      [](Species *s){return &s->getRSpecies();});
            
            if(!_isProtoCompartment) {
                //add dependents
                for(auto &d : getAffectedReactions()) { if(!d->isPassivated()) _dependents.push_back(d); }
                
                for(auto i=0U; i<M; ++i) _rspecies[i]->addAsReactant(this);
                for(auto i=M; i<(M+N); ++i) _rspecies[i]->addAsProduct(this);
            }

        }
        
        /// Implementation of getM()
        inline virtual unsigned short getMImpl() const override {return M;}

        /// Implementation of getN()
        inline virtual unsigned short getNImpl() const override {return N;}

        /// Implementation of size()
        inline virtual unsigned short sizeImpl() const override {return M+N;}

        /// Implementation of the operator==(...)
        virtual bool is_equal(const ReactionBase& other) const override
        {
            const Reaction<M,N> *a = this;
            const Reaction<M,N> *b = static_cast<const Reaction<M,N>*>(&other);
            auto it_pair = mismatch(a->_rspecies.begin(),a->_rspecies.end(),b->_rspecies.begin(),
                                         [](RSpecies* A, RSpecies* B){return A->getSpecies()==B->getSpecies();});
            if(it_pair.first==a->_rspecies.end())
                return true;
            return false;
        }
        
        /// Implementation of getProductOfReactants()
        inline virtual int getProductOfReactantsImpl() const override
        {
            int prod = 1;
            for(auto i=0U; i<M; ++i)
                prod*=_rspecies[i]->getN();
            return prod;
            
        }

        /// Implementation of computePropensity()
        inline virtual float computePropensityImpl() const override
        {
#ifdef TRACK_UPPER_COPY_N
            if(this->Reaction<M,N>::getProductOfProductsImpl()==0){
                return float(0.0);
            }
#endif
            return _rate*Reaction<M,N>::getProductOfReactantsImpl();
        }
        
        /// Implementation of getProductOfProducts()
        inline virtual int getProductOfProductsImpl() const override
        {
#ifdef TRACK_UPPER_COPY_N
            int prod = 1;
            for(auto i=M; i<(M+N); ++i){
                prod*=_rspecies[i]->getN()-_rspecies[i]->getUpperLimitForN();
            }
            return prod;
#else
            return 1;
#endif
        }
    
        /// Implementation of containsSpecies()
        inline virtual bool containsSpeciesImpl(Species *s) const override
        {
            auto it = find_if(_rspecies.begin(), _rspecies.end(),
                                   [s](const RSpecies *rs){return (&rs->getSpecies())==s;});
            return it!=_rspecies.end();
            
        }
    
        /// Implementation of makeStep()
        inline virtual void makeStepImpl() override
        {
            for(auto i=0U; i<M; ++i)
                _rspecies[i]->down();
            for(auto i=M; i<(M+N); ++i)
                _rspecies[i]->up();
        }

        /// Implementation of activateReactionUnconditional()
        virtual void activateReactionUnconditionalImpl() override;

        /// Implementation of passivateReaction()
        virtual void passivateReactionImpl() override;
        
        /// Print information about this reaction to ostream
        virtual void printToStream(ostream& os) const override
        {
            unsigned char i=0;
            auto sit = _rspecies.cbegin();
            auto send = sit+M;
            for (; sit!=send; ++sit)
            {
                os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
                if(i<M-1)
                os << " + ";
                ++i;
            }
            os << " ---> ";
            i=0;
            for (auto sit = send; sit!=_rspecies.cend(); ++sit)
            {
                os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
                if(i<(N-1))
                    os << " + ";
                ++i;
            }
            os << ", " << "curr_rate = " << getRate() << ", a=" << computePropensity() << ", ReactionBase ptr=" << this << "\n";
        }
        
#ifdef RSPECIES_SIGNALING
        /// This method may be called after a makeStep() took place. It iterates over all
        /// reactants and products, and call the emitSignal(delta) of the corresponding
        /// RSpecies objects, where delta=-1 for reactants and delta=+1 for products.
        virtual void broadcastRSpeciesSignals() 
        {
            for(auto i=0U; i<M; ++i)
            {
                if(_rspecies[i]->isSignaling())
                        _rspecies[i]->emitSignal(-1);

            }
            for(auto i=M; i<(M+N); ++i)
            {
                if(_rspecies[i]->isSignaling())
                    _rspecies[i]->emitSignal(+1);
            }
        }
#endif
        
        /// Implementation of  clone()
        virtual Reaction<M,N>* cloneImpl(const SpeciesPtrContainerVector &spcv) override;
    };
    
    /// Partial template speciatialization for Reaction<1,1> to gain efficiency
    template <> inline float Reaction<1,1>::computePropensityImpl() const
    {
#ifdef TRACK_UPPER_COPY_N
        if(_rspecies[1]->getN()>=_rspecies[1]->getUpperLimitForN())
            return 0;
#endif
        return _rate*_rspecies[0]->getN();
    }
    

    
    /// Partial template speciatialization for Reaction<1,1> to gain efficiency
    template <> inline void Reaction<1,1>::makeStepImpl()
    {
        _rspecies[0]->down();
        _rspecies[1]->up();
    }
    
#ifdef RSPECIES_SIGNALING
    /// Partial template speciatialization for Reaction<1,1> to gain efficiency
    template <>  inline void Reaction<1,1>::broadcastRSpeciesSignals()
    {
        if(_rspecies[0]->isSignaling())
            _rspecies[0]->emitSignal(-1);
        
        if(_rspecies[1]->isSignaling())
            _rspecies[1]->emitSignal(+1);
    }
#endif

#endif


