#ifndef MEDYAN_Chemistry_ReactionDy_hpp
#define MEDYAN_Chemistry_ReactionDy_hpp

#include <algorithm>
#include <stdexcept>

#ifdef BOOST_MEM_POOL
    #include <boost/pool/pool.hpp>
    #include <boost/pool/pool_alloc.hpp>
#endif

#include "Chemistry/ChemRNode.h"
#include "Chemistry/ReactionBase.h"
#include "Util/Io/Log.hpp"

class ReactionDy : public ReactionBase {

private:
    std::vector< RSpecies* > rspecies_;
    unsigned short numReactants_ = 0;
    unsigned short numProducts_  = 0;

    void registerInRSpecies_() {
        for(unsigned short i = 0; i < numReactants_; ++i) {
            rspecies_[i]->addAsReactant(this);
        }
        for(unsigned short i = 0; i < numProducts_; ++i) {
            rspecies_[i + numReactants_]->addAsProduct(this);
        }
    }
    void unregisterInRSpecies_() {
        for(unsigned short i = 0; i < numReactants_; ++i) {
            rspecies_[i]->removeAsReactant(this);
        }
        for(unsigned short i = 0; i < numProducts_; ++i) {
            rspecies_[i + numReactants_]->removeAsProduct(this);
        }
    }

public:

    template< typename InputContainer >
    ReactionDy(
        const InputContainer& reactantSpecies,
        const InputContainer& productSpecies,
        float rate,
        float volumeFrac = 1.0,
        int rateVolumeDepExp = 0
    ) : ReactionBase(rate, false, volumeFrac, rateVolumeDepExp) {
        initializeSpecies(reactantSpecies, productSpecies);
    }
        
    ReactionDy(const ReactionDy &) = delete;
        
#ifdef BOOST_MEM_POOL
    /// Advanced memory management
    void* operator new(size_t size) {
        void *ptr = boost::fast_pool_allocator< ReactionDy >::allocate();
        return ptr;
    }

    void operator delete(void* ptr) noexcept {
        boost::fast_pool_allocator< ReactionDy >::deallocate((ReactionDy*)ptr);
    }
#endif

    virtual ~ReactionDy() noexcept override {
        unregisterInRSpecies_();
    }

    ReactionDy& operator=(const ReactionDy &) = delete;


    /// Return a list of reactions which rates would be affected if this
    /// reaction were to be executed.
    virtual vector<ReactionBase*> getAffectedReactions() override {
        unordered_set<ReactionBase*> rxns;

        for(int i = 0; i < rspecies_.size(); i++) {

            for(auto r : rspecies_[i]->reactantReactions()) {
                rxns.insert(r);
            }
        }
        rxns.erase(this);
        return vector<ReactionBase*>(rxns.begin(),rxns.end());
    }

    virtual void updatePropensityImpl() override {
        if(_rnode && !_passivated) _rnode->activateReaction();
    }

    template <typename InputContainer>
    void initializeSpecies(
        const InputContainer& reactantSpecies,
        const InputContainer& productSpecies
    ) {
        using namespace std;

        // Unregister this in all RSpecies
        unregisterInRSpecies_();

        numReactants_ = reactantSpecies.size();
        numProducts_  = productSpecies.size();

        rspecies_.resize(numReactants_ + numProducts_);
        for(unsigned short i = 0; i < numReactants_; ++i) {
            rspecies_[i] = &reactantSpecies[i]->getRSpecies();
        }
        for(unsigned short i = 0; i < numProducts_; ++i) {
            rspecies_[i + numReactants_] = &productSpecies[i]->getRSpecies();
        }
        
        //add dependents
        for(auto d : getAffectedReactions())
            if(!d->isPassivated()) _dependents.insert(d);

        // Register this in all RSpecies
        registerInRSpecies_();
    }
        
    virtual unsigned short getMImpl() const override { return numReactants_; }
    virtual unsigned short getNImpl() const override { return numProducts_; }
    virtual unsigned short sizeImpl() const override { return numReactants_ + numProducts_; }

    /// Implementation of getProductOfReactants()
    virtual floatingpoint getProductOfReactantsImpl() const override
    {
        floatingpoint prod = 1;
        for(unsigned short i = 0; i < numReactants_; ++i)
            prod *= rspecies_[i]->getN();
        return prod;
    }

    /// Implementation of computePropensity()
    virtual floatingpoint computePropensityImpl() const override
    {
        if(isPassivated()) return 0.0;
#ifdef TRACK_UPPER_COPY_N
        if(areEqual(getProductOfProductsImpl(), 0.0)){
            return 0.0;
        }
#endif
        return _rate * getProductOfReactantsImpl();
    }
        
    /// Implementation of getProductOfProducts()
    virtual floatingpoint getProductOfProductsImpl() const override
    {
#ifdef TRACK_UPPER_COPY_N
        floatingpoint prod = 1;
        for(unsigned short i = 0; i < numProducts_; ++i){
            auto& rs = *rspecies_[i + numReactants_];
            prod *= rs.getN() - rs.getUpperLimitForN();
        }
        return prod;
#else
        return 1.0;
#endif
    }


    /// Implementation of makeStep()
    virtual void makeStepImpl() override
    {
        for(unsigned short i = 0; i < numReactants_; ++i) rspecies_[i]->down();
        for(unsigned short i = 0; i < numProducts_; ++i) rspecies_[i + numReactants_]->up();
    }

    /// Implementation of activateReactionUnconditional()
    virtual void activateReactionUnconditionalImpl() override {
#ifdef TRACK_DEPENDENTS
        for(unsigned short i; i < numReactants_; ++i)
        {
            auto& rs = *rspecies_[i];
            for(auto r : rs.reactantReactions()) {
                if(this != r) r->registerNewDependent(this);
            }
            for(auto r : rs.productReactions()) {
                if(this != r) r->registerNewDependent(this);
            }
        }
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
        _passivated = false;
#endif

        if(_rnode) _rnode->activateReaction();
    }

    /// Implementation of passivateReaction()
    virtual void passivateReactionImpl() override {
        if(isPassivated()) return;
#ifdef TRACK_DEPENDENTS
        for(unsigned short i = 0; i < numReactants_; ++i)
        {
            auto& rs = *rspecies_[i];
            for(auto r : rs.reactantReactions()) {
                r->unregisterDependent(this);
            }
            for(auto r : rs.productReactions()) {
                r->unregisterDependent(this);
            }
        }
#endif
#if defined TRACK_ZERO_COPY_N || defined TRACK_UPPER_COPY_N
        _passivated = true;
#endif
        if(_rnode) _rnode->passivateReaction();
    }


    /// Print information about this reaction to ostream
    virtual void printToStream(ostream& os) const override
    {
        unsigned char i=0;
        auto sit = rspecies_.cbegin();
        auto send = sit + numReactants_;
        for (; sit!=send; ++sit)
        {
            os << (*sit)->getFullName() << "{" << (*sit)->getN() << "}";
            if(i < numReactants_ - 1)
                os << " + ";
            ++i;
        }
        os << " ---> ";
        i=0;
        for (auto sit2 = send; sit2 != rspecies_.cend(); ++sit2)
        {
            os << (*sit2)->getFullName() << "{" << (*sit2)->getN() << "}";
            if(i< numProducts_ - 1)
                os << " + ";
            ++i;
        }
        os << ", " << "curr_rate = " << getRate() << ", a="
            << computePropensity() << ", ReactionBase ptr=" << this <<"\n";
    }

    virtual bool updateDependencies() override { return true; }

    // Only for the purpose of overriding
    // If possible, these functions should be removed.
    //---------------------------------

    /// Returns a pointer to the first element of array<RSpecies*, M+N>
    /// The type of the pointer is RSpecies**. In conjunction with getM() and
    /// getN(), this pointer can be used to iterate over RSpecies associated with
    /// this reaction.
    virtual RSpecies** rspecies() override { return &rspecies_[0]; }

    virtual bool is_equal(const ReactionBase& other) const override {
        LOG(ERROR) << "Virtual equality comparison of ReactionDy should not be used";
        throw std::runtime_error("Invalid comparsion of ReactionDy");
    }

    /// Implementation of containsSpecies()
    virtual bool containsSpeciesImpl(Species *s) const override
    {
        auto it = find_if(rspecies_.begin(), rspecies_.end(),
            [s](const RSpecies *rs){ return (&rs->getSpecies()) == s; });
        return it != rspecies_.end();
    }

    /// Implementation of  clone()
    virtual ReactionDy* cloneImpl(const SpeciesPtrContainerVector &spcv) override {
        LOG(ERROR) << "Virtual clone function in ReactionDy is not allowed.";
        throw std::runtime_error("Invalid clone in ReactionDy");
    }

};

#endif
