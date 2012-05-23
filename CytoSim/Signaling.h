//
//  Signaling.h
//  CytoSim
//
//  Created by Garegin Papoian on 5/13/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Signaling_h
#define CytoSim_Signaling_h

#include <memory>
#include <unordered_map>
#include <functional>
#include <cassert>

#include <boost/signals2/signal.hpp>
#include <boost/signals2/connection.hpp>
#include <boost/signals2/shared_connection_block.hpp>
#include "utility.h"

namespace chem {

class RSpecies;
class Reaction;

/// This is a RSpecies signal object that will be called by ChemSignal, usually when requested by the 
/// reaction simulation algorithm
typedef boost::signals2::signal<void (RSpecies *, int)> RSpeciesCopyNChangedSignal;

/// This is a Reaction signal object that will be called by ChemSignal, usually when requested by the 
/// reaction simulation algorithm
typedef boost::signals2::signal<void (Reaction *)> ReactionEventSignal;

/// ChemSignal manages callbacks for RSpecies and Reactions objects. 

/*!  One ChemSignal should be instantiated per chemical network. RSpecies and Reactions that need to be monitored 
 *   can be made to signal, based on the boost:signal2 library. Multiple receiving slots can be connected to signals. 
 *
 *   Slots can be temporarily blocked or completely disconnected. Signals can be destroyed.
 *   Here is an example of usage:
 *   @code 
 ...
 Species* A1 = ...
 ...
 ChemSignal sm;
 A1->makeSignaling(sm);
 PrintRSpecies ps;
 ...
 std::function<void (RSpecies *, int)> psf(ps);
 boost::signals2::shared_connection_block conn_a1(sm.connect(A1, [](RSpecies *s, int i){s->printSelf();}), true);
 // or   sm.connect(A1, print_RSpecies);
 // or  sm.connect(A1, PrintRSpecies());
 // or   boost::signals2::shared_connection_block conn_a1(sm.connect(A1,psf), false);
...
// run simulation for some time
 ...
 conn_a1.block(); // this temporarily stops slots from being called.
...
 conn_a1.disconnect(); // this permanently disconnects the slot
...
 A1->stopSignaling(sm); // this destroys the Signal itself (with associated multiple slots)
...
 @endcode
 */
class ChemSignal {
private:
    std::unordered_map<RSpecies *, std::unique_ptr<RSpeciesCopyNChangedSignal>> _map_RSpecies_signal;///< Keep track of signals corresponding to various RSpecies
    std::unordered_map<Reaction *, std::unique_ptr<ReactionEventSignal>> _map_reaction_signal;///< Keep track of signals corresponding to various Reactions
private:
    /// Search the map for a signal corresponding to parameter, Reaction *r. Throw std::out_of_range exception if not found.
    std::unordered_map<Reaction *, std::unique_ptr<ReactionEventSignal>>::const_iterator 
    findSignal(Reaction *r) const {
        auto sig_it = _map_reaction_signal.find(r);
        if(sig_it==_map_reaction_signal.end())
            throw std::out_of_range("ChemSignal::broadcastReactionSignal(...) key error...");
        return sig_it;
    }
    
    /// Search the map for a signal corresponding to parameter, RSpecies *s. Throw std::out_of_range exception if not found.
    std::unordered_map<RSpecies *, std::unique_ptr<RSpeciesCopyNChangedSignal>>::const_iterator 
    findSignal(RSpecies *s) const {
        auto sig_it = _map_RSpecies_signal.find(s);
        if(sig_it==_map_RSpecies_signal.end())
            throw std::out_of_range("ChemSignal::broadcastRSpeciesSignal(...) key error...");
        return sig_it;
    }
    
    /// Assert that a signal corresponding to Reaction *r does not exist or throw std::runtime_error.
    void assertSignalDoesNotExist (Reaction *r){
        auto sig_it = _map_reaction_signal.find(r);
        if(sig_it!=_map_reaction_signal.end())
            throw std::runtime_error("ChemSignal::addSignalingReaction(...) Reaction already present in the map...");
    }
    
    /// Assert that a signal corresponding to RSpecies *s does not exist or throw std::runtime_error.
    void assertSignalDoesNotExist (RSpecies *s){
        auto sig_it = _map_RSpecies_signal.find(s);
        if(sig_it!=_map_RSpecies_signal.end())
            throw std::runtime_error("ChemSignal::addSignalingRSpecies(...) RSpecies already present in the map...");
    }

public:
    /// Broadcasts signal corresponding to Reaction *r. If the Reaction is not found, std::out_of_range exception will be thrown.
    /// This method is usally only called by the Gillespie-like algorithm, and not the outside code.
    void emitReactionSignal(Reaction *r) const {
        auto sig_it = findSignal(r);
        (*sig_it->second)(r);
    }
    
    /// Broadcasts signal corresponding to RSpecies *s. If the Specis is not found, std::out_of_range exception will be thrown.
    /// This method is usally only called by the Gillespie-like algorithm, and not the outside code.
    void emitRSpeciesSignal(RSpecies *s, int delta) const {
        auto sig_it = findSignal(s);
        (*sig_it->second)(s,delta);
    }
    
    /// Inserts a signal corresponding to RSpecies *s into the map. Should only be called by the RSpecies class.
    /// @note other classes should instead call void RSpecies::makeSignaling(ChemSignal &sm);
    void addSignalingRSpecies(RSpecies *s){
        assertSignalDoesNotExist(s);
        _map_RSpecies_signal.emplace(s,make_unique<RSpeciesCopyNChangedSignal>());
    }

    
    /// Inserts a signal corresponding to Reaction *r into the map. Should only be called by the Reaction class.
    /// @note other classes should instead call void Reaction::makeSignaling(ChemSignal &sm);
    void addSignalingReaction (Reaction *r){
        assertSignalDoesNotExist(r);
        _map_reaction_signal.emplace(r,make_unique<ReactionEventSignal>());
    }
    
    /// Connect the callback, react_callback to a signal corresponding to Reaction *r.
    /// @param Reaction *r - the signal will correspond to this Reaction
    /// @param std::function<void (Reaction *)> const &react_callback - a function object to be called (a slot)
    /// @param int priority - lower priority slots will be called first. Default is 5 Do not use priorities 1 and 2 
    ///                       unless absolutely essential.
    boost::signals2::connection connect(Reaction *r, std::function<void (Reaction *)> const &react_callback, int priority=5) {
        auto sig_it = findSignal(r);
        return sig_it->second->connect(priority, react_callback);
    }

    /// Connect the callback, RSpecies_callback to a signal corresponding to RSpecies *s.
    /// @param RSpecies *s - the signal will correspond to this RSpecies
    /// @param std::function<void (RSpecies *, int)> const &RSpecies_callback - a function object to be called (a slot). 
    ///        int here corresponds to delta, the change in copy number of the RSpecies (for which the signal was emitted)
    /// @param int priority - lower priority slots will be called first. Default is 5 Do not use priorities 1 and 2 
    ///                       unless absolutely essential.
    boost::signals2::connection connect(Species *s, std::function<void (RSpecies *, int)> const &RSpecies_callback, int priority=10) {
        RSpecies *rs = &s->getRSpecies();
        auto sig_it = findSignal(rs);
        return sig_it->second->connect(priority, RSpecies_callback);
    }
    
    /// Destroy signal corresponding to Reaction r. Should only be used by the Reaction class.
    /// @note   Other classes should instead call void Reaction::stopSignaling (ChemSignal &sm);
    void disconnect(Reaction *r){
        auto sig_it = findSignal(r);
        _map_reaction_signal.erase(sig_it);
    }

    /// Destroy signal corresponding to Reaction r. Should only be used by the RSpecies class.
    /// @note   Other classes should instead call void RSpecies::stopSignaling (ChemSignal &sm);
    void disconnect(RSpecies *s){
        auto sig_it = findSignal(s);
        _map_RSpecies_signal.erase(sig_it);
    }
};

} //end of namespace 
#endif
