//
//  Reaction.h
//  CytoSim-Experimenting
//
//  Created by Garegin Papoian on 4/19/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_Experimenting_Reaction_h
#define CytoSim_Experimenting_Reaction_h

#include <iostream>
#include <algorithm>
#include <utility>
#include <array>
#include "common.h"
#include "Species.h"


class ReactionBase {
protected:
    float _forward_rate = 0.0;
    float _backward_rate = 0.0;
    reaction_num_t _id;
    static reaction_num_t _max_id;
public:
    // Constructors
    ReactionBase (float f_rate, float b_rate) : _forward_rate(f_rate), _backward_rate(b_rate), _id(_max_id++){};
    ReactionBase (const ReactionBase &r) = delete; // no copying (including all derived classes)
    ReactionBase& operator=(ReactionBase &r) = delete;  // no assignment (including all derived classes)
    // Setters 
    void setForwardRate(float rate){_forward_rate=rate;}
    void setBackwardRate(float rate){_backward_rate=rate;}
    // Accessors 
    float getForwardRate(float rate){return _forward_rate;}
    float getBackwardRate(float rate){return _backward_rate;}
    reaction_num_t reactionID() const {return _id;}
    //Destructor
    virtual ~ReactionBase() {std::cout << "Reaction, ID=" << _id << " is destroyed" << std::endl;};
    //Fire Reactions
    void doStep(bool forward) {this->doStepImpl(forward);}
    void printSelf() {this->printSelfImpl();}
    // (Implementation: Virtual Functions)
protected:
    virtual void doStepImpl(bool forward) = 0;
    virtual void printSelfImpl() = 0;
};


template <unsigned char M,unsigned char N> 
class Reaction : public ReactionBase {
protected:
    std::array<Species*, M+N> _Species;
public:
    Reaction<M,N>(float f_rate, float b_rate, std::array<Species*, M+N> &species) : 
    ReactionBase(f_rate, b_rate), _Species(std::move(species)) {
        for(int i=0; i<M+N;++i){
            Species *s = this->_Species.at(i);
            bool left = (i<M) ? true : false;
            s->addReaction(this, left);
        }
    }
    ~Reaction<M,N>(){
        for(int i=0; i<M+N;++i)
            this->_Species.at(i)->removeReaction(this);
    }
public:
    void doStepImpl(bool forward){
        species_copy_incr_t delta; 
        for(int i=0; i<M+N;++i){
            delta = (forward && i<M) ? -1 : 1;
            Species *s = this->_Species.at(i);
            //std::cout << "i=" << i << ", delta=" << delta << "\n";
            s->incrementN(delta);
        }
    }
    void printSelfImpl() {
        std::cout << "Reaction, ID=" << _id << "\n";
        for(auto s: this->_Species){
            s->printSelf();
        }
    }
};


//struct ReactionNode {
//    float tau;
//    reaction_num_t rid;
//    bool forward;
//}; 

#endif
