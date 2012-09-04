//
//  ReactionContainer.h
//  CytoSim
//
//  Created by Garegin Papoian on 8/31/12.
//  Copyright (c) 2012 University of Maryland. All rights reserved.
//

#ifndef CytoSim_ReactionContainer_h
#define CytoSim_ReactionContainer_h

#include <iostream>
#include "Reaction.h"

namespace  chem {
        
    class ReactionPtrContainerIFace {
    public:
        virtual Reaction* addReactionUnique(std::unique_ptr<Reaction> &&Reaction) = 0;
        virtual void removeReaction(Reaction* Reaction) = 0;
        virtual Reaction* findReaction (size_t index) = 0;
        virtual void printReaction() {}
    };
    
    
    
    
    class ReactionPtrContainerVector : public  ReactionPtrContainerIFace {
    protected:
        std::vector<std::unique_ptr<Reaction>> _reactions;
    public:
        virtual Reaction* addReactionUnique (std::unique_ptr<Reaction> &&Reaction) {
            _reactions.push_back(std::move(Reaction));
            return _reactions.back().get();
        }
        
        template<typename ...Args>
        Reaction* addReaction( Args&& ...args )
        {
            _reactions.push_back(std::unique_ptr<Reaction>( new Reaction( std::forward<Args>(args)...) ));
            //        _reactions.emplace_back(make_unique(Args...));
            return _reactions.back().get();
        }
        
        
        virtual void removeReaction (Reaction* R) {
            auto child_iter = std::find_if(_reactions.begin(),_reactions.end(),
                                           [R](const std::unique_ptr<Reaction> &element)
                                           {
                                               return element.get()==R ? true : false;
                                           });
            if(child_iter!=_reactions.end())
                _reactions.erase(child_iter);
        }
        
        virtual Reaction* findReaction (size_t index) {
            return _reactions[index].get();
        }
        
        virtual void printReactions() {
            for(auto &r : _reactions)
                std::cout << (*r.get());
        }
        
    };
    
} // end of chem


#endif
