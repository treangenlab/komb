//
// Created by Advait Balaji on 7/1/20.
//

#include "gfa.h"

namespace gfa
{
    Gfa::Gfa() {};
    Gfa::Gfa(uint32_t u1, uint32_t u2, std::string orient1, std::string orient2)
    {
        _u1 = u1;
        _u2 = u2;
        _orient1  = orient1;
        _orient2 = orient2;
    }

    std::string Gfa::getEdge(){
        return _u1+"$"+this->_u2;
    }

    uint32_t Gfa::getFirstUnitig() const{
        return _u1;
    }

    uint32_t Gfa::getSecondUnitig() const {
        return _u2;
    }

    std::string Gfa::getFirstUnitigOrientation() const{
        return _orient1;
    }

    std::string Gfa::getSecondUnitigOrientation() const{
        return _orient2;
    }
}