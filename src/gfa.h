//
// Created by Advait Balaji on 7/1/20.
//

#ifndef KOMB_GFA_H
#define KOMB_GFA_H

#include <string>
namespace gfa
{
    class Gfa
    {
        public:
            Gfa();
            Gfa(uint32_t u1, uint32_t u2, std::string orient1, std::string orient2);
            std::string getEdge();
            uint32_t getFirstUnitig() const ;
            uint32_t getSecondUnitig() const;
            std::string getFirstUnitigOrientation() const;
            std::string getSecondUnitigOrientation() const;

        private:
            uint32_t _u1, _u2;
            std::string _orient1, _orient2;
    };
}

#endif //KOMB_GFA_H
