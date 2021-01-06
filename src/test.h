//
// Created by Advait Balaji on 6/28/20.
//

#ifndef TEST_H
#define TEST_H

#include <iostream>
#include <vector>
#include <memory>
#include <cstdlib>
#include <string>
#include <cstdio>

// vector to hold bases for the genome
std::vector<char> bases{'A','C','G','T'};

// instantiates random genome creator
class RandomGenome
{
    public:
        // constructor
         RandomGenome(uint8_t seed)
        {
            _seed =  seed;
        }
        // destructor
        ~RandomGenome() = default;

        // set seed
         uint8_t getSeed()
        {
            return _seed;
        }

        uint8_t getLength()
        {
             return _length;
        }

        void setLength(uint64_t length)
        {
             _length = length;
        }

        // generate the genome
        std::string generate()
        {
            srand(_seed);
            std::string genome = "";
            for( int i = 0 ; i < _length; i++)
            {
                uint8_t  randn = rand()%4; // randn wil be [0,4)
                genome+=bases[randn];
            }
            return genome;
        }

        void writeFasta(const std::string genome, const std::string fname)
        {
            FILE* f;
            f = fopen(fname.c_str(),"w+");
            fprintf(f,">randomfasta_%u_%lu\n", _seed,_length);
            fprintf(f,"%s",genome.c_str());
            fclose(f);
            fprintf(stdout, "Written random fasta file\n");

        }

    private:
        uint8_t _seed;
        uint64_t _length = 5000000;

};


#endif //TEST_H
