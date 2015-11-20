
//------------------------------------------------------------------
//  **M3SYM** - Simulation Package for the Mechanochemical
//              Dynamics of Active Networks, 3rd Generation
//
//  Copyright (2014) Papoian Lab, University of Maryland
//
//                 ALL RIGHTS RESERVED
//
//  See the Papoian lab page for more information:
//  http://papoian.chem.umd.edu/
//------------------------------------------------------------------

#ifndef M3SYM_Histogram_h
#define M3SYM_Histogram_h

#include <cmath>
#include <fstream>
#include <vector>

#include "common.h"

/// A class to hold frequency of occurences in a general set of data.

/*!
 *  The Histogram class holds frequencies of events in
 *  certain ranges (bins). It contains functions to
 *  increment and store frequencies in a data structure.
 */
class Histogram {
   
private:
    int _numBins; ///< Number of bins in histogram
    
    double _histMin; ///< Minimum histogram value
    double _histMax; ///< Maximum histogram value
    
    double _range; ///< range of bin

    vector<int> _frequencies; ///< The histogram data
    
public:
    Histogram(int numBins, double histMin = 0, double histMax = 0)
        : _numBins(numBins), _histMin(histMin), _histMax(histMax),
          _frequencies(numBins, 0) {

        //calculate range
        _range = (_histMax + _histMin) / _numBins;
    }
    
    ///Adds a value to the histogram by updating the correct bin
    void addValue(double value) {
        
        int bin = (int) (value / _range);
        assert(bin >= 0 && bin < _numBins && "Histogram error - trying to add value outside of bin range.");
        
        _frequencies[bin] += 1;
    }
    
    ///Clears all values from histogram
    void clearValues() {
        _frequencies.assign(_numBins, 0);
    }
    
    ///Print the histogram
    void print(ofstream& outputFile) {
        
        int bin = 0;
        for(auto freq : _frequencies) {
            
            float binMin = bin * _range;
            float binMax = (bin + 1) * _range;
            
            outputFile << binMin << " " << binMax << " " << freq << " ";
            
            bin++;
        }
    }
};


#endif