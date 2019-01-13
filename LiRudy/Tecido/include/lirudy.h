#ifndef PRD_H
#define PRD_H

#include <iostream>
#include <cmath>
#include <vector>
#include <omp.h>
#include "../include/cell.h"

using namespace std;

class LiRudy
{
    static constexpr int NCELLS = 50;               // Number of cells
    static constexpr int NSC = 2;                   // Number of stimulus cells 
    static constexpr int OUT_ID = 0;                // Plot cell identifier

private:
    int nthreads;                                   // Number of threads
    int nbeats;                                     // Number of beats
    int n;                                          // Number of timesteps
    double dt;                                      // Size of the timestep
    double tmax;                                    // Maximum time 
    vector<Cell> cells;                             // Vector of cells

public:
    LiRudy (int argc, char *argv[]);
    void allocMem ();
    void setInitCond ();
    void setStimCells ();
    void setTimeSettings ();
    void solve ();

    // DEBUG
    friend ostream& operator<< (ostream &ost, const LiRudy &lirudy)
    {
        ost << "-------- PARAMETERS ------------" << endl;
        ost << "Dt = " << lirudy.dt << endl;
        ost << "tmax = " << lirudy.tmax << endl;
        ost << "Number of timesteps = " << lirudy.n << endl;
        ost << "--------------------------------" << endl;
        for (int i = 0; i < NCELLS; i++)
            ost << "Cell " << i << " -- v = " << lirudy.cells[i].v << endl;
        return ost;
    }

};

#endif