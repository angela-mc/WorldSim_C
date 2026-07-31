//  Created by Angela on 5/22/20.
//  Copyright © 2020 Angela. All rights reserved.

//  main.cpp
//  cogsims_22mai


#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <time.h>
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <math.h>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <iterator>

#undef RAND_MAX
#define RAND_MAX 2147483647
#define PI 3.14159265

using namespace std;
//--------------------------------------------------------------------------------------------------------------
// *******************************************
//  MODEL PARAMETERS AND PHENOTYPE DEFINITIONS
// *******************************************

const int Runs = 1; // # of replicate simulation runs
const int NumGen = 100000;  // # number of generations per replicate

// mutation and sex parameters
const double mu = 1e-03; // mutation rate for gene and pleiotropy values
const double sdmu = 0.05; // mutation step sizes (as a PERCENTAGE of a feasible trait range)
const double RepThreshold = 0.3; // max Euclidean distance between P vectors allowed for reproduction (0.6 Assumes a difference of ca. 0.2 in every locus)

// DISPERSAL AND COGNITION
const double CogBuff = 0; // cognitive ability [0 = no buffering capability; 1 = complete cognitive buffer]

//const double Pdisp = 0.1; // probability of dispersing to other cell
//const int Dmax = 3; // dispersal ability (in max number of cells from natal cell)
const double lambda_dist = 0.5; // Poisson distribution mean for short dispersal events
const double lambda_dist_long = 2; // Poisson distribution mean for long dispersal events
const double threshold_short=2;
//const double lambda_dist_mean = 7; //normal mean for dipersal tail
//const double lambda_dist_sd = 0.75; //normal sd for dipersal tail

// dispersal parameters
//const double closestlandcell = 0;

//const int InitialCells = 5; // seeds initial population to start the model (only creates individuals in these cells
const int InitialCellID = 1776; //starting cell


// Environment
const int Ncell = 7345; // number of cells in the world
    // prior: (MUST ALLOW FORMATION OF X by 2X matrix with EVEN number of rows), now NO
//delete sd_E from here
//const double sd_E = 0.02; // sd of variation in environmental conditions (each cell has a specific environmental value and this parameter controls the amount of noise around it)

const int CarryingCapacity=100; // maximum number of individuals that are supported by each cell in the world

struct Cell {
    // Environmental value - angela: add sd_E
    double Ebar, ThisE, sd_E; // mean environmental value for the cell and place holder for realized value each generation
    
    // Carrying Capacity
    int MaxOccupants;
    
    // Neighboring cells
    int Neighbor_ids[7]; //last one = closest land cell for water cells
    
};

Cell World[Ncell];

struct Population {
    // number of individuals that are actually alive in population (the vectors below preallocate more individuals in case of dispersal)
    int PopSize;
    
    // Vectors containing the 10 loci genotypes for all individuals
    double Pmom1[CarryingCapacity*5], Pdad1[CarryingCapacity*5];
    double Pmom2[CarryingCapacity*5], Pdad2[CarryingCapacity*5];
    double Pmom3[CarryingCapacity*5], Pdad3[CarryingCapacity*5];
    double Pmom4[CarryingCapacity*5], Pdad4[CarryingCapacity*5];
    double Pmom5[CarryingCapacity*5], Pdad5[CarryingCapacity*5];
    double Pmom6[CarryingCapacity*5], Pdad6[CarryingCapacity*5];
    double Pmom7[CarryingCapacity*5], Pdad7[CarryingCapacity*5];
    double Pmom8[CarryingCapacity*5], Pdad8[CarryingCapacity*5];
    double Pmom9[CarryingCapacity*5], Pdad9[CarryingCapacity*5];
    double Pmom10[CarryingCapacity*5], Pdad10[CarryingCapacity*5];
    
    // relative fitness against others in own populatiom
    double W[CarryingCapacity*5], Wrel[CarryingCapacity*5];
};

Population LocalPopulations[Ncell], Offsp_LocalPopulations[Ncell];

// Costs
const double k = 5; // modulates exponential decay of fitness with mismatch
const double c = 0; // Fitness cost of cognitive ability (should be a small value)

int Generation, Seed, Replicate;

//--------------------------------------------------------------------------------------------------------------
// ******************************
// pseudo-random number generator
// ******************************
int rand_a_state = 0;
void srand_a(int value) {
    rand_a_state = value;
}
int rand_a() {
    rand_a_state = ((rand_a_state * 1103515245) + 12345) & 0x7fffffff;
    return rand_a_state;
}

// random number generator type double (range between 0 and 1)
double ru() { return (double)(rand_a() % 10001) / 10000; }

// random integer generator (between 0 and N)
int rn(int N) { return rand_a() % (N + 1); }

double Poisson(double mean)
{
    double ran = ru();
    double P = exp(-mean);
    double sum = P;
    int i;
    
    if (ran<sum) { return 0; }
    for (i = 1; i<mean * 5; i++)
    {
        P = P*mean / i;
        sum = sum + P;
        if (ran<sum) { return i; }
    }
    return mean * 5; // around a clutch size = 5
}

double Normal(double mean, double stddev)
{//Box muller method
    static double n2 = 0.0;
    static int n2_cached = 0;
    if (!n2_cached)
    {
        double x, y, r;
        do
        {
            x = 2.0*rand_a() / RAND_MAX - 1;
            y = 2.0*rand_a() / RAND_MAX - 1;
            
            r = x*x + y*y;
        } while (r == 0.0 || r > 1.0);
        {
            double d = sqrt(-2.0*log(r) / r);
            double n1 = x*d;
            n2 = y*d;
            double result = n1*stddev + mean;
            n2_cached = 1;
            return result;
        }
    }
    else
    {
        n2_cached = 0;
        return n2*stddev + mean;
    }
}

// find element in vector - can replace it by - find C fxn
int fxnfind(int vector[], int elem, int poz)
{ int found=1;
  for(int k=0;k<=poz;k++)
    {if(elem==vector[k]) found=0;}
    
return(found);}

//--------------------------------------------------------------------------------------------------------------
// ****************
// Useful functions
// ****************

// initialize ouput text files
void WriteHeaders_Pfile(ofstream &ParametersFile) {
    
    ParametersFile << "Runs: " << Runs << endl
    << "NumGen: " << NumGen << endl
    << "mu: " << mu << endl
    << "sdmu: " << sdmu << endl << endl
    
    << "** Model parameters ***" << endl
    << "N cells in world = " << Ncell << endl << endl
    << "Initial seed (number of cells carrying individuals at t = 0): one a priori chosen land & world cell + all nb" << endl
    << "Initial cellid: " << InitialCellID <<endl
    << "Carrying capacity: " << CarryingCapacity << endl
    << "Lambda distance short: "<<lambda_dist<<endl
    << "Lambda distance long: "<<lambda_dist_long<<endl
    << "Threshold short: "<<threshold_short<<endl
    
    << "Exp decay of W (k): " << k << endl
    << "Cost of cognition (c): " << c << endl << endl
    
    << "Cognitive buffer: " << CogBuff << endl
    
    << "*************"
    << "Neighborhood structure"<< endl<<endl
    
    << "CelID"<< '\t' << "UL" << '\t' << "UR" << '\t' << "L" << '\t' << "R" << '\t' << "LL" << '\t' << "LR" << '\t'<< "ClosestLandCell" << endl;
    
    for (int i=0; i<Ncell; i++) {
        ParametersFile << i << '\t' << setprecision(6) << World[i].Neighbor_ids[0]<< '\t' << setprecision(6) << World[i].Neighbor_ids[1]<< '\t'
                                    << setprecision(6) << World[i].Neighbor_ids[2]<< '\t' << setprecision(6) << World[i].Neighbor_ids[3]<< '\t'
                                    << setprecision(6) << World[i].Neighbor_ids[4]<< '\t' << setprecision(6) << World[i].Neighbor_ids[5]<<
                                        '\t'
                                    <<setprecision(6)<<World[i].Neighbor_ids[6]<<
                                    endl;
    }
    
    ParametersFile << "*************" << endl
    << "Run     Seed"
    << endl;
    
 
    
}

void WriteHeaders_Dfile(ofstream &DistFile) {
    
       
    DistFile << "Run" << '\t' << "Generation" << '\t' << "Individual" << '\t'
    << "P1.mom" << '\t' << "P2.mom" << '\t' << "P3.mom" << '\t'
    << "P4.mom" << '\t' << "P5.mom" << '\t' << "P6.mom" << '\t'
    << "P7.mom" << '\t' << "P8.mom" << '\t' << "P9.mom" << '\t'
    << "P10.mom" << '\t'
    
    << "P1.dad" << '\t' << "P2.dad" << '\t' << "P3.dad" << '\t'
    << "P4.dad" << '\t' << "P5.dad" << '\t' << "P6.dad" << '\t'
    << "P7.dad" << '\t' << "P8.dad" << '\t' << "P9.dad" << '\t'
    << "P10.dad" << '\t'
    
    << "Cell ID" << '\t' << endl;
    
}

// ouput entire population to a textfile
void WriteDist(ofstream &DistFile) {
    int counter = 1;
    for (int i = 0; i < Ncell; i++) {
        for (int j = 0; j < LocalPopulations[i].PopSize; j++) {
            DistFile << Replicate
            << '\t' << Generation
            << '\t' << counter
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom1[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom2[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom3[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom4[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom5[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom6[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom7[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom8[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom9[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pmom10[j]
            
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad1[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad2[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad3[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad4[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad5[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad6[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad7[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad8[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad9[j]
            << '\t' << setprecision(6) << LocalPopulations[i].Pdad10[j]
            
            << '\t' << setprecision(6) << i
            << endl;
            
            counter += 1;
        }
    }
}

double LimitCheck(double value, double up, double down)
{
    // checks that value is within bounds (for bounded parameters)
    if (value>up) value = up;
    if (value<down) value = down;
    return value;
}


// Fitness function
double GetW(int cellID, int indID, double E) {
    
    // establish individual phenotype (sum of averages of the two copies at each locus)
    double P = 0;
    P += (LocalPopulations[cellID].Pdad1[indID] + LocalPopulations[cellID].Pmom1[indID])/2;
    P += (LocalPopulations[cellID].Pdad2[indID] + LocalPopulations[cellID].Pmom2[indID])/2;
    P += (LocalPopulations[cellID].Pdad3[indID] + LocalPopulations[cellID].Pmom3[indID])/2;
    P += (LocalPopulations[cellID].Pdad4[indID] + LocalPopulations[cellID].Pmom4[indID])/2;
    P += (LocalPopulations[cellID].Pdad5[indID] + LocalPopulations[cellID].Pmom5[indID])/2;
    P += (LocalPopulations[cellID].Pdad6[indID] + LocalPopulations[cellID].Pmom6[indID])/2;
    P += (LocalPopulations[cellID].Pdad7[indID] + LocalPopulations[cellID].Pmom7[indID])/2;
    P += (LocalPopulations[cellID].Pdad8[indID] + LocalPopulations[cellID].Pmom8[indID])/2;
    P += (LocalPopulations[cellID].Pdad9[indID] + LocalPopulations[cellID].Pmom9[indID])/2;
    P += (LocalPopulations[cellID].Pdad10[indID] + LocalPopulations[cellID].Pmom10[indID])/2;
     
    // calculate fitness based on mismatch
    return (exp(-k*abs(P - CogBuff*(P-E) - E)) - c*CogBuff);
    //-CogBuff ie (P-E) will be smaller than it s actual value
}

// evaluate mating compatibility
double matingCompatibility(int cellID, int male, int female) {
    double compatibility = 0;
    compatibility += pow((LocalPopulations[cellID].Pdad1[male] + LocalPopulations[cellID].Pmom1[male])/2 -
                         (LocalPopulations[cellID].Pdad1[female] + LocalPopulations[cellID].Pmom1[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad2[male] + LocalPopulations[cellID].Pmom2[male])/2 -
                         (LocalPopulations[cellID].Pdad2[female] + LocalPopulations[cellID].Pmom2[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad3[male] + LocalPopulations[cellID].Pmom3[male])/2 -
                         (LocalPopulations[cellID].Pdad3[female] + LocalPopulations[cellID].Pmom3[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad4[male] + LocalPopulations[cellID].Pmom4[male])/2 -
                         (LocalPopulations[cellID].Pdad4[female] + LocalPopulations[cellID].Pmom4[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad5[male] + LocalPopulations[cellID].Pmom5[male])/2 -
                         (LocalPopulations[cellID].Pdad5[female] + LocalPopulations[cellID].Pmom5[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad6[male] + LocalPopulations[cellID].Pmom6[male])/2 -
                         (LocalPopulations[cellID].Pdad6[female] + LocalPopulations[cellID].Pmom6[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad7[male] + LocalPopulations[cellID].Pmom7[male])/2 -
                         (LocalPopulations[cellID].Pdad7[female] + LocalPopulations[cellID].Pmom7[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad8[male] + LocalPopulations[cellID].Pmom8[male])/2 -
                         (LocalPopulations[cellID].Pdad8[female] + LocalPopulations[cellID].Pmom8[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad9[male] + LocalPopulations[cellID].Pmom9[male])/2 -
                         (LocalPopulations[cellID].Pdad9[female] + LocalPopulations[cellID].Pmom9[female])/2, 2);
    
    compatibility += pow((LocalPopulations[cellID].Pdad10[male] + LocalPopulations[cellID].Pmom10[male])/2 -
                         (LocalPopulations[cellID].Pdad10[female] + LocalPopulations[cellID].Pmom10[female])/2, 2);
    
    compatibility = sqrt(compatibility); // Euclidean distance
    
    // calculate fitness based on mismatch
    return (compatibility);
}

// mutates a given allele within certain bounds
double ThisMutation(double currValue, double sd, double LowLimit, double HighLimit) {
    if (ru()<mu)
    {
        currValue += Normal(0, sd);
        currValue = LimitCheck(currValue, HighLimit, LowLimit);
    }
    return currValue;
}

// initializes population at the begining of a new replicate
void Init(ofstream &ParametersFile) {
    // WE NEED A BETTER SEED!!!
    int Seed = static_cast<unsigned int>(time(NULL));
    srand_a(Seed);
    
    // create a world and initialize environmental parameters
    for (int i = 0; i < Ncell; i++) {
             
        // don t limit check
        World[i].ThisE = Normal(World[i].Ebar, World[i].sd_E);
              
        // reset local populations of offspring
        Offsp_LocalPopulations[i].PopSize = 0;
    }
    
    // store this info in case we want to rerun a replicate
    ParametersFile << Replicate << '\t' << Seed << endl;
    
    // initialize population (seed InitialCellID+ its immediate neighbors IF these are not out of the world or water)
    int ncells[7]={-2,-2,-2,-2,-2,-2,-2}; //maximum: initial cell + 6 nbs; initialize with -2 because no cell is -2, so if the ncells vector is not populated with true cells (if nb of InitialCellID are water/out of the world), it will have -2, and when it will check i in Ncell - it won't find it
    
    ncells[0]=InitialCellID;
    int c_ncells=0; //pozition, should go from 0-6 (max)
    
    if(World[InitialCellID].Neighbor_ids[0]!=-1)
        if(World[World[InitialCellID].Neighbor_ids[0]].Ebar!=999 & World[World[InitialCellID].Neighbor_ids[0]].Ebar!=-999)
            {c_ncells=c_ncells+1; ncells[c_ncells]=World[InitialCellID].Neighbor_ids[0];}
    if(World[InitialCellID].Neighbor_ids[1]!=-1)
        if(World[World[InitialCellID].Neighbor_ids[1]].Ebar!=999 & World[World[InitialCellID].Neighbor_ids[1]].Ebar!=-999)
            {c_ncells=c_ncells+1; ncells[c_ncells]=World[InitialCellID].Neighbor_ids[1];}
    if(World[InitialCellID].Neighbor_ids[2]!=-1)
        if(World[World[InitialCellID].Neighbor_ids[2]].Ebar!=999 & World[World[InitialCellID].Neighbor_ids[2]].Ebar!=-999)
            {c_ncells=c_ncells+1; ncells[c_ncells]=World[InitialCellID].Neighbor_ids[2];}
    if(World[InitialCellID].Neighbor_ids[3]!=-1)
        if(World[World[InitialCellID].Neighbor_ids[3]].Ebar!=999 & World[World[InitialCellID].Neighbor_ids[3]].Ebar!=-999)
            {c_ncells=c_ncells+1; ncells[c_ncells]=World[InitialCellID].Neighbor_ids[3];}
    if(World[InitialCellID].Neighbor_ids[4]!=-1)
        if(World[World[InitialCellID].Neighbor_ids[4]].Ebar!=999 & World[World[InitialCellID].Neighbor_ids[4]].Ebar!=-999)
            {c_ncells=c_ncells+1; ncells[c_ncells]=World[InitialCellID].Neighbor_ids[4];}
    if(World[InitialCellID].Neighbor_ids[5]!=-1)
        if(World[World[InitialCellID].Neighbor_ids[5]].Ebar!=999 & World[World[InitialCellID].Neighbor_ids[5]].Ebar!=-999)
            {c_ncells=c_ncells+1; ncells[c_ncells]=World[InitialCellID].Neighbor_ids[5];}

    cout<<endl<<"All initial cells are: ";
    for(int c_ncellsi=0; c_ncellsi<=c_ncells; c_ncellsi++)
        cout<<ncells[c_ncellsi]<<" ";
    /*
    int ncells[7] = {InitialCellID, World[InitialCellID].Neighbor_ids[0],
                     World[InitialCellID].Neighbor_ids[1],
                     World[InitialCellID].Neighbor_ids[2],
                     World[InitialCellID].Neighbor_ids[3],
                     World[InitialCellID].Neighbor_ids[4],
                     World[InitialCellID].Neighbor_ids[5]}; */
    for (int i = 0; i < Ncell; i++) {
        if (std::find(std::begin(ncells), std::end(ncells), i) != std::end(ncells)) {
            LocalPopulations[i].PopSize = CarryingCapacity;
            for (int j=0; j < CarryingCapacity; j++) {
                // assume these guys are already reasonably adapted to local conditions and are a single species
                LocalPopulations[i].Pmom1[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1,-1); //limitcheck has high limit - low limit
                LocalPopulations[i].Pdad1[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom2[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad2[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom3[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad3[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom4[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad4[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom5[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad5[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom6[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad6[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom7[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad7[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom8[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad8[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom9[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad9[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pmom10[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);
                LocalPopulations[i].Pdad10[j] = LimitCheck(Normal(World[i].Ebar, 0.05)/10, 1, -1);

            }
        }
        else {
            LocalPopulations[i].PopSize = 0;
        }
    }
}


// This is where individuals get to use their genes!
void OneGeneration() {
    
    
    // reproduce (only allow reproduction with occupants of the same cell)
    for (int i=0; i<Ncell; i++) {
        if (LocalPopulations[i].PopSize > 0) {
            /*
            // Compute mean fitness of ALL (resident and immigrant) individuals
            double meanW = 0;
            for (int j=0; j<LocalPopulations[i].PopSize; j++) {
                meanW += ( GetW(i, j, World[i].ThisE) / LocalPopulations[i].PopSize );
            }
            */
            double meanW=0.5; // scale N_off to 0.5 value
            
            // compute relative fitness & N_off of all
            int N_off[CarryingCapacity*100]; // if each individ has 100 offsprings
            int ind1[CarryingCapacity*100];
            int ind2[CarryingCapacity*100];
            int cind2=-1; //position of elements
            int cind1=-1; //pozition of elements
            
            for (int j=0; j<LocalPopulations[i].PopSize; j++) {
                // this can be commented out once we will have meanW > 0 bc of function
                if (meanW>0) {
                    // Compute fitness relative to other competitors in this cell
                    LocalPopulations[i].Wrel[j] = GetW(i, j, World[i].ThisE) / meanW;
                    N_off[j] = round(Poisson(LocalPopulations[i].Wrel[j]));
                }
                else {
                    LocalPopulations[i].Wrel[j] = 0;
                    N_off[j] = 0;
                }
                
               // build representation of genes in offspring pool
               if (N_off[j]>0) //if 0, this individual does not reproduce at all and its genotype is      skipped from pool
                  for(int j_noff=0; j_noff<N_off[j]; j_noff++)
                      {cind1++; cind2++;
                       ind1[cind1]=j; ind2[cind2]=j_noff; //pozition of element
                       }
                }
            
            // reproduce proportionally to local success (approximated by N_off[])
            if(cind1>-1) // if any individuals have at least one N_off
              {
               int currindv;//int currpoz;
               int partner; int partnerpoz;
                  
               if(cind1==0) // if cind1==0 ind1 has 1 element, and it must self-reproduction
                  {currindv=ind1[0];
                   partner=ind1[0]; //1 element - currinds = partner  = that element
                      
                    // now they are allowed to have 1 offspring
                    // disperse?
                    int newcell = i;
                    int Dmax = round(Poisson(lambda_dist));
                    if(Dmax>=threshold_short)
                      {Dmax = round(Poisson(lambda_dist_long));
                       if(Dmax<threshold_short) Dmax=threshold_short; // if the long dist is <3 - set it at 3 ie min of long-dist
                        }
                                                                                
                    int poolnb[8000];
                    int cpnb=-1;
                    int c1=0; int c2=0;
                    int oknb[8000];
                    int c_oknb=-1;
                                                                               
                    // if Dmax==0 --> newcell remains i
                    // must add if Dmax>0 otherwise D==0 can produce an error on this compiler...
                    if(Dmax>0)
                      {for (int D=0; D<Dmax; D++)
                           {// first step - nb of i - do not worry about uniqueness (all nbs are inevitably unique)
                            // the nb of i can either be -1; or a cell of out-of-the-world in which case do not add
                            if(D==0)
                                    {                                 if(World[i].Neighbor_ids[0]!=-1)
                                                                        if(World[World[i].Neighbor_ids[0]].Ebar!=999)
                                                                       {cpnb=cpnb+1;
                                                                        poolnb[cpnb] = World[i].Neighbor_ids[0];
                                                                       }
                                                                       if(World[i].Neighbor_ids[1]!=-1)
                                                                         if(World[World[i].Neighbor_ids[1]].Ebar!=999)
                                                                            {cpnb=cpnb+1;
                                                                             poolnb[cpnb] = World[i].Neighbor_ids[1];
                                                                            }
                                                                       if(World[i].Neighbor_ids[2]!=-1)
                                                                          if(World[World[i].Neighbor_ids[2]].Ebar!=999)
                                                                             {cpnb=cpnb+1;
                                                                              poolnb[cpnb] = World[i].Neighbor_ids[2];
                                                                             }
                                                                       if(World[i].Neighbor_ids[3]!=-1)
                                                                          if(World[World[i].Neighbor_ids[3]].Ebar!=999)
                                                                            {cpnb=cpnb+1;
                                                                             poolnb[cpnb] = World[i].Neighbor_ids[3];
                                                                            }
                                                                       if(World[i].Neighbor_ids[4]!=-1)
                                                                          if(World[World[i].Neighbor_ids[4]].Ebar!=999)
                                                                             {cpnb=cpnb+1;
                                                                              poolnb[cpnb] = World[i].Neighbor_ids[4];
                                                                             }
                                                                       if(World[i].Neighbor_ids[5]!=-1)
                                                                          if(World[World[i].Neighbor_ids[5]].Ebar!=999)
                                                                            {cpnb=cpnb+1;
                                                                             poolnb[cpnb] = World[i].Neighbor_ids[5];
                                                                            }
                                c1=0; c2=cpnb; //it will at least add one - don t have cells with no in-this-world nb
                                } // from if (D==0)
                                                                                      
                        // if D>0 ie next step; add nbs of previous stage
                        // add only unique nbs - avoid building ultra large vectors
                        //if(Dmax>1)
                                if(D>0)
                                  {for(int jj=c1;jj<=c2;jj++) //poolnb[[jj]] can be last row, its nb will just be -1
                                     {if (World[poolnb[jj]].Neighbor_ids[0] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[0]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[0], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[0];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[1] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[1]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[1], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[1];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[2] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[2]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[2], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[2];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[3] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[3]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[3], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[3];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[4] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[4]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[4], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[4];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[5] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[5]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[5], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[5];
                                        }

                                       } // from (for jj)
                                                                                         
                                   if(c2!=cpnb) //if cpnb changed ie at least one element was added
                                      {c1=c2+1;
                                       c2=cpnb;}
                                  //if(c2==cpnb) //no new elements were added (is this even possible?)
                                  // cpnb remains same, so do c1 and c2
                                    
                                  }// from if (D>0)
                             }//from for(D)
                                                                                    
                        // select nonwater and non-i cells (all should be unique and not -1 or out of the world)
                        for(int hh=0;hh<=cpnb;hh++)
                           {if(World[poolnb[hh]].Ebar!=-999 & poolnb[hh]!=i)
                              {c_oknb=c_oknb+1;
                               oknb[c_oknb]=poolnb[hh];
                              }
                           } //from for(hh)
                                                                                    
                       // select randomly a cell from oknb if any ok nb (ie not all water)
                       if(c_oknb>-1) newcell = oknb[rn(c_oknb)];
                         else newcell=i; //if no OK nb, only water, then stay put
                      } //from if Dmax>0
                                                        
                    //disperse to new cell: this will be either i (if Dmax==0) or a new cell if (Dmax>0)
                    int j=currindv; //to match old code
                    if (Offsp_LocalPopulations[newcell].PopSize >= CarryingCapacity * 5)
                            {
                                continue;
                            }

                                                                if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom1[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad1[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom1[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad1[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom2[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad2[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom2[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad2[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom3[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad3[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom3[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad3[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom4[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad4[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom4[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad4[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom5[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad5[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom5[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad5[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom6[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad6[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom6[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad6[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom7[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad7[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom7[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad7[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom8[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad8[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom8[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad8[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom9[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad9[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom9[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad9[partner]; }
                                                                                          
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom10[j]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pmom10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad10[j]; }
                                                                  if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom10[partner]; }
                                                                  else { Offsp_LocalPopulations[newcell].Pdad10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad10[partner]; }
                                                                  
                                                                  Offsp_LocalPopulations[newcell].PopSize += 1;

                   } //if cind1==0
               else
                {for(int cc=0; cc<=cind1;cc++)
                  { currindv= ind1[cc]; //cout<<"cc= "<<cc<<" currindv= "<<currindv<<" ind2[cc]= "<<ind2[cc]<<endl;
                    if(ind2[cc]>-1) //has not yet been chosen as partner
                    {ind2[cc]=-1; //this position is not longer available
                     //cout<<"ind2[cc]now= "<<ind2[cc]<<endl;
                  
                    // is there another ind2 that is not -1 at this point?
                    int sumind2=0;
                    for(int l=0;l<=cind2;l++) sumind2=sumind2+ind2[l];
                    //cout<<endl<<"sumind2= "<<sumind2<<endl;
                    if(sumind2==((cind2+1)*(-1))) // all -1 at this point - no partner - must self
                       {partner=currindv;} //cout<<"partner=currdinv= "<<partner<<endl;
                    else // find partner
                       {bool notFound = true;
                        bool notFound2 = true;
                        int RepCount = 1;
                        do{ //find potential partner
                            notFound2 = true;
                          do {partnerpoz = rn(cind2);
                              if(ind2[partnerpoz]>-1) notFound2=false; //if not -1 => good position
                              } while(notFound2);
                          partner= ind1[partnerpoz];
                          //cout<<"candidate partner= "<<partner<<" partnerpoz"<<partnerpoz<<endl;
                          //if (0==0) //check partner do smthing so that sometimes it doesn t check...rn(1)==0??
                          if (matingCompatibility(i, currindv, partner) < RepThreshold)
                             {notFound = false;}
                            //can add a term to say AND it is not itself as self should be separate
                          // include this clause to allow some reproduction in rare or invading genotypes
                          if ((RepCount == 20) & notFound)
                             {partner = currindv;
                              partnerpoz=cc;
                              notFound = false;
                              }
                          RepCount+=1;
                           } while (notFound);
                          //cout<<"oked partner= "<<partner<<" partnerpoz= "<<partnerpoz<<endl;
                        
                          // erase partner from ind1 only if it has not self reproduced
                          if(RepCount!=21) //
                             {ind2[partnerpoz]=-1;
                             }
                         } // from else
                     //cout<<"ind2[partnerpoz]now= "<<ind2[partnerpoz]<<endl;
                     } //if(ind2[cc]>-1), if it is - it has been chosen - move to next cc
                      
                      
                     // now they are allowed to have 1 offspring (either via self OR reproduction)
                     // disperse?
                     int newcell = i;
                     int Dmax = round(Poisson(lambda_dist));
                     if(Dmax>=3)
                       {Dmax = round(Poisson(lambda_dist_long));
                        if(Dmax<3) Dmax=3; // if the long dist is <3 - set it at 3 ie min of long-dist
                        }
                                                           
                    int poolnb[8000];
                    int cpnb=-1;
                    int c1=0; int c2=0;
                    int oknb[8000];
                    int c_oknb=-1;
                                                          
                    // if Dmax==0 --> newcell remains i
                    // must add if Dmax>0 otherwise D==0 can produce an error on this compiler...
                    if(Dmax>0)
                      {for (int D=0; D<Dmax; D++)
                        {// first step - nb of i - do not worry about uniqueness (all nbs are inevitably unique)
                         // the nb of i can either be -1; or a cell of out-of-the-world in which case do not add
                         if(D==0)
                           {if(World[i].Neighbor_ids[0]!=-1)
                             if(World[World[i].Neighbor_ids[0]].Ebar!=999)
                             {cpnb=cpnb+1;
                              poolnb[cpnb] = World[i].Neighbor_ids[0];
                             }
                                                  if(World[i].Neighbor_ids[1]!=-1)
                                                    if(World[World[i].Neighbor_ids[1]].Ebar!=999)
                                                       {cpnb=cpnb+1;
                                                        poolnb[cpnb] = World[i].Neighbor_ids[1];
                                                       }
                                                  if(World[i].Neighbor_ids[2]!=-1)
                                                     if(World[World[i].Neighbor_ids[2]].Ebar!=999)
                                                        {cpnb=cpnb+1;
                                                         poolnb[cpnb] = World[i].Neighbor_ids[2];
                                                        }
                                                  if(World[i].Neighbor_ids[3]!=-1)
                                                     if(World[World[i].Neighbor_ids[3]].Ebar!=999)
                                                       {cpnb=cpnb+1;
                                                        poolnb[cpnb] = World[i].Neighbor_ids[3];
                                                       }
                                                  if(World[i].Neighbor_ids[4]!=-1)
                                                     if(World[World[i].Neighbor_ids[4]].Ebar!=999)
                                                        {cpnb=cpnb+1;
                                                         poolnb[cpnb] = World[i].Neighbor_ids[4];
                                                        }
                                                  if(World[i].Neighbor_ids[5]!=-1)
                                                     if(World[World[i].Neighbor_ids[5]].Ebar!=999)
                                                       {cpnb=cpnb+1;
                                                        poolnb[cpnb] = World[i].Neighbor_ids[5];
                                                       }
                           c1=0; c2=cpnb; //it will at least add one - don t have cells with no in-this-world nb
                           } // from if (D==0)
                                                                 
                           // if D>0 ie next step; add nbs of previous stage
                           // add only unique nbs - avoid building ultra large vectors
                           //if(Dmax>1)
                           if(D>0)
                             {for(int jj=c1;jj<=c2;jj++) //poolnb[[jj]] can be last row, its nb will just be -1
                                {if (World[poolnb[jj]].Neighbor_ids[0] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[0]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[0], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[0];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[1] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[1]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[1], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[1];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[2] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[2]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[2], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[2];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[3] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[3]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[3], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[3];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[4] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[4]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[4], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[4];
                                        }

                                        if (World[poolnb[jj]].Neighbor_ids[5] != -1 &&
                                            World[World[poolnb[jj]].Neighbor_ids[5]].Ebar != 999 &&
                                            fxnfind(poolnb, World[poolnb[jj]].Neighbor_ids[5], cpnb) == 1)
                                        {
                                            cpnb = cpnb + 1;
                                            poolnb[cpnb] = World[poolnb[jj]].Neighbor_ids[5];
                                        }

                                       }
                                                                    
                            if(c2!=cpnb) //if cpnb changed ie at least one element was added
                                {c1=c2+1;
                                 c2=cpnb;}
                            //if(c2==cpnb) //no new elements were added (is this even possible?)
                            // cpnb remains same, so do c1 and c2
                                                      
                            }// from if (D>0)
                         }//from for(D)
                                                               
                        // select nonwater and non-i cells (all should be unique and not -1 or out of the world)
                        for(int hh=0;hh<=cpnb;hh++)
                          {if(World[poolnb[hh]].Ebar!=-999 & poolnb[hh]!=i)
                             {c_oknb=c_oknb+1;
                              oknb[c_oknb]=poolnb[hh];
                              }
                          } //from for(hh)
                                                               
                        // select randomly a cell from oknb if any ok nb (ie not all water)
                        if(c_oknb>-1) newcell = oknb[rn(c_oknb)];
                          else newcell=i; //if no OK nb, only water, then stay put
                      } //from if Dmax>0
                                   
                      //disperse to new cell: this will be either i (if Dmax==0) or a new cell if (Dmax>0)
                      int j=currindv; //to match old code
                                           if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom1[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad1[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom1[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad1[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad1[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom2[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad2[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom2[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad2[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad2[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom3[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad3[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom3[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad3[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad3[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom4[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad4[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom4[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad4[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad4[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom5[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad5[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom5[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad5[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad5[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom6[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad6[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom6[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad6[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad6[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom7[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad7[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom7[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad7[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad7[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom8[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad8[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom8[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad8[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad8[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom9[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad9[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom9[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad9[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad9[partner]; }
                                                                     
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pmom10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom10[j]; }
                                             else { Offsp_LocalPopulations[newcell].Pmom10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad10[j]; }
                                             if (rn(1)==0) { Offsp_LocalPopulations[newcell].Pdad10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pmom10[partner]; }
                                             else { Offsp_LocalPopulations[newcell].Pdad10[Offsp_LocalPopulations[newcell].PopSize] = LocalPopulations[i].Pdad10[partner]; }
                                             
                                             Offsp_LocalPopulations[newcell].PopSize += 1;

                                       
                   } //for (cc) loop
                } // else from if cind1==0
              } //if cind1==-1
                        
        } // from if pop size in cell i >0
    } // from for(i) each cell
    

    // make sure the offspring produce do not exceed carrying capacity
    for (int i=0; i<Ncell; i++) {
        if (Offsp_LocalPopulations[i].PopSize > CarryingCapacity) {
            do {
                // remove extra individuals chosen at ramdom
                int RM = rn(Offsp_LocalPopulations[i].PopSize - 1);
                if (RM < Offsp_LocalPopulations[i].PopSize) {
                    for (int j = RM; j<Offsp_LocalPopulations[i].PopSize-1; j++) {
                        Offsp_LocalPopulations[i].Pmom1[j] = Offsp_LocalPopulations[i].Pmom1[j+1];
                        Offsp_LocalPopulations[i].Pmom2[j] = Offsp_LocalPopulations[i].Pmom2[j+1];
                        Offsp_LocalPopulations[i].Pmom3[j] = Offsp_LocalPopulations[i].Pmom3[j+1];
                        Offsp_LocalPopulations[i].Pmom4[j] = Offsp_LocalPopulations[i].Pmom4[j+1];
                        Offsp_LocalPopulations[i].Pmom5[j] = Offsp_LocalPopulations[i].Pmom5[j+1];
                        Offsp_LocalPopulations[i].Pmom6[j] = Offsp_LocalPopulations[i].Pmom6[j+1];
                        Offsp_LocalPopulations[i].Pmom7[j] = Offsp_LocalPopulations[i].Pmom7[j+1];
                        Offsp_LocalPopulations[i].Pmom8[j] = Offsp_LocalPopulations[i].Pmom8[j+1];
                        Offsp_LocalPopulations[i].Pmom9[j] = Offsp_LocalPopulations[i].Pmom9[j+1];
                        Offsp_LocalPopulations[i].Pmom10[j] = Offsp_LocalPopulations[i].Pmom10[j+1];
                    }
                }
                Offsp_LocalPopulations[i].PopSize -= 1;
            } while (Offsp_LocalPopulations[i].PopSize > CarryingCapacity);
        }
    }

    // now we can allow surviving Offspring to become reproductive
    for (int i=0; i<Ncell; i++) {
        LocalPopulations[i] = Offsp_LocalPopulations[i];
        Offsp_LocalPopulations[i].PopSize = 0;
        
        // and this is where we enable mutation
        for (int j=0; j<LocalPopulations[i].PopSize; j++) {
            LocalPopulations[i].Pmom1[j] = ThisMutation(LocalPopulations[i].Pmom1[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad1[j] = ThisMutation(LocalPopulations[i].Pdad1[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom2[j] = ThisMutation(LocalPopulations[i].Pmom2[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad2[j] = ThisMutation(LocalPopulations[i].Pdad2[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom3[j] = ThisMutation(LocalPopulations[i].Pmom3[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad3[j] = ThisMutation(LocalPopulations[i].Pdad3[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom4[j] = ThisMutation(LocalPopulations[i].Pmom4[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad4[j] = ThisMutation(LocalPopulations[i].Pdad4[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom5[j] = ThisMutation(LocalPopulations[i].Pmom5[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad5[j] = ThisMutation(LocalPopulations[i].Pdad5[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom6[j] = ThisMutation(LocalPopulations[i].Pmom6[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad6[j] = ThisMutation(LocalPopulations[i].Pdad6[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom7[j] = ThisMutation(LocalPopulations[i].Pmom7[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad7[j] = ThisMutation(LocalPopulations[i].Pdad7[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom8[j] = ThisMutation(LocalPopulations[i].Pmom8[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad8[j] = ThisMutation(LocalPopulations[i].Pdad8[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom9[j] = ThisMutation(LocalPopulations[i].Pmom9[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad9[j] = ThisMutation(LocalPopulations[i].Pdad9[j], sdmu, -1, 1);
            
            LocalPopulations[i].Pmom10[j] = ThisMutation(LocalPopulations[i].Pmom10[j], sdmu, -1, 1);
            LocalPopulations[i].Pdad10[j] = ThisMutation(LocalPopulations[i].Pdad10[j], sdmu, -1, 1);
        }
    }
}

// this function runs everything...
int main(int argc, char* argv[])
{
     
    // create and open all output files
    ofstream DataFile, DistFile, ParametersFile;
    
    string fname = argv[1];
    
    /*
    string fname1 = fname;
    fname1.append("_Dist.txt");
    char *fileName1 = (char*)fname1.c_str();
    DistFile.open(fileName1);
    */
    
    string fname2 = fname;
    fname2.append("_Parameters.txt");
    char *fileName2 = (char*)fname2.c_str();
    ParametersFile.open(fileName2);
    
    
    
    //--------------------------------------------------//
    
    //read NB structure of world
       fstream myfilenb("nbmat_matMT_80000000000_justnb_closestland.txt", std::ios_base::in);
       float nba;
       int countnb=0; //celID
       int nbid = -1; //nbid
       while (myfilenb >> nba)
       {
        //printf("%f ", nba);
         nbid++;
         World[countnb].Neighbor_ids[nbid] = nba;
         if(nbid==6) {countnb++; // nextcell, 6 since nbid is initialized at -1
                      nbid=-1; //reset nbid
                     }
        }
    //int nbtest=51;
    //cout<<endl<<World[nbtest].Neighbor_ids[0]<<" "<<World[nbtest].Neighbor_ids[1]<<" "<<World[nbtest].Neighbor_ids[2]<<" "<<
    //World[nbtest].Neighbor_ids[3]<<" "<<World[nbtest].Neighbor_ids[4]<<" "<<World[nbtest].Neighbor_ids[5]<<" ";
    
    //WriteHeaders(DistFile, ParametersFile);
    WriteHeaders_Pfile(ParametersFile);

       
    // read Ebar and sd_e of World
    fstream myfile("MTmatrix_80000000000_scaledCB_999.txt", std::ios_base::in);
     //myfile.open("Trial72MTmatrix.txt");
     float a;
     int count=0;
     while (myfile >> a)
     {
      //printf("%f ", a);
       World[count].Ebar = a;
       count++;
     }

    //myfile.close();
     
     // it works:
     /*for(int lala=0; lala<72;lala++)
        {cout<<"Ebarvali"<<lala<<"  "<<World[lala].Ebar<<" SDEvali "<<World[lala].sd_E<<endl;
         
     } */
     
     // same for sd_E
     fstream myfile2("VTmatrix_80000000000_scaledCB_999.txt", std::ios_base::in);
     float a2;
     int count2=0;
     while (myfile2 >> a2)
     {
       World[count2].sd_E = a2;
       count2++;
     }

    
    for (Replicate = 1; Replicate < Runs + 1; Replicate++)
    {
        Init(ParametersFile);
        
        Generation = 0;
        
        // First do the simulations with fixed population size to determine what to expect
        for (Generation = 1; Generation <= NumGen; Generation++)
        {
            // simulate selection events
            OneGeneration();
           
            // make sure that there not every lineage has gone extinct
            int TotInds = 0;
            for (int i = 0; i < Ncell; i++) {
                TotInds += LocalPopulations[i].PopSize;
            }
            
            // generate some screen output to assess progress
            if (Generation==1 | Generation==2 |Generation==5 |Generation==10 |Generation==20 |Generation==30 |Generation==40 |Generation==50 |Generation==60 |Generation==70 |Generation==80 |Generation==90 |Generation==100 |Generation==200 |Generation%1000==0) {
                
                cout << endl << endl << "***** "
                << "Run: " << Replicate << '\t' << "Generation: " << Generation << '\t' << "WorldPop"<< '\t' << TotInds << endl << endl;

            string fname1 = fname;
            fname1.append("_");
            //if (Generation==NumGen) { //for now - print the last one - much faster
                fname1 = fname1 + std::to_string(Generation);
                fname1.append("_Dist.txt");
                char *fileName1 = (char*)fname1.c_str();
                DistFile.open(fileName1);
                WriteHeaders_Dfile(DistFile);
                WriteDist(DistFile);
                DistFile.close();
              //}
            }

            if (TotInds==0) { break; }
        }
    }
  
     
    // close all output files
    //DistFile.close();
    ParametersFile.close();
   
    // enjoy!!
    
 }



