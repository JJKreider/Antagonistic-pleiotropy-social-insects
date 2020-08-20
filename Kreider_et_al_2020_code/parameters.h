/*
 Antagonistic pleiotropy and the evolution of eusociality
 Jan J. Kreider, Ido Pen, Boris H. Kramer
 */

#ifndef parameters_h
#define parameters_h

// general settings
const int nColony = 1000; // number of colonies
const int tEnd = 500000; // simulation time
const int dataInterval = 100; // data is saved tEnd/dataInterval times

// life history
const int ageMax = 20; // maximum lifespan
const double foragingRisk = 0.0; // extrinsic mortality
const double survivalOffset	= 3.2; // initial allele value for survival
const double fecundityOffset = 0.5; // initial allele value for fecundity
const double maxFecundity = 5.0; // logistic parameter, restricting the increase of fecundity over evolutionary time
const double productivity = 1.0; // colony productivity

// simulation scenarios
const bool gyneReplacesQueen = false; // monogyny = false, serial polygyny = true
const bool workerReplacesQueen = false; // colony inheritance
const bool sterileWorkers = true; // if workers are sterile they are killed when the queen dies, if they are fecund they lay eggs after the queen died

// genetics
const int nGene = 2; // number of genes
const int nCaste = 2; // caste 0 = queens, caste 1 = workers
const int nTrait = 2; // survival, fecundity
const double mutationRate = 0.03; // mutation rate

// partial correlation matrix parameters
const double alpha = 0.0; // adjacent age classes, effect per trait
const double beta = -0.0; // within caste within trait effects
const double ggamma = -0.0; // within caste between trait effects
const double delta = -0.0; // between caste within trait effects
const double eta = -0.0; // between caste between trait effects
const int delay = 0; // delay of antagonistic effects

// mutational variance covariance matrix parameters
const double sigma = 0.2; // mutational effect size over all age classes
const double lambda = -0.2; // mutation bias
#endif
