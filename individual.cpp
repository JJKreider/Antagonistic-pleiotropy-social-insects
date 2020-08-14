/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by Jan J. Kreider, Ido Pen, Boris H. Kramer and G. Sander van Doorn
 */


#include "individual.h"
#include "random.h"
#include <cassert>
#include <iostream>
#include <fstream>
#include <sstream>


std::array<size_t, ageMax + 1u> Female::ageOfDeath; // death ages of Females
std::array<size_t, ageMax + 1u> Queen::ageOfDeath; // death ages of Queens

rnd::MultiNormal const * Gene::geneticArchitecture; // initialise geneticArchitecture with a multinormal distribution

extern unsigned int simulationId;

void Gene::initGeneticArchitecture()
{
	Matrix P(K, K, 0.0); //  partial correlation matrix P with K*K values of 0.0, K = nCaste * nTrait * ageMax

    // diagonal of the matrix
    for (int tr = 0; tr < nTrait; ++tr) {
        for (int cs = 0; cs < nCaste; ++cs) {
            for (int i = 0; i < ageMax; ++i) {
                const int k1 = i + (cs + tr * nCaste) * ageMax;
                P(k1, k1) = 1.0;
                
                // adjacent age classes
                if (i - 1 >= 0) P(k1, k1 - 1) = alpha;
                if (i + 1 < ageMax) P(k1, k1 + 1) = alpha;
            }
        }
    }
    
    // within caste within trait effects
    if (beta < 0) {
        for (int rowQ = 0; rowQ < ageMax; ++rowQ) {
            if (rowQ + delay < ageMax) P(rowQ, rowQ + delay) = beta;
            if (rowQ + delay + 1 < ageMax) P(rowQ, rowQ + delay + 1) = alpha * beta;
            if ((delay > 2) && (rowQ + delay - 1 < ageMax)) P(rowQ, rowQ + delay - 1) = alpha * beta;
            if (rowQ + delay < ageMax) P(rowQ + nTrait * ageMax, rowQ + delay + nTrait * ageMax) = beta;
            if (rowQ + delay + 1 < ageMax) P(rowQ + nTrait * ageMax, rowQ + delay + 1 + nTrait * ageMax) = alpha * beta;
            if ((delay > 2) && (rowQ + delay - 1 < ageMax)) P(rowQ + nTrait * ageMax, rowQ + delay - 1 + nTrait * ageMax) = alpha * beta;
        }
        for (int rowW = ageMax; rowW < (ageMax * 2); ++rowW) {
            if (rowW + delay < (ageMax * 2)) P(rowW, rowW + delay) = beta;
            if (rowW + delay + 1 < (ageMax * 2)) P(rowW, rowW + delay + 1) = alpha * beta;
            if ((delay > 2) && (rowW + delay - 1 < (ageMax * 2))) P(rowW, rowW + delay - 1) = alpha * beta;
            if (rowW + delay < (ageMax * 2)) P(rowW + nTrait * ageMax, rowW + delay + nTrait * ageMax) = beta;
            if (rowW + delay + 1 < (ageMax * 2)) P(rowW + nTrait * ageMax, rowW + delay + 1 + nTrait * ageMax) = alpha * beta;
            if ((delay > 2) && (rowW + delay - 1 < (ageMax * 2))) P(rowW + nTrait * ageMax, rowW + delay - 1 + nTrait * ageMax) = alpha * beta;
        }
    }
    
    // within caste between trait effects
    for (int tr1 = 0; tr1 < nTrait; ++tr1) {
        for (int cs = 0; cs < nCaste; ++cs) {
            for (int i = 0; i < ageMax; ++i) {
                const int k1 = i + (cs + tr1 * nCaste) * ageMax;
                if (!sterileWorkers) {
                    for (int tr2 = tr1 + 1; tr2 < nTrait; ++tr2) {
                        const int k2 = i + (cs + tr2 * nTrait) * ageMax;
                        if (i + delay < ageMax) P(k1, k2 + delay) = ggamma;
                        if (i + delay + 1 < ageMax) P(k1, k2 + delay + 1) = alpha * ggamma;
                        if (i + delay - 1 < ageMax and i + delay - 1 >= 0) P(k1, k2 + delay - 1) = alpha * ggamma;
                        if (i + delay < ageMax) P(k2, k1 + delay) = ggamma;
                        if (i + delay + 1 < ageMax) P(k2, k1 + delay + 1) = alpha * ggamma;
                        if (i + delay - 1 < ageMax and i + delay - 1 >= 0) P(k2, k1 + delay - 1) = alpha * ggamma;
                    }
                } else {
                    for (int tr2 = tr1 + 1; tr2 < nTrait; ++tr2) {
                        const int k2 = i + (cs + tr2 * nTrait) * ageMax;
                        if (k2 < ageMax * 3) {
                            if (i + delay < ageMax) P(k2, k1 + delay) = ggamma;
                            if (i + delay + 1 < ageMax) P(k2, k1 + delay + 1) = alpha * ggamma;
                            if (i + delay - 1 < ageMax and i + delay - 1 >= 0) P(k2, k1 + delay - 1) = alpha * ggamma;
                            if (i + delay < ageMax) P(k1, k2 + delay) = ggamma;
                            if (i + delay + 1 < ageMax) P(k1, k2 + delay + 1) = alpha * ggamma;
                            if (i + delay - 1 < ageMax and i + delay - 1 >= 0) P(k1, k2 + delay - 1) = alpha * ggamma;
                        }
                    }
                }
            }
        }
    }
    
    // between caste within trait effects
    for (int tr = 0; tr < nTrait; ++tr) {
        for (int cs = 0; cs < nCaste; ++cs) {
            for (int i = 0; i < ageMax; ++i) {
                const int k1 = i + (cs + tr * nCaste) * ageMax;
                if (k1 + ageMax < K) {
                    if (!sterileWorkers) {
                        P(k1, k1 + ageMax) = P(k1 + ageMax, k1) = delta;
                        if (i - 1 >= 0) P(k1 - 1, k1 + ageMax) = P(k1 + ageMax, k1 - 1) = delta * alpha;
                        if (i + 1 < ageMax) P(k1 + 1, k1 + ageMax) = P(k1 + ageMax, k1 + 1) = delta * alpha;
                    } else {
                        if (k1 < nCaste * ageMax) {
                            P(k1, k1 + ageMax) = P(k1 + ageMax, k1) = delta;
                            if (i - 1 >= 0) P(k1 - 1, k1 + ageMax) = P(k1 + ageMax, k1 - 1) = delta * alpha;
                            if (i + 1 < ageMax) P(k1 + 1, k1 + ageMax) = P(k1 + ageMax, k1 + 1) = delta * alpha;
                        }
                    }
                }
            }
        }
    }
    
    // between caste between trait effects
    for (int tr = 0; tr < nTrait; ++tr) {
        for (int cs = 0; cs < nCaste; ++cs) {
            for (int i = 0; i < ageMax; ++i) {
                const int k1 = i + (tr + cs * nCaste) * ageMax;
                for (int tr2 = tr + 1; tr2 < nTrait; ++tr2) {
                    const int k2 = i + (tr2 + cs * nCaste) * ageMax;
                     if (k2 < 3 * ageMax) {
                         P(k1 + 2 * ageMax, k2) = eta * P(k1, k1);
                         if (i - 1 >= 0) P(k1 + 2 * ageMax, k2 - 1) = eta * P(k1, k1 - 1);
                         if (i + 1 < ageMax) P(k1 + 2 * ageMax, k2 + 1) = eta * P(k1, k1 + 1);
                         P(k2, k1 + 2 * ageMax) = eta * P(k1, k1);
                         if (i - 1 >= 0) P(k2 - 1, k1 + 2 * ageMax) = eta * P(k1, k1 - 1);
                         if (i + 1 < ageMax) P(k2 + 1, k1 + 2 * ageMax)  = eta * P(k1, k1 + 1);
                         
                         // if workers are not sterile also worker fecundity can be affected by correlations
                         if (!sterileWorkers) {
                             P(k2 + 2 * ageMax, k1) = eta * P(k1, k1);
                             if (i - 1 >= 0) P(k2 - 1 + 2 * ageMax, k1) = eta * P(k1, k1 - 1);
                             if (i + 1 < ageMax) P(k2 + 1 + 2 * ageMax, k1) = eta * P(k1, k1 + 1);
                             P(k1, k2 + 2 * ageMax) = eta * P(k1, k1);
                             if (i - 1 >= 0) P(k1, k2 - 1 + 2 * ageMax) = eta * P(k1, k1 - 1);
                             if (i + 1 < ageMax) P(k1, k2 + 1 + 2 * ageMax) = eta * P(k1, k1 + 1);
                         }
                    } 
                }
            }
        }
    }
    
	std::ostringstream oss; // output stream
	oss << "architecture_" << simulationId << ".csv";  // write file name including simulationID as .csv file
	std::ofstream ofs(oss.str().c_str()); // write output stream in file
	verify(ofs.is_open()); // verify if ofs is opened
	ofs.fill(','); // separate data with commas
	ofs << "partial correlation matrix\n\n" << P << '\n'; // write header partial correlation matrix and the matrix P to the file

	// calculation of mutational variance covariance matrix
    
    // first make all values apart from the diagonal negative
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            P(i, j) = -P(i, j);
        }
        P(i, i) = -P(i, i);
    }
    
    // inverse the matrix
    P = P.inverse();
    
    // apply the tranformation equation
    Matrix Q(K, K, 0.0);
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j < K; ++j) {
            Q(i, j) = P(i, j) / sqrt(P(i, i) * P(j, j));
        }
    }
    
    // create a matrix with the stddev
    Matrix Sig(K, K, 0.0);
    
    for (int i = 0; i < K; ++i) {
        Sig(i, i) = sigma;
    }
   
    // final matrix transformation step
    P = Sig * Q * Sig;
    
    // mutation bias vector
    Vector bias(K);
    for (int i = 0; i < K; ++i) {
        bias[i] = lambda * sigma;
    }

	ofs << "\nmutational variance covariance matrix\n\n" << P << '\n'; // output matrix
    ofs.close();
    
	geneticArchitecture = new rnd::MultiNormal(P, bias); // initialize geneticArchitecture with a multinormal distribution with the new matrix and mutation bias
   
}

Male::Male(Individual const * const mother) // Male inherits genes from mother
{
    for (int i = 0; i < nGene; i += 2) { // loop through genes but in steps of 2, so every second gene
        genome[i / 2] = rnd::uniform() < 0.5 ? mother->genome[i] : mother->genome[i + 1]; // choose one of the diploid genes of mother
    }
    mutate(); // mutate the genes
}

Female::Female(Queen const * const mother) // Female inherits genes from mother
{
	for (int i = 0; i < nGene; i += 2) {  // loop through pairs of genes
		genome[i] = rnd::uniform() < 0.5 ? mother->genome[i] : mother->genome[i + 1]; // take one gene of the diploid gene pair
		genome[i + 1] = mother->sperm.genome[i / 2]; // second part of genome comes from sperm of drones
    }
    mutate(); // mutate the genes
	develop(); // calculate life history values
}


Queen::Queen(Female const * const female, Male drone) : Individual(*female), sperm(drone) // queens are made from a female and carry sperm from a male
{
	age = 0; // new queens, i.e. also when a worker becomes a queen in the colony inheritance scenario, start with age 0
	develop(); // calculate life history values
}

void Individual::mutate()
{
	if (rnd::uniform() < nGene * mutationRate) // check if a mutation occurs
        genome[rnd::integer(nGene)].mutate(); // sample a gene from the genome and increase dx by the values sampled from the multinormal distribution
}

void Male::mutate()
{
	if (rnd::uniform() < nGene * mutationRate / 2) // check if a mutation occurs
		genome[rnd::integer(nGene / 2)].mutate(); // sample a gene from the genome and increase dx by the values sampled from the multinormal distribution
}

void Female::develop() // calculate life history values
{
	for (int a = 0; a < ageMax; ++a) {
        lifeHistory(survival, a) = survivalOffset; // initialize survival values
        lifeHistory(fecundity, a) = fecundityOffset; // initialize fecundity values
        for (int i = 0; i < nGene; ++i) {
            lifeHistory(survival, a) += genome[i](1, survival, a); // add mutation value to respective age class
            lifeHistory(fecundity, a) += genome[i](1, fecundity, a); // add mutation value to respective age class
        }
        lifeHistory(survival, a) = 1 / (1 + exp(-lifeHistory(survival, a))); // transform genotypic to phenotypic survival probability
        lifeHistory(fecundity, a) = maxFecundity / (1 + exp(-lifeHistory(fecundity, a))); // transfrom genotypic to phenotypic fecundity value, limited by maxFecundity parameter
    }
}

void Queen::develop() // calculate life history values
{
	for (int a = 0; a < ageMax; ++a) {
		lifeHistory(survival, a) = survivalOffset; // initialize survival values
        lifeHistory(fecundity, a) = fecundityOffset; // initialize survival values
        for (int i = 0; i < nGene; ++i) {
			lifeHistory(survival, a) += genome[i](0, survival, a); // add mutation value to respective age class
            lifeHistory(fecundity, a) += genome[i](0, fecundity, a); // add mutation value to respective age class
        }
		lifeHistory(survival, a) = 1 / (1 + exp(-lifeHistory(survival, a))); // transform genotypic to phenotypic survival probability
        lifeHistory(fecundity, a) = maxFecundity / (1 + exp(-lifeHistory(fecundity, a))); // transfrom genotypic to phenotypic fecundity value, limited by maxFecundity parameter
    }
}

std::ostream& operator<< (std::ostream &os, const Individual &obj) // output stream for Individual
{
    os << obj.lifeHistory << '\n'; // output lifeHistory table
    return os;
}
 
void Individual::getData(std::vector<double>& sur, std::vector<double>& fec) const
{
    for (int a = 0; a < ageMax; ++a) {
        sur[a] += lifeHistory(survival, a); // write survival probabilities in vector
        fec[a] += lifeHistory(fecundity, a); // write fecundity values in vector
    }
}
