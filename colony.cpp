/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by Jan J. Kreider, Ido Pen, Boris H. Kramer and G. Sander van Doorn
 */


#include <cassert>
#include <cmath>
#include "colony.h"
#include "random.h"

std::vector<Male> Colony::drones; // a vector with objects of the class Male

Colony::Colony(Queen * const q) : queen(q), colonyFitness(0.0) // constructor for Colony. Fitness is set to 0. Colonies have a queen
{
    
}

Colony::~Colony() // destructor for Colony
{
    // delete queen, workers and offspring
	if(queen != nullptr) delete queen;
	while (!workers.empty()) {
		delete workers.back();
		workers.pop_back();
    }
    while (!offspring.empty()) {
        delete offspring.back();
        offspring.pop_back();
    }
}

void Colony::reproduce()
{
    double workerMortality;
    if (isQueenless() && !sterileWorkers) { // if there is no queen and if workers are fecund
        assert(workers.size());
        for (int i = 0u; i < workers.size(); ++i) { // loop through workers
            double tempFecundity = workers[i]->getLifeHistoryTrait(Individual::fecundity); // get fecundity value of every worker
            double optimalWorkerAllocation = 1.0 / (1.0 + sqrt(tempFecundity / productivity)); // proportion of resources allocated to reproduction
            workerMortality = 1.0 - pow(1.0 - foragingRisk, 1.0 - optimalWorkerAllocation); // worker mortality rate, pow(base ,exponent)
            const double eggs = tempFecundity * optimalWorkerAllocation; // number of eggs
            const double resources = productivity * (1.0 - optimalWorkerAllocation); // resources foraged
            double workerFitness = eggs * resources / (eggs + resources); // calculate fitness
            for (size_t h = rnd::poisson(workerFitness); h > 0u; --h) { // a vector with a poisson distribution, lambda = fitness. this introduces extra stochasticity. vector is likely to have a similar size (this is an int) like fitness (this is a double)
                drones.push_back(Male(workers[i])); // Male objects are addec to vector drones
            }
            workerFitness = 0; // set this back to 0 for next worker, for safety
            workers[i]->extrinsicMortality(workerMortality); // every worker gets an extrinsic mortality
        }
    }
    if (!isQueenless()) { // if there is a queen
        double tempFecundity = queen->getLifeHistoryTrait(Individual::fecundity); // get queen fecundity
		double optimalQueenAllocation = size() / (1.0 + sqrt(tempFecundity / productivity)); // allocation depends on size of colony. optimalQueenAllocation is larger when there are more individuals in the colony
		if (optimalQueenAllocation > 1.0) // if optimalQueenAllocation is larger then 1
			optimalQueenAllocation = 1.0; // it is set to 1
		queen->extrinsicMortality(1.0 - pow(1.0 - foragingRisk, 1.0 - optimalQueenAllocation)); // calculate queen death rate. function delivers extrinsic survival
		workerMortality = foragingRisk;
		const double eggs = tempFecundity * optimalQueenAllocation; // number of eggs produced per queen
		const double resources = productivity * (1.0 - optimalQueenAllocation + workers.size()); // foraged resources
		colonyFitness = eggs * resources / (eggs + resources); // calculate fitness
		for (size_t i = rnd::poisson(colonyFitness); i > 0u; --i) { // a vector with a poisson distribution, lambda = fitness. this introduces extra stochasticity. vector is likely to have a similar size (this is an int) like fitness (this is a double)
            Female *tmp = new Female(queen); // new Females are produced by queen
			tmp->extrinsicMortality(0.0); // new Females have extrinsicMortality 0. They do not leave the colony for foraging when they are still offspring
			offspring.push_back(tmp); // new Females are written to the end of the offspring vector
			drones.push_back(Male(queen)); // for every Female the queen also produces one Male which is written to the drone vector
        }
        for (Female * worker : workers) // For females that become members of vector workers
            worker->extrinsicMortality(workerMortality); // extrinsicMortality is set to standard settings from parameter file
    }
}

void Colony::survival()
{
    // survival for workers
	for (size_t i = 0u; i < workers.size(); ) { // loop through workers
		workers[i]->incrementAge(); // increase age by one age class
		if (workers[i]->doesSurvive()) ++i; // if individuals survive, continue loop
		else { // if they die
			delete workers[i]; // delete the respective worker
			workers[i] = workers.back(); // move dead worker to the end of the vector
			workers.pop_back(); // delete vector entry of dead worker
		}
	}
   
    // survival if colony has a queen
	if (!isQueenless()) { // if colony is not queenless
		queen->incrementAge(); // increase age of queen
		if (!queen->doesSurvive()) { // if queen does not survive
			delete queen; // delete queen
			queen = nullptr;
		}
	}

    if (isQueenless()) { // if there is no queen
        if (gyneReplacesQueen && workerReplacesQueen) { // if both colony turnovers are true
            while (isQueenless() && offspring.size() && workers.size() && drones.size()) { // as long as the colony is still queenless, there is offspring, drones, and workers
                if (rnd::uniform() < 0.5) { // 50 % chance for each serial polygyny and colony inheritance
                    queen = new Queen(offspring.back(), getDrone()); // new queen is last data entry in vector offspring. get drone to mate with
                    if (!queen->doesSurvive()) { // if queen does not survive
                        delete queen;  // delete queen
                        queen = nullptr;
                    }
                    delete offspring.back(); // delete offspring which is now the queen
                    offspring.pop_back(); // delete last data entry in offspring
                } else {
                    size_t i = rnd::integer(workers.size()); // i is a random worker from vector workers
                    queen = new Queen(workers[i], getDrone()); // worker i becomes the new queen mated with a random drone
                    delete workers[i]; // worker i is deleted
                    workers[i] = workers.back(); // move worker i to end of vector workers
                    workers.pop_back(); // delete last entry in vector workers
                }
            }
            
            while (isQueenless() && offspring.size() && !workers.size() && drones.size()) { // if there are no workers, then it is automatically serial polygyny
                queen = new Queen(offspring.back(), getDrone()); // new queen is last data entry in vector offspring. get drone to mate with
                if (!queen->doesSurvive()) { // if queen does not survive
                    delete queen;  // delete queen
                    queen = nullptr;
                }
                delete offspring.back(); // delete last data entry in offspring, where queen was stored
                offspring.pop_back();
            }
            
            if (isQueenless() && !offspring.size() && workers.size() && drones.size()) { // if there is no offspring, then it is automatically colony intertance
                size_t i = rnd::integer(workers.size()); // i is a random worker from vector workers
                queen = new Queen(workers[i], getDrone()); // worker i becomes the new queen mated with a random drone
                delete workers[i]; // worker i is deleted
                workers[i] = workers.back(); // move worker i to end of vector workers
                workers.pop_back(); // delete last entry in vector workers
            }
        }
        if (gyneReplacesQueen && !workerReplacesQueen) { // if only serial polygyny
            while (isQueenless() && offspring.size() && drones.size()) {
                queen = new Queen(offspring.back(), getDrone()); // new queen is last data entry in vector offspring. get drone to mate with
                if (!queen->doesSurvive()) { // if queen does not survive
                    delete queen;  // delete queen
                    queen = nullptr;
                }
                delete offspring.back(); // delete last data entry in offspring, where queen was stored
                offspring.pop_back();
            }
        }
        if (!gyneReplacesQueen && workerReplacesQueen) { // if only colony inheritance
            if (workers.size() && drones.size()) {
                size_t i = rnd::integer(workers.size()); // i is a random worker from vector workers
                queen = new Queen(workers[i], getDrone()); // worker i becomes the new queen mated with a random drone
                delete workers[i]; // worker i is deleted
                workers[i] = workers.back(); // move worker i to end of vector workers
                workers.pop_back(); // delete last entry in vector workers
            }
        }
    }
}


void Colony::addWorkers() // transfers objects from vector offspring to vector workers
{
	while (!offspring.empty()) { // while vector offspring is not empty
        if (offspring.back()->doesSurvive()) { // if last element of vector offspring survive
			workers.push_back(offspring.back()); // put last element of vector offspring to worker vector
        } else { // if last element of vector offspring dies
			delete offspring.back(); // delete last element of vector offspring
        }
		offspring.pop_back(); // reduce container position
	}
}


Female* Colony::getGyne() // function selects a new queen
{
	assert(offspring.size());
	Female *tmp = offspring.back(); // tmp is a pointer to the last object in vector offspring
	offspring.pop_back(); // reduce container position
	return tmp; // return last object of vector offspring to be the new queen
}

Male Colony::getDrone() // function selects a drone to mate the queen
{
    assert(drones.size());
    size_t tmp = rnd::integer(drones.size()); // choose a drone from the drone vector
    Male tmpDrone = drones[tmp]; // copy it
    drones.erase(drones.begin() + tmp); // delete it from the vector
    return tmpDrone; // return the copy
}

void Colony::getData(size_t &sumn, size_t& sumq, size_t& sumw, size_t& suma, std::vector<double> &sumxq, std::vector<double> &sumxw, std::vector<double> &sumfeqq, std::vector<double> &sumfeqw) const // get data from colonies
{
	if (size()) { // if there are individuals in the colony
		++sumn;	// increase sumn
		if (!isQueenless()) { // if colony is not queenless
			++sumq; // increase sumq
            suma += queen->getAge(); // sum of ages of queens
			queen->getData(sumxq, sumfeqq); // write life history table of queens, both traits
		}
		sumw += workers.size(); // record number of workers
        for (Female * worker : workers) {
            worker->getData(sumxw, sumfeqw); // write life history tables, both traits
        }
    }
}
