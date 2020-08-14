/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by Jan J. Kreider, Ido Pen, Boris H. Kramer and G. Sander van Doorn
 */

#include "colony.h"
#include "random.h"
#include "utils.h"
#include "parameters.h"
#include <numeric>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <iterator>
#include <fstream>
#include <sstream> 

unsigned int simulationId;
Colony * population[nColony]; //an array population of length nColonies with pointers to Colonies
std::ofstream ofs; // output stream

std::ofstream ofParameters; // output stream

bool is_regular_file(const std::string& filename) // a function to check if an output file exists and is opened
{
	std::fstream f; // output/input stream for files
	f.open(filename.c_str(), std::ios::in); // open file
	return f.is_open(); // returns whether the stream is associated to a file
}

void init() // initialize the simulation
{
	simulationId = static_cast<unsigned int>(std::chrono::high_resolution_clock::now().time_since_epoch().count()); // produce a simulation Id
	rnd::set_seed(simulationId); // use simulation Id to set seed for random number generator

	Gene::initGeneticArchitecture(); // create mutational variance covariance matrix

	Female *q = new Female(); // construct a Female
    for (int i = 0; i < nColony; ++i) { // loop through all colonies
        Queen *tmp = new Queen(q, Male(q)); // every Colony in the population array gets a mated queen with the genetics of the Female q
        population[i] = new Colony(tmp);
    }
    delete q;

	for (int k = 0; k <= ageMax; ++k) // loop through ages
		Queen::ageOfDeath[k] = Female::ageOfDeath[k] = 0u; // initialize arrays ageOfDeath with 0
	
	std::ostringstream str; // output stream str
	str << "data_" << simulationId << ".csv"; // str has a filename data_simulationID.csv
	ofs.open(str.str()); // open the file
	ofs << "t" << ',' // headers for csv file
		<< "proportion occupied" << ','
		<< "proportion queenright" << ','
		<< "avg number of workers/colony" << ','
		<< "avg queen age" << ','
		<< "avg queen lifespan" << ','
        << "avg worker lifespan";
    for (int a = 0; a < ageMax; ++a)
		ofs << ',' << "queen intrinsic survival age class " << a + 1;
	for (int a = 0; a < ageMax; ++a)
		ofs << ',' << "worker intrisic survival age class " << a + 1;
    for (int a = 0; a < ageMax; ++a)
        ofs << ',' << "queen fecundity age class " << a + 1;
    for (int a = 0; a < ageMax; ++a)
        ofs << ',' << "worker fecundity age class " << a + 1;
    ofs << ',' << "avg queen fecundity"
        << ',' << "avg worker fecundity";
    ofs << "\n"; // finish writting
}

void dataAnalysis(int t) // data analysis function
{
	std::vector<double> sumxq(ageMax, 0.0), sumxw(ageMax, 0.0), sumfeqq(ageMax, 0.0), sumfeqw(ageMax, 0.0); // vectors of length ageMax with 0.0 in them
    size_t sumw = 0u, sumq = 0u, sumn = 0u, suma = 0u;

    for (int i = 0; i < nColony; ++i) {// loop through all colonies
		population[i]->getData(sumn, sumq, sumw, suma, sumxq, sumxw, sumfeqq, sumfeqw); // getData from every i in population where i is a colony. sumn = number of occupied colonies, sumq = total number of queens, sumw = total number of workers, suma = sum of ages of queens, sumxq = life history tables of queens survival, sumxw = life history tables of workers survival, sumfeqq = life history table of queens fecundity, sumfeqw = life history tables of workers fecundity
    }

	ofs << t << ',' // write into file ofs
		<< sumn * 1.0 / nColony << ',' // ratio of occupied colonies
		<< sumq * 1.0 / sumn << ',' // ratio of queenright colonies
		<< sumw * 1.0 / sumn << ',' // average number of workers per colony
		<< (sumq ? suma * 1.0 / sumq : nan("")) << ','; // if there are any queens calculate their average age and if there are none print nothing

	sumn = suma = 0u; // set sumn and suma to 0 again
	for (int k = 0; k <= ageMax; ++k) { // loop through ages
		double n = Queen::ageOfDeath[k]; // take all ages of death of dead queens
		sumn += n; // sum the number of dead queens up
		suma += k * n; // sum of the products of the age and the number of deaths in that age class, i.e. death ages
		Queen::ageOfDeath[k] = 0u; // set the value in the age of death array to 0 after it was read
	}
	ofs << suma * 1.0 / sumn << ','; // average queen life span

	sumn = suma = 0u; // set sumn and suma to 0 again
	for (int k = 0; k <= ageMax; ++k) { // loop through ages
		double n = Female::ageOfDeath[k]; // take all ages of death of dead workers
		sumn += n; // sum the number of dead workers up
		suma += k * n; // sum of the products of the age and the number of deaths in that age class, i.e. death ages
		Female::ageOfDeath[k] = 0u; // set the value in the age of death array to 0 after it was read
	}
	ofs << suma * 1.0 / sumn; // average worker life span

	for (int a = 0; a < ageMax; ++a) // loop through ages
		ofs << ',' << (sumq ? sumxq[a] / sumq : nan("")); // if there are queens divide sum of age specific survival probabilities by total number of queens else print nothing. this is queen intrinsic survival
	for (int a = 0; a < ageMax; ++a) // loop through ages
		ofs << ',' << (sumw ? sumxw[a] / sumw : nan("")); // see above but for workers
    for (int a = 0; a < ageMax; ++a)
        ofs << ',' << (sumq ? sumfeqq[a] / sumq : nan("")); // same for fecundity for queens
    for (int a = 0; a < ageMax; ++a)
        ofs << ',' << (sumw ? sumfeqw[a] / sumw : nan("")); // and for workers

    ofs << ',';
    sumn = 0u; // set sumn to 0 again
    double sumf = 0; // this is a local variable to store the sum of all fecundities of all individuals, separately for the castes
    for (int i = 0; i < nColony; ++i) { // loop throught all colonies
        if (population[i]->isQueenless() == false) { // if there is a queen
            ++sumn; // sumn will be number of queenright colonies when loop is finished
            sumf += population[i]->getQueen()->getLifeHistoryTrait(Queen::fecundity); // this is then the sum of all fecundities of all queens
        }
    }

    ofs << sumf * 1.0 / sumn << ','; // average queen fecundity in the population
    
    sumn = 0u; // set this to 0 again
    sumf = 0; // same here
    for (int i = 0; i < nColony; ++i) { // loop through all colonies
        if (population[i]->isQueenless()) { // if the population is queenless
            sumn += population[i]->size(); // add number of workers
        } else {
            sumn += population[i]->size() - 1; // if there is a queen it has to be size - 1 because size was number of workers + queen
        }
        std::vector<Female*> tmpWorkers = population[i]->getWorkers(); // store the workers locally
        for (int j = 0; j < tmpWorkers.size(); ++j) { // loop through the workers
            sumf += tmpWorkers[j]->getLifeHistoryTrait(Female::fecundity); // and get their fecundities summed up
        }
        tmpWorkers.clear(); // clear the local vector
    }
    ofs << sumf * 1.0 / sumn; // calculate average worker fecundity
    ofs << '\n';
    
	ofs.flush(); // synchronize stream buffer with output sequence
}


int main()
{
	init(); // initialize simulation
    for (int t = 0; t < tEnd; ) { // loop trough time steps
		size_t sum = 0u;
		std::vector<size_t> weights(nColony, 0u); // a vector of length nColony with unsigned 0
        
        // this loop is separate from the next one because drones that are produced in colony i are not available in colony i+1, if i becomes queenless
        for (int i = 0; i < nColony; ++i) { // loop through colonies
			if (population[i]->size()) { // if there are individuals in the colony
				population[i]->reproduce(); // colonies reproduce
            }
		}
        
        for (int i = 0; i < nColony; ++i) { // loop through colonies
            if (population[i]->size()) { // if there are individuals in the colony
                population[i]->survival(); // check survival, get new queens if colony inheritance or serial polygyny
                sum += weights[i] = population[i]->nrOffspring(); // write the number of offspring in vector weights and sum up total number of offspring
            }
        }
        
        if (sterileWorkers) { // if workers shall be killed when queen dies, i.e. workers are sterile
            for (int i = 0; i < nColony && sum; ) { // interate i as long as it is smaller than nColony and smaller than sum
                if (population[i]->isQueenless()) { // if colony is queenless
                    size_t j = std::discrete_distribution<size_t>(weights.begin(), weights.end())(rnd::rng); // sample a colony according to its fitness
                    if (population[j]->nrOffspring() > 0) { // if there is offspring in the sampled colony
                        Female * tmp = population[j]->getGyne(); // take a female to be the new queen
                        Queen * q = new Queen(tmp, Colony::getDrone()); // mate the female with a drone and make it the queen
                        delete tmp; // delete pointer
                        if (q->doesSurvive()) {
                            delete population[i];
                            population[i] = new Colony(q);
                            ++i;
                        } else {
                            delete q;
                        }
                        --weights[j];
                        --sum;
                    }
                }
                else ++i; // if there is a queen
            }
        } else { // if workers live on when queen is dead, i.e. they are reproductive
            for (int i = 0; i < nColony && sum; ) { // interate i as long as it is smaller than nColony and smaller than sum
                if (population[i]->size() == 0u && weights[i] == 0u) { // if there are no individuals in a colony and also no offspring
                    size_t j = std::discrete_distribution<size_t>(weights.begin(), weights.end())(rnd::rng); // sample a colony according to its fitness
                    Female * tmp = population[j]->getGyne(); // take a female to be the new queen
                    Queen * q = new Queen(tmp, Colony::getDrone()); // mate the female with a ddrone and make it the queen
                    delete tmp; // delete pointer
                    if (q->doesSurvive()) {
                        delete population[i];
                        population[i] = new Colony(q);
                        ++i;
                    } else {
                        delete q;
                    }
                    --weights[j];
                    --sum;
                }
                else ++i; // if there are individuals or offspring
            }
        }

		bool terminate = true; // create a bool terminate
		for (int i = 0; i < nColony; ++i) { // loop through colonies
			population[i]->addWorkers(); // addWorkers which were offspring before to colony i
			if (population[i]->size()) // if there are individuals in any colony
				terminate = false; // terminate is false
		}
		Colony::clearDrones(); // delete drones

		if (terminate) { // if terminate is true
			std::cout << "population extinct\n"; // output this
			break; // stop
		}
		++t; // iterate simulation
        
        std::cout << t << std::endl; // output the time step
        
		if (t % dataInterval == 0) // if model run is finished
		{
			dataAnalysis(t); // start function dataAnalysis
		}
	}
    
    for (int i = 0; i < nColony; ++i) // loop through colonies
		delete population[i]; // delete colonies
	return 0; // return 0 as exit code
}
