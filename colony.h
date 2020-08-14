/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by Jan J. Kreider, Ido Pen, Boris H. Kramer and G. Sander van Doorn
 */

#ifndef colony_h
#define colony_h

#include "individual.h"
#include <vector>

class Colony {
public:
	Colony(Queen * const); // constructor for Colony
	~Colony(); // destructor for Colony
    Queen const * getQueen() const { return queen; } // return a queen which is a pointer to a class Queen object
    std::vector<Female*> getWorkers() const { return workers; } // same for workers but return a vector of pointers to class Female objects
    void reproduce();
	void survival();
	void addWorkers();
	size_t size() const { return workers.size() + (queen != nullptr); } // function size() gives number of objects in vector workers + queen if there is a queen. So this is colony size in terms of number of individuals of all castes
	size_t nrOffspring() const { return offspring.size(); } // function nrOffspring() gives number of objects in vector offspring. So this is number of offspring
	bool isQueenless() const { return queen == nullptr; } // function isQueenless() is true if there is no queen and false if there is one
	Female *getGyne(); // a pointer to a Female object by getGyne
	void getData(size_t&, size_t&, size_t&, size_t&, std::vector<double>&, std::vector<double>&, std::vector<double>&, std::vector<double>&) const; // function to getData from Colonies
    static Male getDrone(); // get a copy of a Male object by getDrone
    static void clearDrones() { drones.clear(); } // delete data in vector drones
private:
	double colonyFitness; // Colonies have fitness
    std::vector<Female*> workers; // a vector including pointers to objects of class Female
	Queen* queen; // queen is a pointer to a class Queen object
	std::vector<Female*> offspring; // a vector including pointers to objects of class Female
	static std::vector<Male> drones; // a vector of class Male objects
};

#endif 
