/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by Jan J. Kreider, Ido Pen, Boris H. Kramer and G. Sander van Doorn
 */

#ifndef individual_h
#define individual_h

#include<array>
#include<vector>
#include "parameters.h"
#include "corand.h"

class Gene {
public:
	Gene() : dx(K, 0.0) {} // instances of class Gene are constructed with a vector of length K with all entries = 0.0
	Gene(const Gene&) = default;
    void mutate() { dx += geneticArchitecture->operator()(); } // function increases dx by value from multinormal distribution
    double operator()(const int &cs, const int &tr, const int &ag) { return dx[ag + (cs + tr * nCaste) * ageMax]; } // return dx value for age class and trait and caste from the vector dx
	static void initGeneticArchitecture(); // function to initialize genetic architecture
private:
	Vector dx;
	static rnd::MultiNormal const * geneticArchitecture; // geneticArchitecture points to a multinormal distribution
	const static int K = nCaste * nTrait * ageMax; // K determines the size of the partial correlation matrix
};

class Male;
class Individual {
	friend class Male;
public:
    enum Traits { survival, fecundity }; // Traits: 0 = survival, 1 = fecundity
	void incrementAge() { ++age; } // function to increase age
	void extrinsicMortality(double deathRate) { extrinsicSurvival = 1.0 - deathRate; } // function changes the extrinsic survival probability when given the death rate
    bool doesSurvive() const { return age < ageMax && rnd::bernoulli(extrinsicSurvival * lifeHistory(survival, age)); } // function doesSurvive delivers a bool for survival depending on max age and the survival probability
    friend std::ostream& operator<< (std::ostream&, const Individual&); // output stream
	double getLifeHistoryTrait(const Traits &tr) const { return lifeHistory(static_cast<Matrix::size_type>(tr), age); } // this function returns a age-specific trait value
	const Matrix& getLifeHistory() const { return lifeHistory; } // function delivers the lifeHistory matrix
	int getAge() const { return age; } // function to return age
	void getData(std::vector<double>&, std::vector<double>&) const;
protected:
    Individual() : lifeHistory(nTrait, ageMax), age(0) {} // individuals are constructed with a life history matrix and age 0
	Individual(const Individual&) = default; // constructor for inheritance
	std::array<Gene, nGene> genome; // genome array
	int age;
	double extrinsicSurvival;
    Matrix lifeHistory;
	virtual void develop() = 0;
    void mutate();
};

class Queen;
class Female : public Individual { // derived class
	friend class Queen;
public:
	Female() = default;
	Female(const Female&) = delete;
	Female(Queen const * const mother); // inheritance
	~Female() { ++ageOfDeath[age]; } // destructor. write age of death to array
	static std::array<size_t, ageMax + 1u> ageOfDeath; // death ages of Females
private:
	void develop();
};

class Male {
	friend class Female;
public:
	Male(Individual const * const); // interitance
	Male(const Male&) = default;
private:
	Male() = delete;
	void mutate();
	std::array<Gene, nGene / 2> genome; // genome with length nGene/2, i.e. haploid
};

class Queen: public Individual { // derived class
public:
	friend class Female;
	Queen() = delete;
	Queen(Female const * const, Male); // inheritance
    ~Queen() { ++ageOfDeath[age]; } // destructor. write age of death to array
	static std::array<size_t, ageMax + 1u> ageOfDeath; // death ages of Queens
private:
	Male sperm; // sperm is a copy of a Male
	void develop();
};

#endif
