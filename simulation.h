/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by Jan J. Kreider, Ido Pen, Boris H. Kramer and G. Sander van Doorn
 */

#pragma once

#include <fstream>
#include "colony.h"

void saveLifeHistory(const int generation, const Colony * const colonies[], const int sz, std::ofstream& file); // a function that saves the life history
