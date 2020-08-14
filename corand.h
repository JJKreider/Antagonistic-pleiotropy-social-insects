/*
 Publication: Antagonistic pleiotropy and the evolution of eusociality
 Authors: Jan J. Kreider, Ido Pen, Boris H. Kramer
 code written by G. Sander van Doorn
 */

#ifndef corand_h
#define corand_h

#include "random.h"
#include "linalg.h"

namespace rnd
{
    class MultiNormal
    {
    public:
        MultiNormal(const Matrix&, const Vector&);
        MultiNormal(const Matrix&);
        Vector operator()() const;
        friend std::ostream& operator<< (std::ostream&, const MultiNormal&);
    private:
        const Vector::size_type n;
        const Vector mu;
        Vector sigma;
        Matrix U;
    };
}

#endif

