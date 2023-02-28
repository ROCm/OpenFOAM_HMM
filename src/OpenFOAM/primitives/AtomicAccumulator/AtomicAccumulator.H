#ifndef AtomicAccumulator_H
#define AtomicAccumulator_H

#include "VectorSpace.H"
#include "VectorSpaceOps.H"
#include "macros.H"
#include "ops.H"

#include <type_traits>

namespace Foam {

template <typename T>
class AtomicAccumulator {
private:
    T* val;
public:
    static_assert(std::is_scalar<T>::value, "T must be a scalar type");

    AtomicAccumulator(T &ref) : val(&ref) {}
    inline void operator+=(const T& rhs)
    {
        OMP(atomic update)
        *val += rhs;
    }
    inline void operator-=(const T& rhs)
    {
        OMP(atomic update)
        *val -= rhs;
    }
};

template <typename T>
inline typename std::enable_if<std::is_scalar<T>::value, AtomicAccumulator<T>>::type atomicAccumulator(T& ref)
{
    return AtomicAccumulator<T>(ref);
}

template <class Form, class Cmpt, direction Ncmpts>
class AtomicAccumulator<VectorSpace<Form, Cmpt, Ncmpts>> {
private:
    typedef VectorSpace<Form, Cmpt, Ncmpts> vector;

    vector* val;
public:
    AtomicAccumulator(vector &ref) : val(&ref) {}
    inline void operator+=(const vector& rhs)
    {
        VectorSpaceOps<Ncmpts,0>::eqOp(*val, rhs, atomicPlusEqOp<Cmpt>());
    }
    inline void operator-=(const vector& rhs)
    {
        VectorSpaceOps<Ncmpts,0>::eqOp(*val, rhs, atomicMinusEqOp<Cmpt>());
    }
};

template <class Form, class Cmpt, direction Ncmpts>
inline AtomicAccumulator<VectorSpace<Form, Cmpt, Ncmpts>> atomicAccumulator(VectorSpace<Form, Cmpt, Ncmpts>& ref)
{
    return AtomicAccumulator<VectorSpace<Form, Cmpt, Ncmpts>>(ref);
}

}

#endif // AtomicAccumulator_H
