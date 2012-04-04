/***************************************************************************\
|* Function Parser for C++ v4.0.3                                          *|
|*-------------------------------------------------------------------------*|
|* Copyright: Juha Nieminen, Joel Yliluoma                                 *|
\***************************************************************************/

// NOTE:
// This file contains only internal types for the function parser library.
// You don't need to include this file in your code. Include "fparser.hh"
// only.

#ifndef ONCE_FPARSER_AUX_H_
#define ONCE_FPARSER_AUX_H_

#include <cmath>
#include <cstring>

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
#include "mpfr/MpfrFloat.h"
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
#include "mpfr/GmpInt.h"
#endif

#ifdef ONCE_FPARSER_H_
namespace FUNCTIONPARSERTYPES
{
//==========================================================================
// Math funcs
//==========================================================================
    /* fp_pow() is a wrapper for std::pow()
     * that produces an identical value for
     * exp(1) ^ 2.0  (0x4000000000000000)
     * as exp(2.0)   (0x4000000000000000)
     * - std::pow() on x86_64
     * produces 2.0  (0x3FFFFFFFFFFFFFFF) instead!
     * See comments below for other special traits.
     */
    template<typename ValueT>
    inline ValueT fp_pow_with_exp_log(const ValueT& x, const ValueT& y)
    {
        // Requirements: x > 0.
        return fp_exp(fp_log(x) * y);
    }

    template<typename ValueT>
    inline ValueT fp_powi(ValueT x, unsigned long y)
    {
        // Requirements: y is non-negative integer.
        ValueT result(1);
        while(y != 0)
        {
            if(y & 1) { result *= x; y -= 1; }
            else      { x *= x;      y /= 2; }
        }
        return result;
    }

    template<typename ValueT>
    ValueT fp_pow(const ValueT& x, const ValueT& y)
    {
        if(x == ValueT(1)) return ValueT(1);
        // y is now zero or positive

      if(IsIntegerConst(y))
        {
            // Use fast binary exponentiation algorithm
            // See http://en.wikipedia.org/wiki/Exponentiation_by_squaring
            if(y >= ValueT(0))
                return fp_powi(x,              (long)y);
            else
                return ValueT(1) / fp_powi(x, -(long)y);
        }

      if(y >= ValueT(0))
        {
            // y is now positive. Calculate using exp(log(x)*y).
            // See http://en.wikipedia.org/wiki/Exponentiation#Real_powers
            if(x > ValueT(0)) return fp_pow_with_exp_log(x, y);
            if(x == ValueT(0)) return ValueT(0);
            // At this point, y > 0.0 and x is known to be < 0.0,
            // because positive and zero cases are already handled.
            //
            if(!IsIntegerConst(y*ValueT(16)))
                return -fp_pow_with_exp_log(-x, y);
            // ^This is not technically correct, but it allows
            // functions such as cbrt(x^5), that is, x^(5/3),
            // to be evaluated when x is negative.
            // It is too complicated (and slow) to test whether y
            // is a formed from a ratio of an integer to an odd integer.
            // (And due to floating point inaccuracy, pointless too.)
            // For example, x^1.30769230769... is
            // actually x^(17/13), i.e. (x^17) ^ (1/13).
            // (-5)^(17/13) gives us now -8.204227562330453.
            // To see whether the result is right, we can test the given
            // root: (-8.204227562330453)^13 gives us the value of (-5)^17,
            // which proves that the expression was correct.
            //
            // The y*16 check prevents e.g. (-4)^(3/2) from being calculated,
            // as it would confuse functioninfo when pow() returns no error
            // but sqrt() does when the formula is converted into sqrt(x)*x.
            //
            // The errors in this approach are:
            //     (-2)^sqrt(2) should produce NaN
            //                  or actually sqrt(2)I + 2^sqrt(2),
            //                  produces -(2^sqrt(2)) instead.
            //                  (Impact: Neglible)
            // Thus, at worst, we're changing a NaN (or complex)
            // result into a negative real number result.
        }
        else
        {
            // y is negative. Utilize the x^y = 1/(x^-y) identity.
            if(x > ValueT(0)) return fp_pow_with_exp_log(ValueT(1) / x, -y);
            if(x < ValueT(0))
            {
                if(!IsIntegerConst(y*ValueT(-16)))
                    return -fp_pow_with_exp_log(ValueT(-1) / x, -y);
                // ^ See comment above.
            }
            // Remaining case: 0.0 ^ negative number
        }
        // This is reached when:
        //      x=0 and y<0
        //      x<0 and y*16 is either positive or negative integer
        // It is used for producing error values and as a safe fallback.
        return fp_pow_base(x, y);
    }

// -------------------------------------------------------------------------
// double
// -------------------------------------------------------------------------
    inline double fp_abs(double x) { return fabs(x); }
    inline double fp_acos(double x) { return acos(x); }
    inline double fp_asin(double x) { return asin(x); }
    inline double fp_atan(double x) { return atan(x); }
    inline double fp_atan2(double x, double y) { return atan2(x, y); }
#ifdef FP_SUPPORT_CBRT
    inline double fp_cbrt(double x) { return cbrt(x); }
#else
    inline double fp_cbrt(double x) { return x>0 ?  exp(log( x)/3.0)
                                           : x<0 ? -exp(log(-x)/3.0)
                                           : 0.0; }
#endif
    inline double fp_ceil(double x) { return ceil(x); }
    inline double fp_cos(double x) { return cos(x); }
    inline double fp_cosh(double x) { return cosh(x); }
    inline double fp_exp(double x) { return exp(x); }
    inline double fp_floor(double x) { return floor(x); }
    inline double fp_int(double x) { return floor(x + .5); }
    inline double fp_log(double x) { return log(x); }
    inline double fp_log10(double x)
    { return log(x) *
            0.434294481903251827651128918916605082294397005803666566; }
    inline double fp_mod(double x, double y) { return fmod(x, y); }
    inline double fp_sin(double x) { return sin(x); }
    inline double fp_sinh(double x) { return sinh(x); }
    inline double fp_sqrt(double x) { return sqrt(x); }
    inline double fp_tan(double x) { return tan(x); }
    inline double fp_tanh(double x) { return tanh(x); }

#ifndef FP_SUPPORT_ASINH
    inline double fp_asinh(double x) { return log(x + sqrt(x*x + 1.0)); }
    inline double fp_acosh(double x) { return log(x + sqrt(x*x - 1.0)); }
    inline double fp_atanh(double x) { return log((1.0+x) / (1.0-x)) * 0.5; }
#else
    inline double fp_asinh(double x) { return asinh(x); }
    inline double fp_acosh(double x) { return acosh(x); }
    inline double fp_atanh(double x) { return atanh(x); }
#endif // FP_SUPPORT_ASINH

    inline double fp_trunc(double x) { return x<0.0 ? ceil(x) : floor(x); }

    inline double fp_pow_base(double x, double y) { return pow(x, y); }

#ifndef FP_SUPPORT_LOG2
    inline double fp_log2(double x)
    { return log(x) * 1.4426950408889634073599246810018921374266459541529859; }
#else
    inline double fp_log2(double x) { return log2(x); }
#endif // FP_SUPPORT_LOG2

    inline double fp_exp2(double x) { return fp_pow(2.0, x); }

#ifdef FP_EPSILON
    template<typename Value_t>
    inline Value_t fp_epsilon() { return FP_EPSILON; }
#else
    template<typename Value_t>
    inline Value_t fp_epsilon() { return 0.0; }
#endif

#ifdef FP_EPSILON
    inline bool FloatEqual(double a, double b)
    { return fabs(a - b) <= fp_epsilon<double>(); }
#else
    inline bool FloatEqual(double a, double b)
    { return a == b; }
#endif // FP_EPSILON

    inline bool IsIntegerConst(double a)
    { return FloatEqual(a, (double)(long)a); }


// -------------------------------------------------------------------------
// float
// -------------------------------------------------------------------------
#ifdef FP_SUPPORT_FLOAT_TYPE
    inline float fp_abs(float x) { return fabsf(x); }
    inline float fp_acos(float x) { return acosf(x); }
    inline float fp_asin(float x) { return asinf(x); }
    inline float fp_atan(float x) { return atanf(x); }
    inline float fp_atan2(float x, float y) { return atan2f(x, y); }
#ifdef FP_SUPPORT_CBRT
    inline float fp_cbrt(float x) { return cbrtf(x); }
#else
    inline float fp_cbrt(float x) { return x>0 ?  expf(logf( x)/3.0f)
                                         : x<0 ? -expf(logf(-x)/3.0f)
                                         : 0.0f; }
#endif
    inline float fp_ceil(float x) { return ceilf(x); }
    inline float fp_cos(float x) { return cosf(x); }
    inline float fp_cosh(float x) { return coshf(x); }
    inline float fp_exp(float x) { return expf(x); }
    inline float fp_floor(float x) { return floorf(x); }
    inline float fp_int(float x) { return floorf(x + .5F); }
    inline float fp_log(float x) { return logf(x); }
    inline float fp_log10(float x)
    { return logf(x) *
            0.434294481903251827651128918916605082294397005803666566F; }
    inline float fp_mod(float x, float y) { return fmodf(x, y); }
    inline float fp_sin(float x) { return sinf(x); }
    inline float fp_sinh(float x) { return sinhf(x); }
    inline float fp_sqrt(float x) { return sqrtf(x); }
    inline float fp_tan(float x) { return tanf(x); }
    inline float fp_tanh(float x) { return tanhf(x); }

#ifndef FP_SUPPORT_ASINH
    inline float fp_asinh(float x) { return logf(x + sqrt(x*x + 1.0F)); }
    inline float fp_acosh(float x) { return logf(x + sqrt(x*x - 1.0F)); }
    inline float fp_atanh(float x) { return logf((1.0F+x) / (1.0F-x)) * 0.5F; }
#else
    inline float fp_asinh(float x) { return asinhf(x); }
    inline float fp_acosh(float x) { return acoshf(x); }
    inline float fp_atanh(float x) { return atanhf(x); }
#endif // FP_SUPPORT_ASINH

    inline float fp_trunc(float x) { return x<0.0F ? ceilf(x) : floorf(x); }

    inline float fp_pow_base(float x, float y) { return powf(x, y); }

#ifndef FP_SUPPORT_LOG2
    inline float fp_log2(float x)
    { return logf(x) *
            1.4426950408889634073599246810018921374266459541529859F; }
#else
    inline float fp_log2(float x) { return log2f(x); }
#endif // FP_SUPPORT_LOG2

    inline float fp_exp2(float x) { return fp_pow(2.0F, x); }

#ifdef FP_EPSILON
    template<>
    inline float fp_epsilon<float>() { return 1e-6F; }
#else
    template<>
    inline float fp_epsilon<float>() { return 0.0F; }
#endif

#ifdef FP_EPSILON
    inline bool FloatEqual(float a, float b)
    { return fabsf(a - b) <= fp_epsilon<float>(); }
#else
    inline bool FloatEqual(float a, float b)
    { return a == b; }
#endif // FP_EPSILON

    inline bool IsIntegerConst(float a)
    { return FloatEqual(a, (float)(long)a); }
#endif // FP_SUPPORT_FLOAT_TYPE



// -------------------------------------------------------------------------
// long double
// -------------------------------------------------------------------------
#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
    inline long double fp_abs(long double x) { return fabsl(x); }
    inline long double fp_acos(long double x) { return acosl(x); }
    inline long double fp_asin(long double x) { return asinl(x); }
    inline long double fp_atan(long double x) { return atanl(x); }
    inline long double fp_atan2(long double x, long double y)
    { return atan2l(x, y); }
#ifdef FP_SUPPORT_CBRT
    inline long double fp_cbrt(long double x) { return cbrtl(x); }
#else
    inline long double fp_cbrt(long double x)
        { return x>0 ?  expl(logl( x)/3.0l)
               : x<0 ? -expl(logl(-x)/3.0l)
               : 0.0l; }
#endif
    inline long double fp_ceil(long double x) { return ceill(x); }
    inline long double fp_cos(long double x) { return cosl(x); }
    inline long double fp_cosh(long double x) { return coshl(x); }
    inline long double fp_exp(long double x) { return expl(x); }
    inline long double fp_floor(long double x) { return floorl(x); }
    inline long double fp_int(long double x) { return floorl(x + .5L); }
    inline long double fp_log(long double x) { return logl(x); }
    inline long double fp_log10(long double x)
    { return logl(x) *
            0.434294481903251827651128918916605082294397005803666566L; }
    inline long double fp_mod(long double x, long double y)
    { return fmodl(x, y); }
    inline long double fp_sin(long double x) { return sinl(x); }
    inline long double fp_sinh(long double x) { return sinhl(x); }
    inline long double fp_sqrt(long double x) { return sqrtl(x); }
    inline long double fp_tan(long double x) { return tanl(x); }
    inline long double fp_tanh(long double x) { return tanhl(x); }

#ifndef FP_SUPPORT_ASINH
    inline long double fp_asinh(long double x)
    { return logl(x + sqrtl(x*x + 1.0L)); }
    inline long double fp_acosh(long double x)
    { return logl(x + sqrtl(x*x - 1.0L)); }
    inline long double fp_atanh(long double x)
    { return logl((1.0L+x) / (1.0L-x)) * 0.5L; }
#else
    inline long double fp_asinh(long double x) { return asinhl(x); }
    inline long double fp_acosh(long double x) { return acoshl(x); }
    inline long double fp_atanh(long double x) { return atanhl(x); }
#endif // FP_SUPPORT_ASINH

    inline long double fp_trunc(long double x)
    { return x<0.0L ? ceill(x) : floorl(x); }

    inline long double fp_pow_base(long double x, long double y)
    { return powl(x, y); }

#ifndef FP_SUPPORT_LOG2
    inline long double fp_log2(long double x)
    { return log(x) * 1.4426950408889634073599246810018921374266459541529859L; }
#else
    inline long double fp_log2(long double x) { return log2l(x); }
#endif // FP_SUPPORT_LOG2

    inline long double fp_exp2(long double x) { return fp_pow(2.0L, x); }

#ifdef FP_EPSILON
    inline bool FloatEqual(long double a, long double b)
    { return fabsl(a - b) <= fp_epsilon<double>(); }
#else
    inline bool FloatEqual(long double a, long double b)
    { return a == b; }
#endif // FP_EPSILON

    inline bool IsIntegerConst(long double a)
    { return FloatEqual(a, (long double)(long)a); }
#endif // FP_SUPPORT_LONG_DOUBLE_TYPE


// -------------------------------------------------------------------------
// Long int
// -------------------------------------------------------------------------
    inline long fp_abs(long x) { return x < 0 ? -x : x; }
    inline long fp_acos(long) { return 0; }
    inline long fp_asin(long) { return 0; }
    inline long fp_atan(long) { return 0; }
    inline long fp_atan2(long, long) { return 0; }
    inline long fp_cbrt(long) { return 0; }
    inline long fp_ceil(long x) { return x; }
    inline long fp_cos(long) { return 0; }
    inline long fp_cosh(long) { return 0; }
    inline long fp_exp(long) { return 0; }
    inline long fp_floor(long x) { return x; }
    inline long fp_int(long x) { return x; }
    inline long fp_log(long) { return 0; }
    inline long fp_log10(long) { return 0; }
    inline long fp_mod(long x, long y) { return x % y; }
    inline long fp_pow(long, long) { return 0; }
    inline long fp_sin(long) { return 0; }
    inline long fp_sinh(long) { return 0; }
    inline long fp_sqrt(long) { return 1; }
    inline long fp_tan(long) { return 0; }
    inline long fp_tanh(long) { return 0; }
    inline long fp_asinh(long) { return 0; }
    inline long fp_acosh(long) { return 0; }
    inline long fp_atanh(long) { return 0; }
    inline long fp_trunc(long x) { return x; }
    inline long fp_pow_base(long, long) { return 0; }
    inline long fp_log2(long) { return 0; }
    inline long fp_exp2(long) { return 0; }

    template<>
    inline long fp_epsilon<long>() { return 0; }

    inline bool IsIntegerConst(long) { return true; }


// -------------------------------------------------------------------------
// MpfrFloat
// -------------------------------------------------------------------------
#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
    inline MpfrFloat fp_abs(const MpfrFloat& x) { return MpfrFloat::abs(x); }
    inline MpfrFloat fp_acos(const MpfrFloat& x) { return MpfrFloat::acos(x); }
    inline MpfrFloat fp_asin(const MpfrFloat& x) { return MpfrFloat::asin(x); }
    inline MpfrFloat fp_atan(const MpfrFloat& x) { return MpfrFloat::atan(x); }
    inline MpfrFloat fp_atan2(const MpfrFloat& x, const MpfrFloat& y) { return MpfrFloat::atan2(x, y); }
    inline MpfrFloat fp_cbrt(const MpfrFloat& x) { return MpfrFloat::cbrt(x); }
    inline MpfrFloat fp_ceil(const MpfrFloat& x) { return MpfrFloat::ceil(x); }
    inline MpfrFloat fp_cos(const MpfrFloat& x) { return MpfrFloat::cos(x); }
    inline MpfrFloat fp_cosh(const MpfrFloat& x) { return MpfrFloat::cosh(x); }
    inline MpfrFloat fp_exp(const MpfrFloat& x) { return MpfrFloat::exp(x); }
    inline MpfrFloat fp_floor(const MpfrFloat& x) { return MpfrFloat::floor(x); }
    inline MpfrFloat fp_int(const MpfrFloat& x) { return MpfrFloat::round(x); }
    inline MpfrFloat fp_log(const MpfrFloat& x) { return MpfrFloat::log(x); }
    inline MpfrFloat fp_log10(const MpfrFloat& x) { return MpfrFloat::log10(x); }
    inline MpfrFloat fp_mod(const MpfrFloat& x, const MpfrFloat& y) { return x % y; }
    inline MpfrFloat fp_sin(const MpfrFloat& x) { return MpfrFloat::sin(x); }
    inline MpfrFloat fp_sinh(const MpfrFloat& x) { return MpfrFloat::sinh(x); }
    inline MpfrFloat fp_sqrt(const MpfrFloat& x) { return MpfrFloat::sqrt(x); }
    inline MpfrFloat fp_tan(const MpfrFloat& x) { return MpfrFloat::tan(x); }
    inline MpfrFloat fp_tanh(const MpfrFloat& x) { return MpfrFloat::tanh(x); }
    inline MpfrFloat fp_asinh(const MpfrFloat& x) { return MpfrFloat::asinh(x); }
    inline MpfrFloat fp_acosh(const MpfrFloat& x) { return MpfrFloat::acosh(x); }
    inline MpfrFloat fp_atanh(const MpfrFloat& x) { return MpfrFloat::atanh(x); }
    inline MpfrFloat fp_trunc(const MpfrFloat& x) { return MpfrFloat::trunc(x); }

    inline MpfrFloat fp_pow(const MpfrFloat& x, const MpfrFloat& y) { return MpfrFloat::pow(x, y); }
    inline MpfrFloat fp_pow_base(const MpfrFloat& x, const MpfrFloat& y) { return MpfrFloat::pow(x, y); }

    inline MpfrFloat fp_log2(const MpfrFloat& x) { return MpfrFloat::log2(x); }
    inline MpfrFloat fp_exp2(const MpfrFloat& x) { return MpfrFloat::exp2(x); }

    inline bool IsIntegerConst(const MpfrFloat& a) { return a.isInteger(); }

    template<>
    inline MpfrFloat fp_epsilon<MpfrFloat>() { return MpfrFloat::someEpsilon(); }
#endif // FP_SUPPORT_MPFR_FLOAT_TYPE


// -------------------------------------------------------------------------
// GMP int
// -------------------------------------------------------------------------
#ifdef FP_SUPPORT_GMP_INT_TYPE
    inline GmpInt fp_abs(GmpInt x) { return GmpInt::abs(x); }
    inline GmpInt fp_acos(GmpInt) { return 0; }
    inline GmpInt fp_asin(GmpInt) { return 0; }
    inline GmpInt fp_atan(GmpInt) { return 0; }
    inline GmpInt fp_atan2(GmpInt, GmpInt) { return 0; }
    inline GmpInt fp_cbrt(GmpInt) { return 0; }
    inline GmpInt fp_ceil(GmpInt x) { return x; }
    inline GmpInt fp_cos(GmpInt) { return 0; }
    inline GmpInt fp_cosh(GmpInt) { return 0; }
    inline GmpInt fp_exp(GmpInt) { return 0; }
    inline GmpInt fp_floor(GmpInt x) { return x; }
    inline GmpInt fp_int(GmpInt x) { return x; }
    inline GmpInt fp_log(GmpInt) { return 0; }
    inline GmpInt fp_log10(GmpInt) { return 0; }
    inline GmpInt fp_mod(GmpInt x, GmpInt y) { return x % y; }
    inline GmpInt fp_pow(GmpInt, GmpInt) { return 0; }
    inline GmpInt fp_sin(GmpInt) { return 0; }
    inline GmpInt fp_sinh(GmpInt) { return 0; }
    inline GmpInt fp_sqrt(GmpInt) { return 0; }
    inline GmpInt fp_tan(GmpInt) { return 0; }
    inline GmpInt fp_tanh(GmpInt) { return 0; }
    inline GmpInt fp_asinh(GmpInt) { return 0; }
    inline GmpInt fp_acosh(GmpInt) { return 0; }
    inline GmpInt fp_atanh(GmpInt) { return 0; }
    inline GmpInt fp_trunc(GmpInt x) { return x; }
    inline GmpInt fp_pow_base(GmpInt, GmpInt) { return 0; }
    inline GmpInt fp_log2(GmpInt) { return 0; }
    inline GmpInt fp_exp2(GmpInt) { return 0; }

    template<>
    inline GmpInt fp_epsilon<GmpInt>() { return 0; }

    inline bool IsIntegerConst(GmpInt) { return true; }
#endif // FP_SUPPORT_GMP_INT_TYPE


  
// -------------------------------------------------------------------------
// Comparison
// -------------------------------------------------------------------------
#ifdef FP_EPSILON
    template<typename Value_t>
    inline bool fp_equal(Value_t x, Value_t y)
    { return fp_abs(x - y) <= fp_epsilon<Value_t>(); }

    template<typename Value_t>
    inline bool fp_nequal(Value_t x, Value_t y)
    { return fp_abs(x - y) > fp_epsilon<Value_t>(); }

    template<typename Value_t>
    inline bool fp_less(Value_t x, Value_t y)
    { return x < y - fp_epsilon<Value_t>(); }

    template<typename Value_t>
    inline bool fp_lessOrEq(Value_t x, Value_t y)
    { return x <= y + fp_epsilon<Value_t>(); }
#else // FP_EPSILON
    template<typename Value_t>
    inline bool fp_equal(Value_t x, Value_t y) { return x == y; }

    template<typename Value_t>
    inline bool fp_nequal(Value_t x, Value_t y) { return x != y; }

    template<typename Value_t>
    inline bool fp_less(Value_t x, Value_t y) { return x < y; }

    template<typename Value_t>
    inline bool fp_lessOrEq(Value_t x, Value_t y) { return x <= y; }
#endif // FP_EPSILON

    template<>
    inline bool fp_equal<long>(long x, long y) { return x == y; }

    template<>
    inline bool fp_nequal<long>(long x, long y) { return x != y; }

    template<>
    inline bool fp_less<long>(long x, long y) { return x < y; }

    template<>
    inline bool fp_lessOrEq<long>(long x, long y) { return x <= y; }

#ifdef FP_SUPPORT_GMP_INT_TYPE
    template<>
    inline bool fp_equal<GmpInt>(GmpInt x, GmpInt y) { return x == y; }

    template<>
    inline bool fp_nequal<GmpInt>(GmpInt x, GmpInt y) { return x != y; }

    template<>
    inline bool fp_less<GmpInt>(GmpInt x, GmpInt y) { return x < y; }

    template<>
    inline bool fp_lessOrEq<GmpInt>(GmpInt x, GmpInt y) { return x <= y; }
#endif
} // namespace FUNCTIONPARSERTYPES

#endif // ONCE_FPARSER_H_


#ifndef FP_DISABLE_DOUBLE_TYPE
#define FUNCTIONPARSER_INSTANTIATE_DOUBLE \
    template class FunctionParserBase<double>;
#else
#define FUNCTIONPARSER_INSTANTIATE_DOUBLE
#endif

#ifdef FP_SUPPORT_FLOAT_TYPE
#define FUNCTIONPARSER_INSTANTIATE_FLOAT \
    template class FunctionParserBase<float>;
#else
#define FUNCTIONPARSER_INSTANTIATE_FLOAT
#endif

#ifdef FP_SUPPORT_LONG_DOUBLE_TYPE
#define FUNCTIONPARSER_INSTANTIATE_LONG_DOUBLE \
    template class FunctionParserBase<long double>;
#else
#define FUNCTIONPARSER_INSTANTIATE_LONG_DOUBLE
#endif

#ifdef FP_SUPPORT_LONG_INT_TYPE
#define FUNCTIONPARSER_INSTANTIATE_LONG_INT \
    template class FunctionParserBase<long>;
#else
#define FUNCTIONPARSER_INSTANTIATE_LONG_INT
#endif

#ifdef FP_SUPPORT_MPFR_FLOAT_TYPE
#define FUNCTIONPARSER_INSTANTIATE_MPFR_FLOAT \
    template class FunctionParserBase<MpfrFloat>;
#else
#define FUNCTIONPARSER_INSTANTIATE_MPFR_FLOAT
#endif

#ifdef FP_SUPPORT_GMP_INT_TYPE
#define FUNCTIONPARSER_INSTANTIATE_GMP_INT \
    template class FunctionParserBase<GmpInt>;
#else
#define FUNCTIONPARSER_INSTANTIATE_GMP_INT
#endif

/* Add 'FUNCTIONPARSER_INSTANTIATE_TYPES' at the end of all .cc files
   containing FunctionParserBase implementations.
 */
#define FUNCTIONPARSER_INSTANTIATE_TYPES \
    FUNCTIONPARSER_INSTANTIATE_DOUBLE \
    FUNCTIONPARSER_INSTANTIATE_FLOAT \
    FUNCTIONPARSER_INSTANTIATE_LONG_DOUBLE \
    FUNCTIONPARSER_INSTANTIATE_LONG_INT \
    FUNCTIONPARSER_INSTANTIATE_MPFR_FLOAT \
    FUNCTIONPARSER_INSTANTIATE_GMP_INT

#endif // ONCE_FPARSER_AUX_H_
