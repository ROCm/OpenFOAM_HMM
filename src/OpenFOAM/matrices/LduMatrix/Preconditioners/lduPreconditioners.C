#include "NoPreconditioner.H"
#include "DiagonalPreconditioner.H"
#include "TDILUPreconditioner.H"
#include "fieldTypes.H"

#define makeLduPreconditioners(Type, DType, LUType)                           \
                                                                              \
    makeLduPreconditioner(NoPreconditioner, Type, DType, LUType);             \
    makeLduSymPreconditioner(NoPreconditioner, Type, DType, LUType);          \
    makeLduAsymPreconditioner(NoPreconditioner, Type, DType, LUType);         \
                                                                              \
    makeLduPreconditioner(DiagonalPreconditioner, Type, DType, LUType);       \
    makeLduSymPreconditioner(DiagonalPreconditioner, Type, DType, LUType);    \
    makeLduAsymPreconditioner(DiagonalPreconditioner, Type, DType, LUType);   \
                                                                              \
    makeLduPreconditioner(TDILUPreconditioner, Type, DType, LUType);          \
    makeLduAsymPreconditioner(TDILUPreconditioner, Type, DType, LUType);

namespace Foam
{
    makeLduPreconditioners(scalar, scalar, scalar);
    makeLduPreconditioners(vector, scalar, scalar);
    makeLduPreconditioners(sphericalTensor, scalar, scalar);
    makeLduPreconditioners(symmTensor, scalar, scalar);
    makeLduPreconditioners(tensor, scalar, scalar);
};
