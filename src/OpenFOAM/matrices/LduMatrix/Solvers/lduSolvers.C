#include "TPCG.H"
#include "TPBiCG.H"
#include "SmoothSolver.H"
#include "fieldTypes.H"

#define makeLduSolvers(Type, DType, LUType)                                   \
                                                                              \
    makeLduSolver(DiagonalSolver, Type, DType, LUType);                       \
    makeLduSymSolver(DiagonalSolver, Type, DType, LUType);                    \
    makeLduAsymSolver(DiagonalSolver, Type, DType, LUType);                   \
                                                                              \
    makeLduSolver(TPCG, Type, DType, LUType);                                 \
    makeLduSymSolver(TPCG, Type, DType, LUType);                              \
                                                                              \
    makeLduSolver(TPBiCG, Type, DType, LUType);                               \
    makeLduAsymSolver(TPBiCG, Type, DType, LUType);                           \
                                                                              \
    makeLduSolver(SmoothSolver, Type, DType, LUType);                         \
    makeLduSymSolver(SmoothSolver, Type, DType, LUType);                      \
    makeLduAsymSolver(SmoothSolver, Type, DType, LUType);

namespace Foam
{
    makeLduSolvers(scalar, scalar, scalar);
    makeLduSolvers(vector, scalar, scalar);
    makeLduSolvers(sphericalTensor, scalar, scalar);
    makeLduSolvers(symmTensor, scalar, scalar);
    makeLduSolvers(tensor, scalar, scalar);
};
