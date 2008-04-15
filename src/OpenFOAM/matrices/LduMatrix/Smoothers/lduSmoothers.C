#include "TGaussSeidelSmoother.H"
#include "fieldTypes.H"

#define makeLduSmoothers(Type, DType, LUType)                                 \
                                                                              \
    makeLduSmoother(TGaussSeidelSmoother, Type, DType, LUType);               \
    makeLduSymSmoother(TGaussSeidelSmoother, Type, DType, LUType);            \
    makeLduAsymSmoother(TGaussSeidelSmoother, Type, DType, LUType);

namespace Foam
{
    makeLduSmoothers(scalar, scalar, scalar);
    makeLduSmoothers(vector, scalar, scalar);
    makeLduSmoothers(sphericalTensor, scalar, scalar);
    makeLduSmoothers(symmTensor, scalar, scalar);
    makeLduSmoothers(tensor, scalar, scalar);
};
