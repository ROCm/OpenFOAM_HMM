#include "LduMatrix.H"
#include "fieldTypes.H"

namespace Foam
{
    makeLduMatrix(scalar, scalar, scalar);
    makeLduMatrix(vector, scalar, scalar);
    makeLduMatrix(sphericalTensor, scalar, scalar);
    makeLduMatrix(symmTensor, scalar, scalar);
    makeLduMatrix(tensor, scalar, scalar);
};
