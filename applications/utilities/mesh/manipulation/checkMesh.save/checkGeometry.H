#include "label.H"

namespace Foam
{
    class polyMesh;

    label checkGeometry
    (
        const polyMesh& mesh,
        bool checkPointNearness,
        bool checkCellDeterminant
    );
}
