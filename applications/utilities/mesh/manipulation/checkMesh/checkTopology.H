#include "labelList.H"
#include "autoPtr.H"

namespace Foam
{
    // Forward Declarations
    class polyMesh;
    class pointSet;
    class coordSetWriter;
    class surfaceWriter;

    template<class PatchType>
    void checkPatch
    (
        const bool allGeometry,
        const word& name,
        const polyMesh& mesh,
        const PatchType& pp,
        const labelList& meshFaces,
        const labelList& meshEdges,
        pointSet& points
    );

    label checkTopology
    (
        const polyMesh& mesh,
        const bool allTopology,
        const bool allGeometry,
        autoPtr<surfaceWriter>& surfWriter,
        autoPtr<coordSetWriter>& setWriter
    );
}
