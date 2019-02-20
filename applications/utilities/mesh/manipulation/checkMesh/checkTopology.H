#include "label.H"
#include "autoPtr.H"
#include "writer.H"

namespace Foam
{
    class polyMesh;
    class pointSet;
    class surfaceWriter;

    template<class PatchType>
    void checkPatch
    (
        const bool allGeometry,
        const word& name,
        const PatchType& pp,
        pointSet& points
    );

    label checkTopology
    (
        const polyMesh& mesh,
        const bool allTopology,
        const bool allGeometry,
        autoPtr<surfaceWriter>& surfWriter,
        const autoPtr<writer<scalar>>& setWriter
    );
}
