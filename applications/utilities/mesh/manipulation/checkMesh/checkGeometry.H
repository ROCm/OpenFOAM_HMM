#include "label.H"
#include "HashSet.H"
#include "labelVector.H"

namespace Foam
{
    // Forward Declarations
    class polyMesh;
    class wedgePolyPatch;
    template<class T> class writer;
    class surfaceWriter;

    label findOppositeWedge(const polyMesh&, const wedgePolyPatch&);

    //- Check wedge orientation
    bool checkWedges
    (
        const polyMesh&,
        const bool report,
        const Vector<label>&,
        labelHashSet*
    );

    //- Check 0th vertex on coupled faces
    bool checkCoupledPoints(const polyMesh&, const bool report, labelHashSet*);

    label checkGeometry
    (
        const polyMesh& mesh,
        const bool allGeometry,
        autoPtr<surfaceWriter>& surfWriter,
        const autoPtr<writer<scalar>>& setWriter
    );
}
