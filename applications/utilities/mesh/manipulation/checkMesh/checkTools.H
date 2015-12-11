#include "scalar.H"

namespace Foam
{
    class polyMesh;
    class surfaceWriter;
    class faceSet;
    class cellSet;

    void printMeshStats(const polyMesh& mesh, const bool allTopology);

    //- Write vtk representation of (assembled) faceSet to vtk file in
    //  postProcessing/ directory
    void mergeAndWrite(const surfaceWriter&, const faceSet&);

    //- Write vtk representation of (assembled) cellSet to vtk file in
    //  postProcessing/ directory
    void mergeAndWrite(const surfaceWriter&, const cellSet&);
}
