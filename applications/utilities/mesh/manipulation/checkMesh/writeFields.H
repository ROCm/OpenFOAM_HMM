#include "fvMesh.H"

namespace Foam
{
    void writeFields
    (
        const fvMesh&,
        const HashSet<word>& selectedFields
    );
}
