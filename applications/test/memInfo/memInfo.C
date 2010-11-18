#include "memInfo.H"
#include "IOstreams.H"
#include "List.H"
#include "vector.H"

using namespace Foam;

int main()
{
    memInfo m;

    Info<< m << endl;

    List<vector> l(10000000, vector::one);

    Info<< m.update() << endl;

    return 0;
}
