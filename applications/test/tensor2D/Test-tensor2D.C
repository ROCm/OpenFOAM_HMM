#include "tensor2D.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    vector2D v1(1, 2), v2(3, 4);
    tensor2D t3 = v1*v2;

    Info<< v1 << "*" << v2 << " = " << t3 << endl;

    {
        Info<< "rows:" << nl;
        for (direction i=0; i < 2; ++i)
        {
            Info<< "  (" << i << ") = " << t3.row(i) << nl;
        }
    }

    {
        Info<< "cols:" << nl;
        for (direction i=0; i < 2; ++i)
        {
            Info<< "  (" << i << ") = " << t3.col(i) << nl;
        }
        Info<< "col<0> = " << t3.col<0>() << nl;
        Info<< "col<1> = " << t3.col<1>() << nl;
        // Compilation error:  Info << "col<3> = " << t3.col<3>() << nl;

        t3.col<0>({0, 2});
        Info<< "replaced col<0> = " << t3.col<0>() << nl;
        Info<< "tensor " << t3 << nl;

        t3.row<1>(Zero);
        Info<< "replaced row<1> = " << t3.row<1>() << nl;
        Info<< "tensor " << t3 << nl;
    }
    Info<< nl;

    return 0;
}
