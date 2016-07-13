/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fft.H"

#include <fftw3.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void fft::transform
(
    complexField& field,
    const UList<int>& nn,
    transformDirection dir
)
{
    forAll(nn, idim)
    {
        // Check for power of two
        unsigned int dimCount = nn[idim];
        if (!dimCount || (dimCount & (dimCount - 1)))
        {
            FatalErrorInFunction
                << "number of elements in direction " << idim
                << " is not a power of 2" << endl
                << "    Number of elements in each direction = " << nn
                << abort(FatalError);
        }
    }

    // Copy field into fftw containers
    label N = field.size();
    fftw_complex in[N], out[N];

    forAll(field, i)
    {
        in[i][0] = field[i].Re();
        in[i][1] = field[i].Im();
    }

    // Create the plan
    // FFTW_FORWARD = -1
    // FFTW_BACKWARD = 1

    // 1-D plan
    // fftw_plan plan =
    //    fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Generic 1..3-D plan
    label rank = nn.size();
    fftw_plan plan =
        fftw_plan_dft(rank, nn.begin(), in, out, dir, FFTW_ESTIMATE);

    // Compute the FFT
    fftw_execute(plan);

    forAll(field, i)
    {
        field[i].Re() = out[i][0];
        field[i].Im() = out[i][1];
    }

    fftw_destroy_plan(plan);

/*
    fftw_real in[N], out[N];
    // Create a plan for real data input
    // - generates 1-sided FFT
    // - direction given by FFTW_R2HC or FFTW_HC2R.
    // - in[N] is the array of real input val ues
    // - out[N] is the array of N/2 real valuea followed by N/2 complex values
    // - 0 component = DC component
    fftw_plan plan = fftw_plan_r2r_2d(N, in, out, FFTW_R2HC, FFTW_ESTIMATE)

    // Compute the FFT
    fftw_execute(plan);

    fftw_destroy_plan(plan);
*/
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

tmp<complexField> fft::forwardTransform
(
    const tmp<complexField>& tfield,
    const UList<int>& nn
)
{
    tmp<complexField> tfftField(new complexField(tfield));

    transform(tfftField.ref(), nn, FORWARD_TRANSFORM);

    tfield.clear();

    return tfftField;
}


tmp<complexField> fft::reverseTransform
(
    const tmp<complexField>& tfield,
    const UList<int>& nn
)
{
    tmp<complexField> tifftField(new complexField(tfield));

    transform(tifftField.ref(), nn, REVERSE_TRANSFORM);

    tfield.clear();

    return tifftField;
}


tmp<complexVectorField> fft::forwardTransform
(
    const tmp<complexVectorField>& tfield,
    const UList<int>& nn
)
{
    tmp<complexVectorField> tfftVectorField
    (
        new complexVectorField
        (
            tfield().size()
        )
    );

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        tfftVectorField.ref().replace
        (
            cmpt,
            forwardTransform(tfield().component(cmpt), nn)
        );
    }

    tfield.clear();

    return tfftVectorField;
}


tmp<complexVectorField> fft::reverseTransform
(
    const tmp<complexVectorField>& tfield,
    const UList<int>& nn
)
{
    tmp<complexVectorField> tifftVectorField
    (
        new complexVectorField
        (
            tfield().size()
        )
    );

    for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
    {
        tifftVectorField.ref().replace
        (
            cmpt,
            reverseTransform(tfield().component(cmpt), nn)
        );
    }

    tfield.clear();

    return tifftVectorField;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
