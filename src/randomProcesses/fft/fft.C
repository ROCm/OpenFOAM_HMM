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

Foam::tmp<Foam::complexField>
Foam::fft::realTransform1D(const scalarField& field)
{
    const label n = field.size();
    const label nBy2 = n/2;

    // Copy of input field for use by fftw
    // - require non-const access to input and output
    // - use double to avoid additional libfftwf for single-precision

    List<double> in(n);
    List<double> out(n);

    for (label i=0; i < n; ++i)
    {
        in[i] = field[i];
    }

    // Using real to half-complex fftw 'kind'
    fftw_plan plan = fftw_plan_r2r_1d
    (
        n,
        in.data(),
        out.data(),
        FFTW_R2HC,
        FFTW_ESTIMATE
    );

    fftw_execute(plan);

    // field[0] = DC component
    tmp<complexField> tresult(new complexField(nBy2 + 1));
    complexField& result = tresult.ref();

    result[0].Re() = out[0];
    result[nBy2].Re() = out[nBy2];
    for (label i = 1; i < nBy2; ++i)
    {
        result[i].Re() = out[i];
        result[i].Im() = out[n - i];
    }

    fftw_destroy_plan(plan);

    return tresult;
}


Foam::tmp<Foam::complexField> Foam::fft::realTransform1D
(
    const tmp<scalarField>& tfield
)
{
    tmp<complexField> tresult = realTransform1D(tfield());
    tfield.clear();
    return tresult;
}


void Foam::fft::transform
(
    complexField& field,
    const UList<int>& nn,
    transformDirection dir
)
{
    // Copy field into fftw containers
    const label N = field.size();
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
    const label rank = nn.size();
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

Foam::tmp<Foam::complexField> Foam::fft::forwardTransform
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


Foam::tmp<Foam::complexField> Foam::fft::reverseTransform
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


Foam::tmp<Foam::complexVectorField> Foam::fft::forwardTransform
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


Foam::tmp<Foam::complexVectorField> Foam::fft::reverseTransform
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


// ************************************************************************* //
