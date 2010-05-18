/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "processorPolyPatch.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::processorPolyPatch::transform(Field<T>& l) const
{
    if (l.size() != size())
    {
        FatalErrorIn("processorPolyPatch::transform(Field<T>&) const")
            << "Size of field " << l.size() << " differs from patch size "
            << size() << abort(FatalError);
    }

    forAll(patchIDs_, subI)
    {
        label patchI = patchIDs_[subI];

        if (patchI != -1)
        {
            // Get field on patch
            typename Field<T>::subField subFld(subSlice(l, subI));

            refCast<const coupledPolyPatch>
            (
                boundaryMesh()[patchI]
            ).transform(subFld);
        }
    }
}


// ************************************************************************* //
