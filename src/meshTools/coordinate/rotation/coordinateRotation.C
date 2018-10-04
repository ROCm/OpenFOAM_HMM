/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "coordinateRotation.H"
#include "dictionary.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coordinateRotation, 0);
    defineRunTimeSelectionTable(coordinateRotation, dictionary);
    defineRunTimeSelectionTable(coordinateRotation, objectRegistry);
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::vector Foam::coordinateRotation::findOrthogonal(const vector& axis)
{
    direction maxCmpt = 0;
    scalar maxVal = mag(axis[maxCmpt]);

    for (direction cmpt=1; cmpt < vector::nComponents; ++cmpt)
    {
        const scalar val = mag(axis[cmpt]);

        if (maxVal < val)
        {
            maxVal  = val;
            maxCmpt = cmpt;
        }
    }

    direction cmpt = ((maxCmpt == vector::nComponents-1) ? 0 : (maxCmpt+1));

    vector dirn(Zero);
    dirn.component(cmpt) = ((axis[maxCmpt] < 0) ? -1 : 1);

    return dirn;
}


Foam::symmTensor Foam::coordinateRotation::transformPrincipal
(
    const tensor& tt,
    const vector& v
)
{
    return symmTensor
    (
        tt.xx()*v.x()*tt.xx()
      + tt.xy()*v.y()*tt.xy()
      + tt.xz()*v.z()*tt.xz(),

        tt.xx()*v.x()*tt.yx()
      + tt.xy()*v.y()*tt.yy()
      + tt.xz()*v.z()*tt.yz(),

        tt.xx()*v.x()*tt.zx()
      + tt.xy()*v.y()*tt.zy()
      + tt.xz()*v.z()*tt.zz(),

        tt.yx()*v.x()*tt.yx()
      + tt.yy()*v.y()*tt.yy()
      + tt.yz()*v.z()*tt.yz(),

        tt.yx()*v.x()*tt.zx()
      + tt.yy()*v.y()*tt.zy()
      + tt.yz()*v.z()*tt.zz(),

        tt.zx()*v.x()*tt.zx()
      + tt.zy()*v.y()*tt.zy()
      + tt.zz()*v.z()*tt.zz()
    );
}


void Foam::coordinateRotation::write(Ostream& os) const
{
     os.writeEntry("e1", e1());
     os.writeEntry("e2", e2());
     os.writeEntry("e3", e3());
}


// ************************************************************************* //
