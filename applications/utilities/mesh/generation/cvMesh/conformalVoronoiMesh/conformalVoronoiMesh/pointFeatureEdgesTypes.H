/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
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

    As a special exception, you have permission to link this program with the
    CGAL library and distribute executables, as long as you follow the
    requirements of the GNU GPL in regard to all of the software in the
    executable aside from CGAL.

Class
    pointFeatureEdgesTypes

Description
    struct for holding information on the types of feature edges attached to
    feature points

\*---------------------------------------------------------------------------*/

#ifndef pointFeatureEdgesTypes_H
#define pointFeatureEdgesTypes_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class pointFeatureEdgesTypes Declaration
\*---------------------------------------------------------------------------*/

//- Hold the types of feature edges attached to the point.
struct pointFeatureEdgesTypes
{
    label ptI;
    label nExternal, nInternal;
    label nFlat, nOpen, nMultiple, nNonFeature;

    pointFeatureEdgesTypes(const label pointLabel)
    {
        ptI = pointLabel;
        nExternal = nInternal = 0;
        nFlat = nOpen = nMultiple = nNonFeature = 0;
    }

    friend Ostream& operator<<(Ostream& os, const pointFeatureEdgesTypes& p)
    {
        os  << "E " << "I " << "F " << "O " << "M " << "N "
            << "Pt" << nl
            << p.nExternal << " " << p.nInternal << " "
            << p.nFlat << " " << p.nOpen << " "
            << p.nMultiple << " " << p.nNonFeature << " "
            << p.ptI
            << endl;

        return os;
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
