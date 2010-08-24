/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

Class
    Foam::topoDistanceData

Description
    For use with FaceCellWave. Determines topological distance to starting faces

SourceFiles
    topoDistanceDataI.H
    topoDistanceData.C

\*---------------------------------------------------------------------------*/

#ifndef topoDistanceData_H
#define topoDistanceData_H

#include "point.H"
#include "tensor.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyPatch;
class polyMesh;

/*---------------------------------------------------------------------------*\
                           Class topoDistanceData Declaration
\*---------------------------------------------------------------------------*/

class topoDistanceData
{
    // Private data

        //- Starting data
        label data_;

        //- Distance
        label distance_;


public:

    // Constructors

        //- Construct null
        inline topoDistanceData();

        //- Construct from count
        inline topoDistanceData
        (
            const label data,
            const label distance
        );


    // Member Functions

        // Access


            inline label data() const
            {
                return data_;
            }
            inline label distance() const
            {
                return distance_;
            }


        // Needed by FaceCellWave

            //- Check whether origin has been changed at all or
            //  still contains original (invalid) value.
            inline bool valid() const;

            //- Check for identical geometrical data. Used for cyclics checking.
            inline bool sameGeometry
            (
                const polyMesh&,
                const topoDistanceData&,
                const scalar
            ) const;

            //- Convert any absolute coordinates into relative to (patch)face
            //  centre
            inline void leaveDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Reverse of leaveDomain
            inline void enterDomain
            (
                const polyMesh&,
                const polyPatch&,
                const label patchFaceI,
                const point& faceCentre
            );

            //- Apply rotation matrix to any coordinates
            inline void transform
            (
                const polyMesh&,
                const tensor&
            );

            //- Influence of neighbouring face.
            inline bool updateCell
            (
                const polyMesh&,
                const label thisCellI,
                const label neighbourFaceI,
                const topoDistanceData& neighbourInfo,
                const scalar tol
            );

            //- Influence of neighbouring cell.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const label neighbourCellI,
                const topoDistanceData& neighbourInfo,
                const scalar tol
            );

            //- Influence of different value on same face.
            inline bool updateFace
            (
                const polyMesh&,
                const label thisFaceI,
                const topoDistanceData& neighbourInfo,
                const scalar tol
            );

    // Member Operators

        // Needed for List IO
        inline bool operator==(const topoDistanceData&) const;

        inline bool operator!=(const topoDistanceData&) const;


    // IOstream Operators

        friend Ostream& operator<<(Ostream&, const topoDistanceData&);
        friend Istream& operator>>(Istream&, topoDistanceData&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "topoDistanceDataI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
