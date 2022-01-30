/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2017 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

Class
    Foam::cylinderAnnulusToFace

Description
    A \c topoSetFaceSource to select all faces whose face
    centre inside a given bounding cylinder annulus.

    Operands:
    \table
      Operand   | Type    | Location
      output    | faceSet | $FOAM_CASE/constant/polyMesh/sets/\<set\>
    \endtable

Usage
    Minimal example by using \c system/topoSetDict.actions:
    \verbatim
    {
        // Mandatory (inherited) entries
        name           <name>;
        type           faceSet;
        action         <action>;

        // Mandatory entries
        source         cylinderAnnulusToFace;
        p1             (<p1X> <p1Y> <p1Z>);
        p2             (<p2X> <p2Y> <p2Z>);
        outerRadius    <radius1>;

        // Optional entries
        innerRadius    <radius2>;
    }
    \endverbatim

    where the entries mean:
    \table
      Property   | Description                         | Type | Req'd | Dflt
      name       | Name of faceSet                     | word |  yes  | -
      type       | Type name: faceSet                  | word |  yes  | -
      action     | Action applied on faces - see below | word |  yes  | -
      source     | Source name: cylinderAnnulusToFace  | word |  yes  | -
      p1          | Coordinate of one of the endpoints | vector |  yes  | -
      p2          | Coordinate of the other endpoint   | vector |  yes  | -
      outerRadius | Cylinder outer radius              | scalar |  yes  | -
      innerRadius | Cylinder inner radius              | scalar |  no   | 0
    \endtable

    Options for the \c action entry:
    \verbatim
      new      | Create a new faceSet from selected faces
      add      | Add selected faces into this faceSet
      subtract | Remove selected faces from this faceSet
    \endverbatim

See also
    - Foam::topoSetSource
    - Foam::topoSetFaceSource
    - Foam::cylinderToFace

SourceFiles
    cylinderAnnulusToFace.C

\*---------------------------------------------------------------------------*/

#ifndef cylinderAnnulusToFace_H
#define cylinderAnnulusToFace_H

#include "topoSetFaceSource.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class cylinderAnnulusToFace Declaration
\*---------------------------------------------------------------------------*/

class cylinderAnnulusToFace
:
    public topoSetFaceSource
{
    // Private Data

        //- Add usage string
        static addToUsageTable usage_;

        //- First point on cylinder axis
        point point1_;

        //- Second point on cylinder axis
        point point2_;

        //- Outer radius
        scalar radius_;

        //- Inner radius
        scalar innerRadius_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("cylinderAnnulusToFace");


    // Constructors

        //- Construct from components
        cylinderAnnulusToFace
        (
            const polyMesh& mesh,
            const point& point1,
            const point& point2,
            const scalar radius,
            const scalar innerRadius = 0
        );

        //- Construct from dictionary
        cylinderAnnulusToFace(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        cylinderAnnulusToFace(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~cylinderAnnulusToFace() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
