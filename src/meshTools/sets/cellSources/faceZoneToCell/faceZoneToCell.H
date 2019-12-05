/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::faceZoneToCell

Description
    A topoSetCellSource to select cells based on side of faceZone.

    \heading Dictionary parameters
    \table
        Property    | Description                           | Required | Default
        option      | Selection type (master / slave)       | yes   |
        zones       | The face zone names or regexs         | possibly |
        zone        | The face zone name or regex           | possibly |
        name        | Older specification for 'zone'        | no    |
    \endtable

Note
    Must specify "zones", "zone" or "name" (highest to lowest precedence).

SourceFiles
    faceZoneToCell.C

\*---------------------------------------------------------------------------*/

#ifndef faceZoneToCell_H
#define faceZoneToCell_H

#include "topoSetCellSource.H"
#include "wordRes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class faceZoneToCell Declaration
\*---------------------------------------------------------------------------*/

class faceZoneToCell
:
    public topoSetCellSource
{
public:
        //- Enumeration defining the valid options
        enum faceAction
        {
            MASTER,
            SLAVE
        };

private:

    // Private data

        static const Enum<faceAction> faceActionNames_;

        //- Add usage string
        static addToUsageTable usage_;

        //- Matcher for face zones
        wordRes selectedZones_;

        //- Option
        faceAction option_;


    // Private Member Functions

        void combine(topoSet& set, const bool add) const;


public:

    //- Runtime type information
    TypeName("faceZoneToCell");

    // Constructors

        //- Construct from components
        faceZoneToCell
        (
            const polyMesh& mesh,
            const wordRe& zoneName,
            const faceAction option
        );

        //- Construct from dictionary
        faceZoneToCell(const polyMesh& mesh, const dictionary& dict);

        //- Construct from Istream
        faceZoneToCell(const polyMesh& mesh, Istream& is);


    //- Destructor
    virtual ~faceZoneToCell() = default;


    // Member Functions

        virtual void applyToSet
        (
            const topoSetSource::setAction action,
            topoSet& set
        ) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
