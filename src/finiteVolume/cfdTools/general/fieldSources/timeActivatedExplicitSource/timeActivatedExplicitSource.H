/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
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
    Foam::timeActivatedExplicitSourceNew

Description
    Creates a cell set source that is activated at a given time, and remains
    active for a specified duration. The source value is either in specific
    (XX/m3) or absolute (XX) units.

SourceFiles
    timeActivatedExplicitSourceNew.C

\*---------------------------------------------------------------------------*/

#ifndef timeActivatedExplicitSource_H
#define timeActivatedExplicitSource_H

#include "IOdictionary.H"
#include "autoPtr.H"
#include "topoSetSource.H"
#include "cellSet.H"
#include "fvMesh.H"
#include "Time.H"
#include "NamedEnum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Class timeActivatedExplicitSource Declaration
\*---------------------------------------------------------------------------*/

class timeActivatedExplicitSource
:
    public IOdictionary
{
public:

    enum volumeType
    {
        vtSpecific,
        vtAbsolute
    };

    static const NamedEnum<volumeType, 2> volumeTypeNames_;


private:

    // Private member functions

        //- Update the cell set that the source is acting on
        void updateCellSet();


protected:

    // Protected data

        //- Name of the source
        word sourceName_;

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to time database
        const Time& runTime_;


        // Source properties

            //- Dimensions
            const dimensionSet dimensions_;

            //- Volume type
            volumeType volumeType_;

            //- Time start [s]
            scalar timeStart_;

            //- Duration [s]
            scalar duration_;

            //- Value when "on"
            scalar onValue_;

            //- Value when "off"
            scalar offValue_;

            //- Current source value
            scalar currentValue_;


        // Cell set

            //- Cell volumes
            scalarList V_;

            //- Name of cell source (XXXToCell)
            word cellSource_;

            //- Method by which the cells will be selected
            autoPtr<topoSetSource> cellSelector_;

            //- Set of selected cells
            cellSet selectedCellSet_;


    // Protected Member Functions

        //- Disallow default bitwise copy construct
        timeActivatedExplicitSource(const timeActivatedExplicitSource&);

        //- Disallow default bitwise assignment
        void operator=(const timeActivatedExplicitSource&);


public:

    // Constructors

        //- Construct from explicit source name and mesh
        timeActivatedExplicitSource
        (
            const word&,
            const fvMesh&,
            const dimensionSet&
        );


    // Member Functions

        // Access

            //- Return the start time
            scalar timeStart() const;

            //- Return the duration
            scalar duration() const;

            //- Return the current value of the source
            const dimensionedScalar currentValue() const;

            //- Return a tmp field of the source
            virtual tmp<DimensionedField<scalar, volMesh> > Su() const;


        //- Read properties dictionary
        virtual bool read();

        //- Update
        virtual void update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
