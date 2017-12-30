/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::shellSurfaces

Description
    Encapsulates queries for volume refinement ('refine all cells within
    shell').

SourceFiles
    shellSurfaces.C

\*---------------------------------------------------------------------------*/

#ifndef shellSurfaces_H
#define shellSurfaces_H

#include "searchableSurface.H"
#include "Enum.H"
#include "Tuple2.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class searchableSurfaces;

/*---------------------------------------------------------------------------*\
                           Class shellSurfaces Declaration
\*---------------------------------------------------------------------------*/

class shellSurfaces
{
public:

    // Public data types

        //- Volume refinement controls
        enum refineMode
        {
            INSIDE,         // Refine all inside shell
            OUTSIDE,        // ,,         outside
            DISTANCE        // Refine based on distance to shell
        };


private:

    // Private data

        //- Reference to all geometry.
        const searchableSurfaces& allGeometry_;

        //- Indices of surfaces that are shells
        labelList shells_;

        //- Per shell whether to refine inside or outside
        List<refineMode> modes_;

        //- Per shell the list of ranges
        List<scalarField> distances_;

        //- Per shell per distance the refinement level
        labelListList levels_;


        // Gap level refinement

            //- Per shell, per region the small-gap level specification
            List<List<FixedList<label, 3>>> extendedGapLevel_;

            //- Per shell, per region the small-gap level specification
            List<List<volumeType>> extendedGapMode_;


    // Private data

        //- refineMode names
        static const Enum<refineMode> refineModeNames_;


    // Private Member Functions

        //- Helper function for initialisation of levels
        void setAndCheckLevels
        (
            const label shellI,
            const List<Tuple2<scalar, label>>&
        );

        //- Helper function for checking of gap information
        void checkGapLevels
        (
            const dictionary&,
            const label shellI,
            const List<FixedList<label, 3>>& levels
        );

        void orient();

        //- Find first shell with a level higher than maxLevel
        void findHigherLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& maxLevel
        ) const;

        //- Update highest min gap level
        void findHigherGapLevel
        (
            const pointField& pt,
            const labelList& ptLevel,
            const label shellI,
            labelList& gapShell,
            List<FixedList<label, 3>>& gapInfo,
            List<volumeType>& gapMode
        ) const;

        //- Find first shell with a level lower or equal to minLevel. Update
        //  minLevel and shell.
        void findLevel
        (
            const pointField& pt,
            const label shellI,
            labelList& minLevel,
            labelList& shell
        ) const;

public:

    // Constructors

        //- Construct from geometry and dictionary
        shellSurfaces
        (
            const searchableSurfaces& allGeometry,
            const dictionary& shellsDict
        );


    // Member Functions

        // Access

            //- Indices of surfaces that are shells
            const labelList& shells() const
            {
                return shells_;
            }


        // Query

            //- Highest shell level
            label maxLevel() const;

            //- Highest shell gap level
            labelList maxGapLevel() const;

            //- Find shell level higher than ptLevel
            void findHigherLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& maxLevel
            ) const;

            //- Find a shell whose minimum gap level is >= ptLevel
            void findHigherGapLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& gapShell,
                List<FixedList<label, 3>>& gapInfo,
                List<volumeType>& gapMode
            ) const;

            //- Find a shell whose minimum gap level is >= ptLevel. gapInfo
            //  is (0 0 0) if no shell found
            void findHigherGapLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                List<FixedList<label, 3>>& gapInfo,
                List<volumeType>& gapMode
            ) const;

            //- Find first shell (or -1) with level equal or lower than ptLevel.
            void findLevel
            (
                const pointField& pt,
                const labelList& ptLevel,
                labelList& shell
            ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
