/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

Class
    Foam::vtkTools

Description
    Various functions for collecting and writing binary data.

SourceFiles
    vtkTools.C

\*---------------------------------------------------------------------------*/

#ifndef vtkTools_H
#define vtkTools_H

#include "floatScalar.H"
#include "DynamicList.H"
#include "volFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class vtkTools Declaration
\*---------------------------------------------------------------------------*/

class vtkTools
{
    // Private Member Functions

        // Swap halves of word.

            static void swapWord(label& word32);
            static void swapWords(const label nWords, label* words32);


public:

    // Write ascii or binary. If binary optionally in-place swaps argument

        static void write(std::ostream&, const bool, List<floatScalar>&);
        static void write(std::ostream&, const bool, DynamicList<floatScalar>&);
        static void write(std::ostream&, const bool, labelList&);
        static void write(std::ostream&, const bool, DynamicList<label>&);


    // Write header

        static void writeHeader
        (
            std::ostream&,
            const bool isBinary,
            const std::string& title
        );
        static void writeCellDataHeader
        (
            std::ostream&,
            const label nCells,
            const label nFields
        );
        static void writePointDataHeader
        (
            std::ostream&,
            const label nPoints,
            const label nFields
        );


    // Convert to VTK and store

        static void insert(const scalar, DynamicList<floatScalar>&);
        static void insert(const point&, DynamicList<floatScalar>&);
        static void insert(const sphericalTensor&, DynamicList<floatScalar>&);
        static void insert(const symmTensor&, DynamicList<floatScalar>&);
        static void insert(const tensor&, DynamicList<floatScalar>&);


    //- Append elements to DynamicList
    static void insert(const labelList&, DynamicList<label>&);

    template<class Type>
    static void insert(const List<Type>&, DynamicList<floatScalar>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "vtkToolsTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
