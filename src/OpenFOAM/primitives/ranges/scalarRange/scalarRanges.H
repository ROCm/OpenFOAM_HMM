/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2018 OpenCFD Ltd.
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
    Foam::scalarRanges

Description
    A collection of scalar bounds to be used as a unary predicate.

SourceFiles
    scalarRanges.C

\*---------------------------------------------------------------------------*/

#ifndef scalarRanges_H
#define scalarRanges_H

#include "scalarRange.H"
#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class scalarRanges Declaration
\*---------------------------------------------------------------------------*/

class scalarRanges
:
    public List<scalarRange>
{
public:

    // Constructors

        //- Construct null
        inline scalarRanges();

        //- Construct by parsing a string for scalar ranges
        //  The individual items are space, comma or semicolon delimited.
        scalarRanges(const std::string& str);


    // Member Operators

        //- Unary predicate to test if the value is within any of the ranges
        inline bool operator()(const scalar& value) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "scalarRangesI.H"

#endif

// ************************************************************************* //
