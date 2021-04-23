/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2017-2021 OpenCFD Ltd.
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
    Foam::metisDecomp

Description
    Metis domain decomposition

    When run in parallel will collect the entire graph on to the master,
    decompose and send back.

    Coefficients dictionary: \a metisCoeffs, \a coeffs.

    \verbatim
    numberOfSubdomains   N;
    method               metis;

    metisCoeffs
    {
        method           recursive;  // k-way
        options          ( ...);
        processorWeights ( ... );
    }
    \endverbatim

    Method coefficients:
    \table
        Property  | Description                 | Required | Default
        method    | recursive / k-way           | no | recursive
        options   | metis options               | no
        processorWeights | list of weighting per partition  | no
    \endtable

SourceFiles
    metisDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef metisDecomp_H
#define metisDecomp_H

#include "metisLikeDecomp.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class metisDecomp Declaration
\*---------------------------------------------------------------------------*/

class metisDecomp
:
    public metisLikeDecomp
{
protected:

    // Protected Member Functions

        //- Decompose non-parallel
        virtual label decomposeSerial
        (
            const labelList& adjncy,
            const labelList& xadj,
            const List<scalar>& cellWeights,
            labelList& decomp
        ) const;


        //- No copy construct
        metisDecomp(const metisDecomp&) = delete;

        //- No copy assignment
        void operator=(const metisDecomp&) = delete;


public:

    //- Runtime type information
    TypeName("metis");


    // Constructors

        //- Construct given decomposition dictionary and optional region name
        explicit metisDecomp
        (
            const dictionary& decompDict,
            const word& regionName = ""
        );


    //- Destructor
    virtual ~metisDecomp() = default;


    // Member Functions

        virtual bool parallelAware() const
        {
            return true;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
