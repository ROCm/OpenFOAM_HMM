/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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
    Foam::distributionModels::normal

Description
    A normal distribution model

    \verbatim
        model = strength * exp(-0.5*((x - expectation)/variance)^2 )
    \endverbatim

    strength only has meaning if there's more than one distribution model

SourceFiles
    normal.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_normal_H
#define distributionModels_normal_H

#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                           Class normal Declaration
\*---------------------------------------------------------------------------*/

class normal
:
    public distributionModel
{
    // Private Data

        //- Mean of the parent general normal distribution
        scalar mu_;

        //- Standard deviation of the parent general normal distribution
        scalar sigma_;


public:

    //- Runtime type information
    TypeName("normal");


    // Constructors

        //- Construct from components
        normal(const dictionary& dict, Random& rndGen);

        //- Copy construct
        normal(const normal& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new normal(*this));
        }

        //- No copy assignment
        void operator=(const normal&) = delete;


    //- Destructor
    virtual ~normal() = default;


    // Member Functions

        //- Sample the distribution
        virtual scalar sample() const;

        //- Return the theoretical mean of the distribution
        virtual scalar meanValue() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
