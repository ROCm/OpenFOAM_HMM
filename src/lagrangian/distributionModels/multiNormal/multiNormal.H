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
    Foam::distributionModels::multiNormal

Description
    A multiNormal distribution model

    \verbatim
        model = sum_i strength_i * exp(-0.5*((x - expectation_i)/variance_i)^2 )
    \endverbatim


SourceFiles
    multiNormal.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_multiNormal_H
#define distributionModels_multiNormal_H

#include "distributionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                         Class multiNormal Declaration
\*---------------------------------------------------------------------------*/

class multiNormal
:
    public distributionModel
{
    // Private Data

        //- List of means of the parent general normal distributions
        List<scalar> mu_;

        //- List of standard deviations of
        //- the parent general normal distributions
        List<scalar> sigma_;

        //- List of weights of a given distribution in the mixture
        List<scalar> weight_;


public:

    //- Runtime type information
    TypeName("multiNormal");


    // Constructors

        //- Construct from components
        multiNormal(const dictionary& dict, Random& rndGen);

        //- Copy construct
        multiNormal(const multiNormal& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new multiNormal(*this));
        }

        //- No copy assignment
        void operator=(const multiNormal&) = delete;


    //- Destructor
    virtual ~multiNormal() = default;


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
