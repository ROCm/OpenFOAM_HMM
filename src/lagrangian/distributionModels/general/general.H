/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::distributionModels::general

Description
    general distribution model

SourceFiles
    general.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_general_H
#define distributionModels_general_H

#include "distributionModel.H"
#include "Vector.H"
#include "VectorSpace.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

namespace distributionModels
{
    class general;
}

// Forward declaration of friend functions and operators
Istream& operator>>(Istream&, distributionModels::general&);
Ostream& operator<<(Ostream&, const distributionModels::general&);

namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                           Class general Declaration
\*---------------------------------------------------------------------------*/

class general
:
    public distributionModel
{
    // Private data

        typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

        // List of (bin probability)
        List<pair> xy_;

        scalar meanValue_;

        List<scalar> integral_;


    // Private member functions

        //- Initialise the distribution parameters
        void initialise();


public:

    //- Runtime type information
    TypeName("general");


    // Constructors

        //- Construct from components
        general(const dictionary& dict, Random& rndGen);

        //- Construct from components
        general
        (
            const UList<scalar>& sampleData,
            const scalar binWidth,
            Random& rndGen
        );

        //- Construct copy
        general(const general& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new general(*this));
        }


    //- Destructor
    virtual ~general();


    // Member Functions

        //- Bin boundaries
        virtual tmp<Field<scalar>> x() const;

        //- Probabilities
        virtual tmp<Field<scalar>> y() const;

        //- Sample the distributionModel
        virtual scalar sample() const;

        //- Return the minimum value
        virtual scalar minValue() const;

        //- Return the maximum value
        virtual scalar maxValue() const;

        //- Return the mean value
        virtual scalar meanValue() const;

        //- Write data to stream
        virtual void writeData(Ostream& os) const;

        //- Read data from stream
        virtual void readData(Istream& os);

        //- Write data in dictionary format
        virtual dictionary writeDict(const word& dictName) const;

        //- Read data from dictionary
        virtual void readDict(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace distributionModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
