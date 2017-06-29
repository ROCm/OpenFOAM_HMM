/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenCFD Ltd.
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
    Foam::distributionModels::binned

Description
    A binned distribution model where the random sample is taken from a
    discrete set of bins

SourceFiles
    binned.C

\*---------------------------------------------------------------------------*/

#ifndef distributionModels_binned_H
#define distributionModels_binned_H

#include "distributionModel.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

namespace distributionModels
{
    class binned;
}

// Forward declaration of friend functions and operators
Istream& operator>>(Istream&, distributionModels::binned&);
Ostream& operator<<(Ostream&, const distributionModels::binned&);

namespace distributionModels
{

/*---------------------------------------------------------------------------*\
                           Class binned Declaration
\*---------------------------------------------------------------------------*/

class binned
:
    public distributionModel
{
    // Private data

        typedef VectorSpace<Vector<scalar>, scalar, 2> pair;

        // List of (bin probability)
        List<pair> xy_;

        //- Distribution mean value
        scalar meanValue_;


    // Private member functions

        //- Initialise the distribution parameters
        void initialise();


public:

    //- Runtime type information
    TypeName("binned");

    static const char* header;


    // Constructors

        //- Construct from dictionary
        binned(const dictionary& dict, Random& rndGen);

        //- Construct from components
        binned
        (
            const UList<scalar>& sampleData,
            const scalar binWidth,
            Random& rndGen
        );

        //- Construct copy
        binned(const binned& p);

        //- Construct and return a clone
        virtual autoPtr<distributionModel> clone() const
        {
            return autoPtr<distributionModel>(new binned(*this));
        }


    //- Destructor
    virtual ~binned();


    // Member Functions

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
