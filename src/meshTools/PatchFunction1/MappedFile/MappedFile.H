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
    Foam::PatchFunction1Types::MappedFile

Description

SourceFiles
    MappedFile.C

\*---------------------------------------------------------------------------*/

#ifndef PatchFunction1Types_MappedFile_H
#define PatchFunction1Types_MappedFile_H

#include "PatchFunction1.H"
#include "pointToPointPlanarInterpolation.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace PatchFunction1Types
{

/*---------------------------------------------------------------------------*\
                     Class MappedFile Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class MappedFile
:
    public PatchFunction1<Type>
{
    // Private data

        //- Whether constructed from dictionary
        const bool dictConstructed_;

        //- Name of the field data table, defaults to the name of the field
        word fieldTableName_;

        //- If true adjust the mapped field to maintain average value
        Switch setAverage_;

        //- Fraction of perturbation (fraction of bounding box) to add
        scalar perturb_;

        //- Name of points file; default = "points"
        word pointsName_;

        //- Interpolation scheme to use
        word mapMethod_;

        //- 2D interpolation (for 'planarInterpolation' mapMethod)
        mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

        //- List of boundaryData time directories
        mutable instantList sampleTimes_;

        //- Current starting index in sampleTimes
        mutable label startSampleTime_;

        //- Interpolated values from startSampleTime
        mutable Field<Type> startSampledValues_;

        //- If setAverage: starting average value
        mutable Type startAverage_;

        //- Current end index in sampleTimes
        mutable label endSampleTime_;

        //- Interpolated values from endSampleTime
        mutable Field<Type> endSampledValues_;

        //- If setAverage: end average value
        mutable Type endAverage_;

        //- Time varying offset values to interpolated data
        autoPtr<Function1<Type>> offset_;


    // Private Member Functions

        void checkTable(const scalar t) const;

        //- No copy assignment
        void operator=(const MappedFile<Type>&) = delete;


public:

    // Runtime type information
    TypeName("mappedFile");


    // Constructors

        //- Construct from entry name and dictionary
        MappedFile
        (
            const polyPatch& pp,
            const word& type,
            const word& entryName,
            const dictionary& dict,
            const bool faceValues = true
        );

        //- Construct from entry name and dictionary
        MappedFile
        (
            const polyPatch& pp,
            const word& entryName,
            const dictionary& dict,
            const word& fieldTableName,
            const bool faceValues
        );

        //- Copy constructor
        explicit MappedFile(const MappedFile<Type>& ut);

        //- Copy constructor setting patch
        explicit MappedFile
        (
            const MappedFile<Type>& ut,
            const polyPatch& pp
        );

        //- Construct and return a clone
        virtual tmp<PatchFunction1<Type>> clone() const
        {
            return tmp<PatchFunction1<Type>>
            (
                new MappedFile<Type>(*this)
            );
        }

        //- Construct and return a clone setting patch
        virtual tmp<PatchFunction1<Type>> clone(const polyPatch& pp) const
        {
            return tmp<PatchFunction1<Type>>
            (
                new MappedFile<Type>(*this, pp)
            );
        }


    //- Destructor
    virtual ~MappedFile() = default;


    // Member Functions

        // Evaluation

            //- Return MappedFile value
            virtual tmp<Field<Type>> value(const scalar) const;

            //- Is value constant (i.e. independent of x)
            virtual bool constant() const
            {
                return sampleTimes_.size() == 1;
            }

            //- Is value uniform (i.e. independent of coordinate)
            virtual bool uniform() const
            {
                return PatchFunction1<Type>::uniform();
            }

            //- Integrate between two values
            virtual tmp<Field<Type>> integrate
            (
                const scalar x1,
                const scalar x2
            ) const;


        // Mapping

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const FieldMapper& mapper);

            //- Reverse map the given PatchFunction1 onto this PatchFunction1
            virtual void rmap
            (
                const PatchFunction1<Type>& pf1,
                const labelList& addr
            );


        // I-O

            //- Write in dictionary format
            virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace PatchFunction1Types
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MappedFile.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
