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
    Foam::coordinateScaling

Description
    Helper class to wrap coordinate system and component-wise scaling

See also
    Foam::Function1Types

SourceFiles
    coordinateScaling.C

\*---------------------------------------------------------------------------*/

#ifndef coordinateScaling_H
#define coordinateScaling_H

#include "coordinateSystem.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class coordinateScaling Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class coordinateScaling
{
    // Private data

        //- Local co-ordinate system
        const autoPtr<coordinateSystem> coordSys_;

        //- In local co-ordinate system component-wise scaling
        PtrList<Function1<Type>> scale_;

        //- Cached whether any scaling or coordinate system
        bool active_;


    // Private Member Functions

        //- No copy assignment
        void operator=(const coordinateScaling<Type>&) = delete;

public:

    // Constructors

        //- Construct null
        coordinateScaling();

        //- Construct from registry and dictionary
        coordinateScaling
        (
            const objectRegistry&,
            const dictionary&
        );

        //- Construct copy
        coordinateScaling(const coordinateScaling&);


    //- Destructor
    virtual ~coordinateScaling();


    // Member functions

        //- Has any scaling or coordinate transformation
        bool active() const
        {
            return active_;
        }

        //- Access to optional coordinate system
        const autoPtr<coordinateSystem>& coordSys() const
        {
            return coordSys_;
        }

        //- Evaluate
        virtual tmp<Field<Type>> transform
        (
            const pointField& pos,
            const Field<Type>& local
        ) const;

        //- Write dictionary entry
        virtual void writeEntry(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<Field<label>> coordinateScaling<label>::transform
(
    const pointField& pos,
    const Field<label>& local
) const;

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "coordinateScaling.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
