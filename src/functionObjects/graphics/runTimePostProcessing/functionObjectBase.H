/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::runTimePostPro::functionObjectBase

Description
    Base class for function object visualisation

SourceFiles
    functionObjectBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_functionObjectBase_H
#define functionObjects_runTimePostPro_functionObjectBase_H

#include "fieldVisualisationBase.H"
#include "stateFunctionObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
namespace runTimePostPro
{

/*---------------------------------------------------------------------------*\
                     Class functionObjectBase Declaration
\*---------------------------------------------------------------------------*/

class functionObjectBase
:
    public fieldVisualisationBase
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        functionObjectBase(const functionObjectBase&);

        //- Disallow default bitwise assignment
        void operator=(const functionObjectBase&);


protected:

    // Protected data

        //- Reference to the state
        const stateFunctionObject& state_;

        //- Function object name
        word functionObjectName_;

        //- Flag to indicate that source data should be cleared after use
        bool clearObjects_;


    // Protected Member Functions

        //- Retrieve file used to create the scene object
        fileName getFileName(const word& keyword, const word& fieldName) const;

        //- Remove file used to create the scene object
        bool removeFile(const word& keyword, const word& fieldName);


public:

    //- Run-time type information
    TypeName("functionObjectBase");


    // Constructors

        //- Construct from dictionary
        functionObjectBase
        (
            const stateFunctionObject& state,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    //- Destructor
    virtual ~functionObjectBase();


    // Member Functions

        //- Clear files used to create the object(s)
        virtual bool clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
