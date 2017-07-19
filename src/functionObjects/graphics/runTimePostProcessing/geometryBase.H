/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
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
    Foam::functionObjects::runTimePostPro::geometryBase

Description
    Base class for surface handling

SourceFiles
    geometryBase.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_runTimePostPro_geometryBase_H
#define functionObjects_runTimePostPro_geometryBase_H

#include "dictionary.H"
#include "vector.H"
#include "Function1.H"
#include "HashPtrTable.H"
#include "Enum.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class vtkRenderer;
class vtkActor;

namespace Foam
{
namespace functionObjects
{

class runTimePostProcessing;

namespace runTimePostPro
{


/*---------------------------------------------------------------------------*\
                        Class geometryBase Declaration
\*---------------------------------------------------------------------------*/

class geometryBase
{

public:

    // Public enumerations

        enum renderModeType
        {
            rmFlat,             //!< Flat shading
            rmGouraud,          //!< Gouraud shading
            rmPhong             //!< Phong shading
        };

        static const Enum<renderModeType> renderModeTypeNames;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        geometryBase(const geometryBase&);

        //- Disallow default bitwise assignment
        void operator=(const geometryBase&);


protected:

    // Protected data

        //- Reference to the parent function object
        const runTimePostProcessing& parent_;

        //- Name
        word name_;

        //- Visible flag
        bool visible_;

        //- Render mode
        renderModeType renderMode_;

        //- Opacity
        autoPtr<Function1<scalar>> opacity_;

        //- Reference to the colours
        const HashPtrTable<Function1<vector>, word>& colours_;


    // Protected functions

        //- Initialse actor
        void initialiseActor(vtkActor* actor) const;


public:

    // Constructors

        //- Construct from dictionary
        geometryBase
        (
            const runTimePostProcessing& parent_,
            const dictionary& dict,
            const HashPtrTable<Function1<vector>, word>& colours
        );


    //- Destructor
    virtual ~geometryBase();


    // Member Functions

        // Access

            //- Return the reference to the parent function object
            const runTimePostProcessing& parent() const;

            //- Return the name
            const word& name() const;

            //- Return the visible flag
            bool visible() const;

            //- Return the opacity
            scalar opacity(const scalar position) const;

            //- Return reference to the colours
            const HashPtrTable<Function1<vector>, word>& colours() const;


        //- Add geometry to scene
        virtual void addGeometryToScene
        (
            const scalar position,
            vtkRenderer* renderer
        ) = 0;

        //- Update the actors
        virtual void updateActors(const scalar position) = 0;

        //- Clear files used to create the object(s)
        virtual bool clear() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace runTimePostPro
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
