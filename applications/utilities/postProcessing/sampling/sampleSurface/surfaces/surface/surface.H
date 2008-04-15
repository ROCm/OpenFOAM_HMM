/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::surface

Description

SourceFiles
    surface.C

\*---------------------------------------------------------------------------*/

#ifndef surface_H
#define surface_H

#include "pointField.H"
#include "word.H"
#include "labelList.H"
#include "faceList.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "autoPtr.H"
#include "volFieldsFwd.H"
#include "fieldsCache.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class meshSearch;
class volPointInterpolation;
template<class Type> class fieldsCache;

/*---------------------------------------------------------------------------*\
                           Class surface Declaration
\*---------------------------------------------------------------------------*/

class surface
{
    // Private data

        //- Reference to mesh
        const polyMesh& mesh_;

        //- Reference to mesh searching class
        meshSearch& searchEngine_;

        //- Name of set
        const word name_;

protected:

    // Protected static functions

        //- Read bool from dictionary. Return provided value if not found
        static bool getBool(const dictionary& dict, const word&, const bool);


public:

    //- Runtime type information
    TypeName("surface");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            surface,
            word,
            (
                const polyMesh& mesh,
                meshSearch& searchEngine,
                const dictionary& dict
            ),
            (mesh, searchEngine, dict)
        );


        //- Class used for the read-construction of
        //  PtrLists of surface
        class iNew
        {
            const polyMesh& mesh_;
            meshSearch& searchEngine_;

        public:

            iNew(const polyMesh& mesh, meshSearch& searchEngine)
            :
                mesh_(mesh),
                searchEngine_(searchEngine)
            {}

            autoPtr<surface> operator()(Istream& is) const
            {
                word sampleType(is);
                dictionary dict(is);                    
                return surface::New(sampleType, mesh_, searchEngine_, dict);
            }
        };


    // Constructors

        //- Construct from mesh, name
        surface
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const word& name
        );

        //- Construct from dictionary
        surface
        (
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict
        );

        //- Clone
        autoPtr<surface> clone() const
        {
            notImplemented("autoPtr<surface> clone() const");
            return autoPtr<surface>(NULL);
        }


    // Selectors

        //- Return a reference to the selected surface
        static autoPtr<surface> New
        (
            const word& sampleType,
            const polyMesh& mesh,
            meshSearch& searchEngine,
            const dictionary& dict
        );


    // Destructor

        virtual ~surface();


    // Member Functions

        const polyMesh& mesh() const
        {
            return mesh_;
        }

        meshSearch& searchEngine() const
        {
            return searchEngine_;
        }

        //- Name of surface
        const word& name() const
        {
            return name_;
        }

        //- Points of surface
        virtual const pointField& points() const = 0;

        //- Faces of surface
        virtual const faceList& faces() const = 0;

        //- Correct for mesh movement and/or field changes
        virtual void correct
        (
            const bool meshChanged,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes,
            const fieldsCache<scalar>&
        ) = 0;

        //- interpolate field to surface
        virtual tmp<scalarField> interpolate
        (
            const word&,
            const fieldsCache<scalar>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const = 0;

        //- interpolate field to surface
        virtual tmp<vectorField> interpolate
        (
            const word&,
            const fieldsCache<vector>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const = 0;

        //- interpolate field to surface
        virtual tmp<sphericalTensorField> interpolate
        (
            const word&,
            const fieldsCache<sphericalTensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const = 0;

        //- interpolate field to surface
        virtual tmp<symmTensorField> interpolate
        (
            const word&,
            const fieldsCache<symmTensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const = 0;

        //- interpolate field to surface
        virtual tmp<tensorField> interpolate
        (
            const word&,
            const fieldsCache<tensor>&,
            const volPointInterpolation&,
            const dictionary& interpolationSchemes
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
