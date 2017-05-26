/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::meshSubsetHelper

Description
    Simple helper to hold a mesh or mesh-subset and provide uniform access.

SourceFiles
    meshSubsetHelper.C
    meshSubsetHelperTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef meshSubsetHelper_H
#define meshSubsetHelper_H

#include "fvMeshSubset.H"
#include "zeroGradientFvPatchField.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Time;

/*---------------------------------------------------------------------------*\
                      Class meshSubsetHelper Declaration
\*---------------------------------------------------------------------------*/

class meshSubsetHelper
{
public:

        //- Internal book-keeping for subset type
        enum subsetType
        {
            NONE,       //<! Not a subset
            SET,        //<! Using a cellSet for the subset
            ZONE        //<! Using a cellZone for the subset
        };


    // Private data

        //- Reference to mesh
        fvMesh& baseMesh_;

        //- Subsetting engine + sub-fvMesh
        fvMeshSubset subsetter_;

        //- Subset type if any.
        const subsetType type_;

        //- Name of current cellSet/cellZone (or empty)
        const word name_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        meshSubsetHelper(const meshSubsetHelper&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const meshSubsetHelper&) = delete;


public:

    // Constructors

        //- Construct from components
        meshSubsetHelper(fvMesh& baseMesh);

        //- Construct from components
        meshSubsetHelper
        (
            fvMesh& baseMesh,
            const subsetType,
            const word& name
        );


    // Member Functions

    // Access

        //- The entire base mesh
        inline const fvMesh& baseMesh() const
        {
            return baseMesh_;
        }

        //- The mesh subsetter
        inline const fvMeshSubset& subsetter() const
        {
            return subsetter_;
        }

        //- Check if running a sub-mesh is being used
        inline bool useSubMesh() const
        {
            return type_ != NONE;
        }

        //- Access either mesh or submesh
        inline const fvMesh& mesh() const
        {
            if (useSubMesh())
            {
                return subsetter_.subMesh();
            }
            else
            {
                return baseMesh_;
            }
        }


    // Edit

        //- Update mesh subset
        void correct(bool verbose = false);

        //- Read mesh. Correct on topo-change
        polyMesh::readUpdateState readUpdate();


        //- Construct volField (with zeroGradient) from an internal field
        template<class Type>
        static tmp<GeometricField<Type, fvPatchField, volMesh>>
        zeroGradientField
        (
            const typename GeometricField
            <
                Type,
                fvPatchField,
                volMesh
            >::Internal&
        );


        //- Wrapper for field or the subsetted field.
        //  Map volume field (does in fact do very little interpolation;
        //  just copied from fvMeshSubset)
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>&
        ) const;


        //- Convert an internal field to a volume field
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const typename GeometricField
            <
                Type,
                fvPatchField,
                volMesh
            >::Internal&
        ) const;


        //- Map volume field (does in fact do very little interpolation;
        //  just copied from fvMeshSubset)
        template<class GeoField>
        tmp<GeoField> interpolate(const GeoField& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "meshSubsetHelperTemplates.C"
#endif


#endif

// ************************************************************************* //
