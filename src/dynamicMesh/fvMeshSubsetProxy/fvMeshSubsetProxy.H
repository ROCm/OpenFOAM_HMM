/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2018 OpenCFD Ltd.
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
    Foam::fvMeshSubsetProxy

Description
    Simple proxy for holding a mesh, or mesh-subset.
    The subMeshes are currently limited to cellSet or cellZone definitions.

SourceFiles
    fvMeshSubsetProxy.C
    fvMeshSubsetProxyTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshSubsetProxy_H
#define fvMeshSubsetProxy_H

#include "fvMeshSubset.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class Time;

/*---------------------------------------------------------------------------*\
                      Class fvMeshSubsetProxy Declaration
\*---------------------------------------------------------------------------*/

class fvMeshSubsetProxy
{
public:

        //- Internal book-keeping for subset type
        enum subsetType
        {
            NONE = 0,   //<! Not a subset
            SET,        //<! Using a cellSet for the subset
            ZONE        //<! Using a cellZone for the subset
        };


private:

    // Private data

        //- Reference to mesh
        fvMesh& baseMesh_;

        //- Subsetting engine + sub-fvMesh
        fvMeshSubset subsetter_;

        //- Subset type if any.
        const subsetType type_;

        //- Name of current cellSet/cellZone (or empty)
        const word name_;

        //- Patch ID for exposed internal faces
        const label exposedPatchId_;


    // Private Member Functions

        //- No copy construct
        fvMeshSubsetProxy(const fvMeshSubsetProxy&) = delete;

        //- No copy assignment
        void operator=(const fvMeshSubsetProxy&) = delete;


public:

    // Constructors

        //- Construct from components
        fvMeshSubsetProxy(fvMesh& baseMesh);

        //- Construct from components
        fvMeshSubsetProxy
        (
            fvMesh& baseMesh,
            const subsetType,
            const word& name,
            const label exposedPatchId = -1
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

        //- Check if a sub-mesh is being used
        inline bool useSubMesh() const
        {
            return type_ != NONE;
        }

        //- Access either base-mesh or sub-mesh
        inline const fvMesh& mesh() const
        {
            if (useSubMesh())
            {
                return subsetter_.subMesh();
            }

            return baseMesh_;
        }

        //- The associated (set or zone) name if any.
        inline const word& name() const
        {
            return name_;
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
            >::Internal& fld
        );


        //- Wrapper for field or the subsetted field.
        //  Map volume field (does in fact do very little interpolation;
        //  just copied from fvMeshSubset)
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        interpolate
        (
            const GeometricField<Type, fvPatchField, volMesh>& fld
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
            >::Internal& fld
        ) const;


        //- Map volume field. Just forwards to fvMeshSubset
        template<class GeoField>
        tmp<GeoField> interpolate(const GeoField& fld) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "fvMeshSubsetProxyTemplates.C"
#endif

#endif

// ************************************************************************* //
