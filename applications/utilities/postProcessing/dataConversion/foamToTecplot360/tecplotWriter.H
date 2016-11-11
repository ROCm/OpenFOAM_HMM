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
    Foam::tecplotWriter

Description
    Write binary tecplot files using tecio.

Note
    The tecplot API uses pass by reference for all routines.
    Its standard integer is defined as INTEGER4 (ie, int32_t),
    which is also used when passing boolean values.

SourceFiles
    tecplotWriter.C
    tecplotWriterTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef tecplotWriter_H
#define tecplotWriter_H

#include "Time.H"
#include "indirectPrimitivePatch.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

/*---------------------------------------------------------------------------*\
                        Class tecplotWriter Declaration
\*---------------------------------------------------------------------------*/

class tecplotWriter
{
    //- Tecplot ZoneTypes
    enum tecplotZoneType
    {
        ZONE_ORDERED         = 0,
        ZONE_FELINESEG       = 1,
        ZONE_FETRIANGLE      = 2,
        ZONE_FEQUADRILATERAL = 3,
        ZONE_FETETRAHEDRON   = 4,
        ZONE_FEBRICK         = 5,
        ZONE_FEPOLYGON       = 6,
        ZONE_FEPOLYHEDRON    = 7
    };


    // Static data members

        // Values commonly used internally
        static const int32_t tecConst_0;
        static const int32_t tecConst_1;
        static const int32_t tecConst_False;
        static const int32_t tecConst_True;


    // Private data

        //- Time reference. Used for paths and the solution time.
        const Time& time_;


    // Private Member Functions

        template<class GeoField>
        static wordList getNames(const PtrList<const GeoField>&);

        //- Disallow default bitwise copy construct
        tecplotWriter(const tecplotWriter&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const tecplotWriter&) = delete;

public:

    //- Data location types
    enum dataLocation
    {
        CELL_CENTERED = 0,
        NODE_CENTERED = 1
    };


    //- Data file type
    enum dataFileType
    {
        FILETYPE_FULL = 0,
        FILETYPE_GRID = 1,
        FILETYPE_SOLUTION = 2
    };


    // Static data members

        //- Commonly used "X Y Z" string
        static const string XYZ;


    // Constructors

        //- Construct from components
        tecplotWriter(const Time&);


    // Member Functions

        //- Initialize writing
        void writeInit
        (
            const word& name,
            const string& varNames,
            const fileName&,
            const dataFileType fileType
        ) const;

        //- Write mesh as polyhedral zone
        void writePolyhedralZone
        (
            const word& zoneName,
            const int32_t strandID,
            const fvMesh& mesh,
            const UList<int32_t>& varLocArray,
            const int32_t NumFaceNodes
        ) const;

        //- Write surface as polygonal zone
        void writePolygonalZone
        (
            const word& zoneName,
            const int32_t strandID,
            const indirectPrimitivePatch& pp,
            const UList<int32_t>& varLocArray
        ) const;

        //- Write unordered data (or rather 1D ordered)
        void writeOrderedZone
        (
            const word& zoneName,
            const int32_t strandID,
            const label n,
            const UList<int32_t>& varLocArray
        ) const;

        //- Write mesh
        void writeConnectivity(const fvMesh&) const;

        //- Write surface
        void writeConnectivity(const indirectPrimitivePatch& pp) const;

        //- Finalize writing
        void writeEnd() const;


        //- Write generic Field, component-wise
        template<class Type>
        void writeField(const Field<Type>&) const;

        //- Write generic Field, component-wise
        template<class Type>
        void writeField(const tmp<Field<Type>>&) const;

        //- Write all fields listed
        template<class GeoField>
        void writeFields(const PtrList<const GeoField>&) const;


        //- Get either fvPatchField or patchInternalField
        template<class Type>
        static tmp<Field<Type>> getPatchField
        (
            const bool nearCellValue,
            const GeometricField<Type, fvPatchField, volMesh>&,
            const label patchi
        );

        //- Get mixed field: fvsPatchField for boundary faces and
        //  internalField for internal faces.
        template<class Type>
        static tmp<Field<Type>> getFaceField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>&,
            const labelUList& faceLabels
        );


        //- Fill in tecplot names/locations for the given input names
        template<class Type>
        static void getTecplotNames
        (
            const wordList& names,
            const int32_t loc,
            string& varNames,
            DynamicList<int32_t>& varLocation
        );

        //- Fill in tecplot names/locations for the given input fields
        template<class GeoField>
        static void getTecplotNames
        (
            const PtrList<GeoField>& flds,
            const int32_t loc,
            string& varNames,
            DynamicList<int32_t>& varLocation
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "tecplotWriterTemplates.C"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
