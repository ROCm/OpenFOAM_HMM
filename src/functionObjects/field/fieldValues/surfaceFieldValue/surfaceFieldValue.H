/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015-2017 OpenCFD Ltd.
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
    Foam::functionObjects::fieldValues::surfaceFieldValue

Group
    grpFieldFunctionObjects

Description
    Provides a 'face regionType' variant of the fieldValues function object.

    Given a list of user-specified fields and a selection of mesh (or general
    surface) faces, a number of operations can be performed, such as sums,
    averages and integrations.

    For example, to calculate the volumetric or mass flux across a patch,
    apply the 'sum' operator to the flux field (typically \c phi)

Usage
    Examples of function object specification:
    \verbatim
    movingWallPatch
    {
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        regionType      patch;
        name            movingWall;

        operation       areaAverage;
        fields          (p phi U);
    }

    surfaceFieldValue1
    {
        type            surfaceFieldValue;
        libs            ("libfieldFunctionObjects.so");

        log             true;
        writeControl    writeTime;
        writeFields     true;

        surfaceFormat   none;
        regionType      faceZone;
        name            f0;

        operation       sum;
        weightField     alpha1;
        fields          (p phi U);
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: surfaceFieldValue   | yes         |
        log          | write data to standard output | no    | no
        writeFields  | Write the region field values  | yes     |
        writeArea    | Write the area of the surfaceFieldValue | no |
        surfaceFormat | output value format    | no          |
        regionType   | face regionType: see below  | yes         |
        name         | name of face regionType if required  | no |
        operation    | operation to perform    | yes         |
        postOperation | post-operation to perform    | no | none |
        weightField  | name of field to apply weighting | no |
        scaleFactor  | scale factor            | no          | 1
        fields       | list of fields to operate on | yes    |
    \endtable

    Where \c regionType is defined by
    \plaintable
        faceZone     | requires a 'name' entry to specify the faceZone
        patch        | requires a 'name' entry to specify the patch
        surface      | requires a 'name' entry to specify the surfMesh
        sampledSurface | requires a 'sampledSurfaceDict' sub-dictionary
    \endplaintable

    The \c operation is one of:
    \plaintable
       none          | no operation
       min           | minimum
       max           | maximum
       sum           | sum
       sumMag        | sum of component magnitudes
       sumDirection  | sum values that are positive in given direction
       sumDirectionBalance | sum of balance of values in given direction
       average       | ensemble average
       areaAverage   | area-weighted average
       areaIntegrate | area integral
       CoV           | coefficient of variation: standard deviation/mean
       areaNormalAverage| area-weighted average in face normal direction
       areaNormalIntegrate | area-weighted integral in face normal directon
       weightedSum   | weighted sum
       weightedAverage | weighted average
       weightedAreaAverage | weighted area average
       weightedAreaIntegrate | weighted area integral
       absWeightedSum | sum using absolute weighting
       absWeightedAverage | average using absolute weighting
       absWeightedAreaAverage | area average using absolute weighting
       absWeightedAreaIntegrate | area integral using absolute weighting
    \endplaintable

Note
    - The values reported by the areaNormalAverage and areaNormalIntegrate
      operations are written as the first component of a field with the same
      rank as the input field.
    - faces on empty patches get ignored
    - if the field is a volField the \c faceZone can only consist of boundary
      faces
    - Using \c surface:
      - The keyword %subRegion should not be used to select surfaces.
        Specify instead the regionType 'surface' and provide the surface name.
    - using \c sampledSurface:
        - not available for surface fields
        - if interpolate=true they use \c interpolationCellPoint
          otherwise they use cell values
        - each triangle in \c sampledSurface is logically only in one cell
          so interpolation will be wrong when triangles are larger than
          cells.  This can only happen for sampling on a \c triSurfaceMesh
        - take care when using isoSurfaces - these might have duplicate
          triangles and so integration might be wrong

See also
    Foam::fieldValues
    Foam::functionObject

SourceFiles
    surfaceFieldValue.C
    surfaceFieldValueTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_surfaceFieldValue_H
#define functionObjects_surfaceFieldValue_H

#include "fieldValue.H"
#include "Enum.H"
#include "meshedSurf.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "volFieldsFwd.H"
#include "surfFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sampledSurface;
class surfaceWriter;

namespace functionObjects
{
namespace fieldValues
{

/*---------------------------------------------------------------------------*\
                      Class surfaceFieldValue Declaration
\*---------------------------------------------------------------------------*/

class surfaceFieldValue
:
    public fieldValue
{
public:

    // Public data types

        //- Region type enumeration
        enum regionTypes
        {
            stFaceZone,         //!< Calculate on a faceZone
            stPatch,            //!< Calculate on a patch
            stSurface,          //!< Calculate with fields on a surfMesh
            stSampledSurface    //!< Sample onto surface and calculate
        };

        //- Region type names
        static const Enum<regionTypes> regionTypeNames_;

        //- Bitmask values for operation variants
        enum operationVariant
        {
            typeBase = 0,         //!< Base operation
            typeWeighted = 0x100, //!< Operation using weighting
            typeAbsolute = 0x200, //!< Operation using mag (eg, for weighting)
        };

        //- Operation type enumeration
        enum operationType
        {
            // Normal operations

            opNone = 0,             //!< No operation
            opMin,                  //!< Minimum value
            opMax,                  //!< Maximum value
            opSum,                  //!< Sum of values
            opSumMag,               //!< Sum of component magnitudes
            opSumDirection,         //!< Sum in a given direction

            //! Sum of balance of values in given direction
            opSumDirectionBalance,

            opAverage,              //!< Ensemble average
            opAreaAverage,          //!< Area average
            opAreaIntegrate,        //!< Area integral
            opCoV,                  //!< Coefficient of variation
            opAreaNormalAverage,    //!< Area average in normal direction
            opAreaNormalIntegrate,  //!< Area integral in normal direction

            // Weighted variants

            //! Weighted sum
            opWeightedSum = (opSum | typeWeighted),

            //! Weighted average
            opWeightedAverage = (opAverage | typeWeighted),

            //! Weighted area average
            opWeightedAreaAverage = (opAreaAverage | typeWeighted),

            //! Weighted area integral
            opWeightedAreaIntegrate = (opAreaIntegrate | typeWeighted),

            // Variants using absolute weighting

            //! Sum using abs weighting
            opAbsWeightedSum = (opWeightedSum | typeAbsolute),

            //! Average using abs weighting
            opAbsWeightedAverage = (opWeightedAverage | typeAbsolute),

            //! Area average using abs weighting
            opAbsWeightedAreaAverage = (opWeightedAreaAverage | typeAbsolute),

            //! Area integral using abs weighting
            opAbsWeightedAreaIntegrate =
                (opWeightedAreaIntegrate | typeAbsolute),
        };

        //- Operation type names
        static const Enum<operationType> operationTypeNames_;


        //- Post-operation type enumeration
        enum postOperationType
        {
            postOpNone,   //!< No additional operation after calculation
            postOpSqrt    //!< Perform sqrt after normal operation
        };

        //- Operation type names
        static const Enum<postOperationType> postOperationTypeNames_;


private:

    // Private Member Functions

        //- Set faces to evaluate based on a face zone
        void setFaceZoneFaces();

        //- Set faces to evaluate based on a patch
        void setPatchFaces();

        //- Combine mesh faces and points from multiple processors
        void combineMeshGeometry
        (
            faceList& faces,
            pointField& points
        ) const;

        //- Combine surface faces and points from multiple processors
        void combineSurfaceGeometry
        (
            faceList& faces,
            pointField& points
        ) const;

        //- Calculate and return total area of the surfaceFieldValue: sum(magSf)
        scalar totalArea() const;


protected:

    // Protected data

        //- Region type
        regionTypes regionType_;

        //- Operation to apply to values
        operationType operation_;

        //- Optional post-evaluation operation
        postOperationType postOperation_;

        //- Weight field name - optional
        word weightFieldName_;

        //- Total area of the surfaceFieldValue
        scalar totalArea_;

        //- Optionally write the area of the surfaceFieldValue
        bool writeArea_;

        //- Global number of faces
        label nFaces_;


        // If operating on mesh faces (faceZone, patch)

            //- Local list of face IDs
            labelList faceId_;

            //- Local list of patch ID per face
            labelList facePatchId_;

            //- List representing the face flip map
            //  (false: use as-is, true: negate)
            boolList faceFlip_;


        //- Underlying sampledSurface (if operating on sampledSurface)
        autoPtr<sampledSurface> surfacePtr_;

        //- Surface writer
        autoPtr<surfaceWriter> surfaceWriterPtr_;


    // Protected Member Functions

        //- The volume mesh or surface registry being used
        const objectRegistry& obr() const;

        //- Return the local list of face IDs
        inline const labelList& faceId() const;

        //- Return the local list of patch ID per face
        inline const labelList& facePatch() const;

        //- Return the local true/false list representing the face flip map
        inline const boolList& faceFlip() const;

        //- True if the operation needs a surface Sf
        bool usesSf() const;

        //- True if the operation variant uses mag
        inline bool usesMag() const;

        //- True if the operation variant uses a weight-field
        inline bool usesWeight() const;

        //- True if operation variant uses a weight-field that is available.
        //  Checks for availability on any processor.
        template<class WeightType>
        inline bool canWeight(const Field<WeightType>& weightField) const;

        //- Initialise, e.g. face addressing
        void initialise(const dictionary& dict);

        //- Return true if the field name is known and a valid type
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Return field values by looking up field name
        template<class Type>
        tmp<Field<Type>> getFieldValues
        (
            const word& fieldName,
            const bool mustGet = false
        ) const;

        //- Apply the 'operation' to the values. Operation must preserve Type.
        template<class Type, class WeightType>
        Type processSameTypeValues
        (
            const Field<Type>& values,
            const vectorField& Sf,
            const Field<WeightType>& weightField
        ) const;

        //- Apply the 'operation' to the values. Wrapper around
        //  processSameTypeValues. See also template specialisation below.
        template<class Type, class WeightType>
        Type processValues
        (
            const Field<Type>& values,
            const vectorField& Sf,
            const Field<WeightType>& weightField
        ) const;


        //- Filter a surface field according to faceIds
        template<class Type>
        tmp<Field<Type>> filterField
        (
            const GeometricField<Type, fvsPatchField, surfaceMesh>& field
        ) const;

        //- Filter a volume field according to faceIds
        template<class Type>
        tmp<Field<Type>> filterField
        (
            const GeometricField<Type, fvPatchField, volMesh>& field
        ) const;


        //- Weighting factor.
        //  Possibly applies mag() depending on the operation type.
        template<class WeightType>
        tmp<scalarField> weightingFactor
        (
            const Field<WeightType>& weightField
        ) const;

        //- Weighting factor, weight field with the area.
        //  Possibly applies mag() depending on the operation type.
        //  Reverts to mag(Sf) if the weight field is not available.
        template<class WeightType>
        tmp<scalarField> weightingFactor
        (
            const Field<WeightType>& weightField,
            const vectorField& Sf
        ) const;


        //- Templated helper function to output field values
        template<class WeightType>
        label writeAll
        (
            const vectorField& Sf,
            const Field<WeightType>& weightField,
            const meshedSurf& surfToWrite
        );

        //- Templated helper function to output field values
        template<class Type, class WeightType>
        bool writeValues
        (
            const word& fieldName,
            const vectorField& Sf,
            const Field<WeightType>& weightField,
            const meshedSurf& surfToWrite
        );


        //- Output file header information
        virtual void writeFileHeader(Ostream& os) const;


public:

    //- Run-time type information
    TypeName("surfaceFieldValue");


    // Constructors

        //- Construct from name, Time and dictionary
        surfaceFieldValue
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from name, objectRegistry and dictionary
        surfaceFieldValue
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );


    //- Destructor
    virtual ~surfaceFieldValue();


    // Public Member Functions

        //- Return the region type
        inline const regionTypes& regionType() const;

        //- Return the output directory
        inline fileName outputDir() const;

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate and write
        virtual bool write();
};


//- Specialisation for scalar fields
template<>
scalar surfaceFieldValue::processValues
(
    const Field<scalar>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const;


//- Specialisation for vector fields
template<>
vector surfaceFieldValue::processValues
(
    const Field<vector>& values,
    const vectorField& Sf,
    const scalarField& weightField
) const;


//- Specialisation for scalar - pass through
template<>
tmp<scalarField> surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField
) const;


//- Specialisation for scalar - scalar * Area
template<>
tmp<scalarField> surfaceFieldValue::weightingFactor
(
    const Field<scalar>& weightField,
    const vectorField& Sf
) const;


//- Specialisation for vector - vector (dot) Area
template<>
tmp<scalarField> surfaceFieldValue::weightingFactor
(
    const Field<vector>& weightField,
    const vectorField& Sf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fieldValues
} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "surfaceFieldValueI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "surfaceFieldValueTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
