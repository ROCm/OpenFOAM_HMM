/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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
        writeFields     false;

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
        Property     | Description                          | Required | Default
        type         | Type name: surfaceFieldValue         | yes |
        log          | Write data to standard output        | no  | no
        regionType   | Face regionType: see below           | yes |
        name         | Name for regionType                  | yes |
        operation    | Operation to perform                 | yes |
        postOperation | Post-operation to perform           | no  | none
        fields       | List of fields to operate on         | yes |
        weightField  | Name of field to apply weighting     | no  |
        scaleFactor  | Output value scaling factor          | no  | 1
        writeArea    | Write the surface area               | no  |
        writeFields  | Write the region field values        | yes |
        surfaceFormat | Output value format                 | no  | none
    \endtable

    Where \c regionType is defined by
    \plaintable
        faceZone     | The \b name entry to specify the faceZone
        patch        | The \b name entry to specify the patch
        functionObjectSurface | The \b name entry to specify a polySurface
        sampledSurface | A \b sampledSurfaceDict sub-dictionary and \b name
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
       areaNormalAverage | area-weighted average in face normal direction
       areaNormalIntegrate | area-weighted integral in face normal directon
       uniformity    | uniformity index
       weightedSum           | weighted sum
       weightedAverage       | weighted average
       weightedAreaAverage   | weighted area average
       weightedAreaIntegrate | weighted area integral
       weightedUniformity    | weighted uniformity index
       absWeightedSum           | sum using absolute weighting
       absWeightedAverage       | average using absolute weighting
       absWeightedAreaAverage   | area average using absolute weighting
       absWeightedAreaIntegrate | area integral using absolute weighting
       absWeightedUniformity    | uniformity index using absolute weighting
    \endplaintable

Note
    - The values reported by the areaNormalAverage and areaNormalIntegrate
      operations are written as the first component of a field with the same
      rank as the input field.
    - Faces on empty patches get ignored
    - If the field is a volField the \c faceZone can only consist of boundary
      faces
    - Using \c functionObjectSurface:
      - The keyword %subRegion should not be used to select surfaces.
        Instead specify the regionType 'functionObjectSurface' and provide
        the name.
    - Using \c sampledSurface:
        - not available for surface fields
        - if interpolate=true they use \c interpolationCellPoint
          otherwise they use cell values
        - each triangle in \c sampledSurface is logically only in one cell
          so interpolation will be wrong when triangles are larger than
          cells.  This can only happen for sampling on a \c triSurfaceMesh
        - take care when using isoSurfaces - these might have duplicate
          triangles and so integration might be wrong

    Uniformity
    \f[
        UI(\phi) = 1 - \frac{1}{2 \overline{\phi} A}
        \int{\left| W \phi \cdot \hat{n} - \bar{W} \bar{\phi}\right| d\vec{A}}
        \,,\;
        \bar{\phi} = \frac{\int{W \phi \cdot d\vec{A}}}{\int{W \cdot d\vec{A}}}
    \f]

    A velocity uniformity index is calculated with no weighting (W=1) and
    \f$ \phi = \vec{U} \f$.

    A scalar concentration uniformity index is calculated with either
    \f$ \rho \vec{U} \f$ or \f$ \vec{U} \f$ for weighting and
    \f$ \phi = conc \f$.

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
#include "surfaceMesh.H"
#include "polySurface.H"
#include "fvsPatchField.H"
#include "volFieldsFwd.H"
#include "polySurfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
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

    // Public Data Types

        //- Region type enumeration
        enum regionTypes
        {
            stFaceZone = 0x01,   //!< Calculate on a faceZone
            stPatch    = 0x02,   //!< Calculate on a patch
            stObject   = 0x11,   //!< Calculate with function object surface
            stSampled  = 0x12    //!< Sample onto surface and calculate
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
            opUniformity,           //!< Uniformity index

            // Weighted variants

            //! Weighted sum
            opWeightedSum = (opSum | typeWeighted),

            //! Weighted average
            opWeightedAverage = (opAverage | typeWeighted),

            //! Weighted area average
            opWeightedAreaAverage = (opAreaAverage | typeWeighted),

            //! Weighted area integral
            opWeightedAreaIntegrate = (opAreaIntegrate | typeWeighted),

            //! Weighted uniformity index
            opWeightedUniformity = (opUniformity | typeWeighted),

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

            //! Uniformity index using abs weighting
            opAbsWeightedUniformity =
                (opWeightedUniformity | typeAbsolute),
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

    // Protected Data

        //- Region type
        regionTypes regionType_;

        //- Operation to apply to values
        operationType operation_;

        //- Optional post-evaluation operation
        postOperationType postOperation_;

        //- Weight field name - optional
        word weightFieldName_;

        //- Track if the surface needs an update
        bool needsUpdate_;

        //- Optionally write the area of the surfaceFieldValue
        bool writeArea_;

        //- Total area of the surfaceFieldValue
        scalar totalArea_;

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

        //- The sampledSurface (when operating on sampledSurface)
        autoPtr<sampledSurface> sampledPtr_;

        //- Surface writer
        autoPtr<surfaceWriter> surfaceWriterPtr_;


    // Protected Member Functions

        //- The volume mesh or surface registry being used
        const objectRegistry& obr() const;

        //- Can the surface definition sample surface-fields?
        inline bool withSurfaceFields() const;

        //- Can use mesh topological merge?
        inline bool withTopologicalMerge() const;

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

        //- Update the surface and surface information as required.
        //  Do nothing (and return false) if no update was required
        bool update();

        //- Return true if the field name is known and a valid type
        template<class Type>
        bool validField(const word& fieldName) const;

        //- Return field values by looking up field name
        template<class Type>
        tmp<Field<Type>> getFieldValues
        (
            const word& fieldName,
            const bool mandatory = false
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
            const pointField& points,
            const faceList& faces
        );

        //- Templated helper function to output field values
        template<class Type, class WeightType>
        bool writeValues
        (
            const word& fieldName,
            const vectorField& Sf,
            const Field<WeightType>& weightField,
            const pointField& points,
            const faceList& faces
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
    virtual ~surfaceFieldValue() = default;


    // Public Member Functions

        //- Return the region type
        inline regionTypes regionType() const;

        //- Return the output directory
        inline fileName outputDir() const;

        //- Read from dictionary
        virtual bool read(const dictionary& dict);

        //- Calculate and write
        virtual bool write();

        //- Update for changes of mesh
        virtual void updateMesh(const mapPolyMesh& mpm);

        //- Update for changes of mesh
        virtual void movePoints(const polyMesh& mesh);
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
