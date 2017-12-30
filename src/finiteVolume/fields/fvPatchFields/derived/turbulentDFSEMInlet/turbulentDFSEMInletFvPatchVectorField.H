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
    Foam::turbulentDFSEMInletFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    Velocity boundary condition including synthesised eddies for use with LES
    and DES turbulent flows.

    Reference:
    \verbatim
        Poletto, R., Craft, T., and Revell, A.,
        "A New Divergence Free Synthetic Eddy Method for the Reproduction
         of Inlet Flow Conditions for LES",
        Flow Turbulence Combust (2013) 91:519-539
    \endverbatim

    Reynolds stress, velocity and turbulence length scale values can either
    be sepcified directly, or mapped.  If mapping, the values should be
    entered in the same form as the timeVaryingMappedFixedValue condition,
    except that no interpolation in time is supported.  These should be
    located in the directory:

        \$FOAM_CASE/constant/boundaryData/<patchName>/points
        \$FOAM_CASE/constant/boundaryData/<patchName>/0/\{R|U|L\}


Usage
    \table
        Property  | Description         | Required         | Default value
        value     | Restart value       | yes              |
        delta     | Local limiting length scale  | yes     |
        R         | Reynolds stress field | no             |
        U         | Velocity field      | no               |
        L         | Turbulence length scale field | no     |
        d         | Eddy density(fill fraction) | no       | 1
        kappa     | Von Karman constant | no               | 0.41
        mapMethod | Method to map reference values | no | planarInterpolation
        perturb   | Point perturbation for interpolation | no | 1e-5
        writeEddies | Flag to write eddies as OBJ file | no | no
    \endtable

Note
    - The \c delta value typically represents a channel half-height
    - For R, U and L specification: if the entry is not user input, it is
      assumed that the data will be mapped

SeeAlso
    timeVaryingMappedFixedValueFvPatchField

SourceFiles
    turbulentDFSEMInletFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef turbulentDFSEMInletFvPatchVectorField_H
#define turbulentDFSEMInletFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Random.H"
#include "eddy.H"
#include "pointIndexHit.H"
#include "instantList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class cartesianCS;
class pointToPointPlanarInterpolation;

/*---------------------------------------------------------------------------*\
           Class turbulentDFSEMInletFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class turbulentDFSEMInletFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Maximum number of attempts when seeding eddies
        static label seedIterMax_;

        //- Typical length scale, e.g. half channel height
        const scalar delta_;

        //- Ratio of sum of eddy volumes to eddy box volume; default = 1
        const scalar d_;

        //- Von Karman constant
        const scalar kappa_;

        //- Global numbering for faces
        mutable autoPtr<globalIndex> globalFacesPtr_;


        // Table reading for patch inlet flow properties

            //- Fraction of perturbation (fraction of bounding box) to add
            scalar perturb_;

            //- Interpolation scheme to use
            word mapMethod_;

            //- 2D interpolation (for 'planarInterpolation' mapMethod)
            mutable autoPtr<pointToPointPlanarInterpolation> mapperPtr_;

            //- Flag to identify to interpolate the R field
            bool interpolateR_;

            //- Reynolds stress tensor
            symmTensorField R_;

            //- Flag to identify to interpolate the L field
            bool interpolateL_;

            //- Length scale
            scalarField L_;

            //- Flag to identify to interpolate the U field
            bool interpolateU_;

            //- Inlet velocity
            vectorField U_;

            //- Mean inlet velocity
            vector UMean_;


        // Patch information

            //- Patch area - total across all processors
            scalar patchArea_;

            //- Decomposed patch faces as a list of triangles
            faceList triFace_;

            //- Addressing from per triangle to patch face
            labelList triToFace_;

            //- Cumulative triangle area per triangle face
            scalarList triCumulativeMagSf_;

            //- Cumulative area fractions per processor
            scalarList sumTriMagSf_;


        //- List of eddies
        List<eddy> eddies_;

        //- Minimum number of cells required to resolve an eddy
        label nCellPerEddy_;

        //- Patch normal into the domain
        vector patchNormal_;

        //- Eddy box volume
        scalar v0_;

        //- Random number generator
        Random rndGen_;

        //- Length scale per patch face
        scalarField sigmax_;

        //- Maximum length scale (across all processors)
        scalar maxSigmaX_;

        //- Global number of eddies
        label nEddy_;

        //- Current time index (used for updating)
        label curTimeIndex_;

        //- Patch bounds (local processor)
        boundBox patchBounds_;

        //- Single processor contains all eddies (flag)
        bool singleProc_;

        //- Flag to write the eddies to file
        bool writeEddies_;


    // Private Member Functions

        //- Initialise info for patch point search
        void initialisePatch();

        //- Initialise the eddy box
        void initialiseEddyBox();

        //- Set a new eddy position
        pointIndexHit setNewPosition(const bool global);

        //- Initialise eddies
        void initialiseEddies();

        //- Convect the eddies
        void convectEddies(const scalar deltaT);

        //- Calculate the velocity fluctuation at a point
        vector uDashEddy(const List<eddy>& eddies, const point& globalX) const;

        //- Helper function to interpolate values from the boundary data or
        //  read from dictionary
        template<class Type>
        tmp<Field<Type>> interpolateOrRead
        (
            const word& fieldName,
            const dictionary& dict,
            bool& interpolateField
        ) const;

        //- Helper function to interpolate values from the boundary data
        template<class Type>
        tmp<Field<Type>> interpolateBoundaryData
        (
            const word& fieldName
        ) const;

        //- Write Lumley coefficients to file
        void writeLumleyCoeffs() const;

        //- Write eddy info in OBJ format
        void writeEddyOBJ() const;

        //- Return a reference to the patch mapper object
        const pointToPointPlanarInterpolation& patchMapper() const;

        //- Return eddies from remote processors that interact with local
        //  processor
        void calcOverlappingProcEddies
        (
            List<List<eddy>>& overlappingEddies
        ) const;


public:

   //- Runtime type information
   TypeName("turbulentDFSEMInlet");


   // Constructors

        //- Construct from patch and internal field
        turbulentDFSEMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        turbulentDFSEMInletFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given turbulentDFSEMInletFvPatchVectorField
        //  onto a new patch
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDFSEMInletFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        turbulentDFSEMInletFvPatchVectorField
        (
            const turbulentDFSEMInletFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new turbulentDFSEMInletFvPatchVectorField(*this, iF)
            );
        }

    virtual ~turbulentDFSEMInletFvPatchVectorField();


    // Member functions


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap(const fvPatchFieldMapper& m);

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField& ptf,
                const labelList& addr
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "turbulentDFSEMInletFvPatchVectorFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
