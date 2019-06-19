/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2007-2019 PCOpt/NTUA
                            | Copyright (C) 2013-2019 FOSS GP
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
    Foam::sensitivity

Description
    Abstract base class for adjoint sensitivities

SourceFiles
    sensitivity.C

\*---------------------------------------------------------------------------*/

#ifndef sensitivity_H
#define sensitivity_H

#include "volFields.H"
#include "surfaceFields.H"
#include "dictionary.H"
#include "volPointInterpolation.H"

#include "pointMesh.H"
#include "pointPatchField.H"
#include "pointPatchFieldsFwd.H"
#include "fixedValuePointPatchField.H"
#include "boundaryFieldsFwd.H"
#include "createZeroField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sensitivity Declaration
\*---------------------------------------------------------------------------*/

class sensitivity
{
protected:

    // Protected data

        const fvMesh& mesh_;
        dictionary dict_;

        // Cleaner option to go for a labelHashSet. Kept this way for
        // compatibility
        labelList sensitivityPatchIDs_;
        word adjointSolverName_;
        word surfaceFieldSuffix_;
        bool writeAllSurfaceFiles_;

        // autoPtrs for fields holding sensitivities.
        // Not all of them are required for each case

        // Boundary sensitivities at faces. Shape opt & flow control
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Wall face sens w.r.t. (x,y.z)
        autoPtr<boundaryVectorField> wallFaceSensVecPtr_;

        //- Wall face sens projected to normal
        autoPtr<boundaryScalarField> wallFaceSensNormalPtr_;

        //- Normal sens as vectors
        autoPtr<boundaryVectorField> wallFaceSensNormalVecPtr_;

        // Boundary sensitivities at points. Shape opt
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        //- Wall point sens w.r.t. (x,y.z)
        autoPtr<pointBoundaryVectorField> wallPointSensVecPtr_;

        //- Wall point sens projected to normal
        autoPtr<pointBoundaryScalarField> wallPointSensNormalPtr_;

        //- Normal sens as vectors
        autoPtr<pointBoundaryVectorField> wallPointSensNormalVecPtr_;

        //field sensitivities. Topology optimisation
        //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        autoPtr<volScalarField> fieldSensPtr_;


    // Protected Member Functions

        //- Constructs volField based on boundaryField and writes it
        template<class Type>
        void constructAndWriteSensitivityField
        (
            const autoPtr
            <
                typename GeometricField<Type, fvPatchField, volMesh>::Boundary
            >& sensFieldPtr,
            const word& name
        ) const;

        //- Constructs pointField based on boundaryField and writes it
        template<class Type>
        void constructAndWriteSensitivtyPointField
        (
            const autoPtr<List<Field<Type>>>& sensFieldPtr,
            const word& name
        ) const;

        //- Constructs volField based on boundaryField and writes it
        template<class Type>
        tmp<GeometricField<Type, fvPatchField, volMesh>>
        constructVolSensitivtyField
        (
            const autoPtr
            <
                typename GeometricField<Type, fvPatchField, volMesh>::Boundary
            >& sensFieldPtr,
            const word& name
        ) const;

        //- Write face-based sensitivities, if present
        void writeFaceBasedSens() const;

        //- Write point-based sensitivities, if present
        void writePointBasedSens() const;


private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        sensitivity(const sensitivity&) = delete;

        //- Disallow default bitwise assignment
        void operator=(const sensitivity&) = delete;


public:

    //- Runtime type information
    TypeName("sensitivity");

    // Constructors

        //- Construct from components
        sensitivity
        (
            const fvMesh& mesh,
            const dictionary& dict,
            const word& adjointSolverName
        );

    //- Destructor
    virtual ~sensitivity() = default;


    // Member Functions

        //- Return the construction dictionary
        const dictionary& dict() const;

        //- Read dictionary if changed
        virtual bool readDict(const dictionary& dict);

        //- Get patch IDs on which sensitivities are computed
        const labelList& sensitivityPatchIDs() const;

        //- Overwrite sensitivityPatchIDs
        void setSensitivityPatchIDs(const labelList& sensPatchIDs);

        //- Compute design variables number. Does nothing in the base
        //  Used to get the correct design variables number when
        //  setSensitivityPatchIDs are not set in the constructor
        virtual void computeDerivativesSize();

        //- Calculates and returns sensitivity fields.
        //  Used with optimisation libraries
        virtual const scalarField& calculateSensitivities() = 0;

        //- Write sensitivity fields.
        //  If valid, copies boundaryFields to volFields and writes them.
        //  Virtual to be reimplemented by control points-based methods
        //  (Bezier, RBF) which do not need to write fields
        virtual void write(const word& baseName = word::null);

        //- Get wall face sensitivity vectors field
        tmp<volVectorField> getWallFaceSensVec();

        //- Get wall face sensitivity projected to normal field
        tmp<volScalarField> getWallFaceSensNormal();

        //- Get wall face normal sens as vectors field
        tmp<volVectorField> getWallFaceSensNormalVec();

        //- Get wall point sensitivity vectors field
        //  Uses volPointInterpolation
        tmp<pointVectorField> getWallPointSensVec();

        //- Get wall point sensitivity projected to normal field
        //  Uses volPointInterpolation
        tmp<pointScalarField> getWallPointSensNormal();

        //- Get wall point sens as vectors field
        //  Uses volPointInterpolation
        tmp<pointVectorField> getWallPointSensNormalVec();

        //- Get wall face sensitivity vectors field
        virtual const boundaryVectorField& getWallFaceSensVecBoundary() const;

        //- Get wall face sensitivity projected to normal field
        virtual const boundaryScalarField&
            getWallFaceSensNormalBoundary() const;

        //- Get wall face normal sens as vectors field
        virtual const boundaryVectorField&
            getWallFaceSensNormalVecBoundary() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "sensitivityTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
