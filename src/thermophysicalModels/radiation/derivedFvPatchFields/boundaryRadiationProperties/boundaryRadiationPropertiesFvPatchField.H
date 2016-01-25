/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd
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
    Foam::boundaryRadiationPropertiesFvPatchField

Description
    Common functions to emissivity. It gets supplied from lookup into a
    dictionary or calculated by the solidThermo:

    - 'lookup' : lookup volScalarField with name
    - 'solidThermo' : use solidThermo
    - 'model' : use a model

SourceFiles
    boundaryRadiationPropertiesFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRadiationPropertiesFvPatchField_H
#define boundaryRadiationPropertiesFvPatchField_H

#include "scalarField.H"
#include "NamedEnum.H"
#include "fvPatch.H"
#include "calculatedFvPatchFields.H"
#include "transmissivityModel.H"
#include "absorptionEmissionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
/*---------------------------------------------------------------------------*\
             Class boundaryRadiationPropertiesFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class boundaryRadiationPropertiesFvPatchField
:
    public calculatedFvPatchScalarField
{
public:

    // - Type of look up
    enum methodType
    {
        SOLIDRADIATION,
        LOOKUP,
        MODEL
    };


private:

    // Private data

        static const NamedEnum<methodType, 3> methodTypeNames_;

        //- How to get property
        const methodType method_;

        //- Dictionary
        const dictionary dict_;

        //- Absorption/emission model
        autoPtr<absorptionEmissionModel> absorptionEmission_;

        //- transmissivity model
        autoPtr<transmissivityModel> transmissivity_;


    // Private functions

        //- Return nbr patch index
        label nbrPatchIndex() const;

        //- Return nbr mesh
        const fvMesh& nbrRegion() const;


public:

    //- Runtime type information
    TypeName("boundaryRadiation");

    // Constructors


        //- Construct from patch and internal field
        boundaryRadiationPropertiesFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        boundaryRadiationPropertiesFvPatchField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given a
        // greyDiffusiveViewFactorFixedValueFvPatchScalarField onto a new patch
        boundaryRadiationPropertiesFvPatchField
        (
            const boundaryRadiationPropertiesFvPatchField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        boundaryRadiationPropertiesFvPatchField
        (
            const boundaryRadiationPropertiesFvPatchField&
        );

        //- Construct as copy setting internal field reference
        boundaryRadiationPropertiesFvPatchField
        (
            const boundaryRadiationPropertiesFvPatchField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new boundaryRadiationPropertiesFvPatchField(*this)
            );
        }

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new boundaryRadiationPropertiesFvPatchField(*this, iF)
            );
        }


    // Member functions

        //- Method to obtain properties
        word method() const
        {
            return methodTypeNames_[method_];
        }


        //- Return absorptionEmissionModel
        const absorptionEmissionModel& absorptionEmission() const;

        //- Return transmissivityModel
        const transmissivityModel& transmissiveModel() const;

        //- Calculate corresponding emissivity field for bandI
        tmp<scalarField> emissivity(const label bandI = 0) const;

        //- Calculate corresponding absorptivity field for bandI
        tmp<scalarField> absorptivity(const label bandI = 0) const;

        //- Calculate corresponding transmissivity field for bandI
        tmp<scalarField> transmissivity(const label bandI = 0) const;

        //- Calculate corresponding reflectivity field
        tmp<scalarField> reflectivity(const label bandI = 0) const;

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
