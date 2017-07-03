/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::boundaryRadiationPropertiesPatch

Description
    Common functions to emissivity. It gets supplied from lookup into a
    dictionary or calculated by the solidThermo:

    - 'lookup' : lookup volScalarField with name
    - 'solidThermo' : use solidThermo
    - 'model' : use a model

SourceFiles
    boundaryRadiationPropertiesPatch.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryRadiationPropertiesPatch_H
#define boundaryRadiationPropertiesPatch_H

#include "scalarField.H"
#include "Enum.H"
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
             Class boundaryRadiationPropertiesPatch Declaration
\*---------------------------------------------------------------------------*/

class boundaryRadiationPropertiesPatch
{
public:

    // - Type of method
    enum methodType
    {
        SOLIDRADIATION,
        LOOKUP,
        MODEL
    };


private:

    // Private data

        static const Enum<methodType> methodTypeNames_;

        //- How to get property
        const methodType method_;

        //- Dictionary
        const dictionary dict_;

        //- Absorption/emission model
        autoPtr<absorptionEmissionModel> absorptionEmission_;

        //- transmissivity model
        autoPtr<transmissivityModel> transmissivity_;

        //- reference to fvPatch
        const polyPatch& patch_;


    // Private functions

        //- Return nbr patch index
        label nbrPatchIndex() const;

        //- Return nbr mesh
        const fvMesh& nbrRegion() const;


public:


    // Constructors


        //- Construct from patch, internal field and dictionary
        boundaryRadiationPropertiesPatch
        (
            const polyPatch&,
            const dictionary&
        );

         //- Destructor
        virtual ~boundaryRadiationPropertiesPatch();


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
        void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam
} // End namespace radiation

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
