/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2011 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "radiationCoupledBase.H"
#include "volFields.H"
#include "basicSolidThermo.H"

#include "directMappedPatchBase.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

namespace Foam
{
    template<>
    const char* Foam::NamedEnum
    <
        Foam::radiationCoupledBase::emissivityMethodType,
        2
    >::names[] =
    {
        "solidThermo",
        "lookup"
    };
}


const Foam::NamedEnum<Foam::radiationCoupledBase::emissivityMethodType, 2>
    Foam::radiationCoupledBase::emissivityMethodTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const word& calculationType
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_[calculationType]),
    emissivity_(patch.size(), 0.0)
{}


Foam::radiationCoupledBase::radiationCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(emissivityMethodTypeNames_.read(dict.lookup("emissivityMode")))
{
    switch (method_)
    {
        case SOLIDTHERMO:
        {
            if (!isA<directMappedPatchBase>(patch_.patch()))
            {
                FatalErrorIn
                (
                    "radiationCoupledBase::radiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "\n    patch type '" << patch_.type()
                    << "' not type '" << directMappedPatchBase::typeName << "'"
                    << "\n    for patch " << patch_.name()
                    << exit(FatalError);
            }

            const directMappedPatchBase& mpp = refCast
            <
                const directMappedPatchBase
            >
            (
                patch_.patch()
            );

            const polyMesh& nbrMesh = mpp.sampleMesh();

            if
            (
                !nbrMesh.foundObject<basicSolidThermo>
                (
                    "solidThermophysicalProperties"
                )
            )
            {
                FatalErrorIn
                (
                    "radiationCoupledBase::radiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "\n solidThermophysicalProperties does not exist "
                    << "\n in mesh ' " << nbrMesh.name() << "'"
                    << exit(FatalError);
            }
        }
        break;

        case LOOKUP:
        {
            if(!dict.found("emissivity"))
            {
                FatalErrorIn
                (
                    "radiationCoupledBase::radiationCoupledBase\n"
                    "(\n"
                    "    const fvPatch& p,\n"
                    "    const dictionary& dict\n"
                    ")\n"
                )   << "\n    emissivity key"
                    << "\n    does not exist "
                    << "\n    for patch " << patch_.name()
                    << exit(FatalError);
            }
            else
            {
                emissivity_ = scalarField("emissivity", dict, patch_.size());
            }
        }
        break;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::radiationCoupledBase::emissivity() const
{
    switch (method_)
    {
        case SOLIDTHERMO:
        {

            // Get the coupling information from the directMappedPatchBase
            const directMappedPatchBase& mpp =
                refCast<const directMappedPatchBase>
                (
                    patch_.patch()
                );

            const polyMesh& nbrMesh = mpp.sampleMesh();
            const fvPatch& nbrPatch = refCast<const fvMesh>
            (
                nbrMesh
            ).boundary()[mpp.samplePolyPatch().index()];

            scalarField emissivity
            (
                nbrPatch.lookupPatchField<volScalarField, scalar>("emissivity")
            );

            // Force recalculation of mapping and schedule
            const mapDistribute& distMap = mpp.map();

            distMap.distribute(emissivity);

            return tmp<scalarField>
            (
                new scalarField(emissivity)
            );

        }
        break;

        case LOOKUP:
        {
            // return local value
            return emissivity_;
        }

        default:
        {
            FatalErrorIn
            (
                "radiationCoupledBase::emissivity(const scalarField&)"
            )
                << "Unimplemented method " << method_ << endl
                << "Please set 'emissivity' to one of "
                << emissivityMethodTypeNames_.toc()
                << " and 'emissivityName' to the name of the volScalar"
                << exit(FatalError);
        }
        break;
    }
    return scalarField(0);
}


void Foam::radiationCoupledBase::write(Ostream& os) const
{
    os.writeKeyword("emissivityMode") << emissivityMethodTypeNames_[method_]
        << token::END_STATEMENT << nl;
    os.writeKeyword("emissivity") << emissivity_
        << token::END_STATEMENT << nl;
}


// ************************************************************************* //
