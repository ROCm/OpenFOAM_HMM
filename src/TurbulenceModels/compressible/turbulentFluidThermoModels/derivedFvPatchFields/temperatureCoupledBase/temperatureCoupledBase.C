/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017-2018 OpenCFD Ltd.
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

#include "temperatureCoupledBase.H"
#include "volFields.H"
#include "fluidThermo.H"
#include "solidThermo.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * * * //

const Foam::Enum
<
    Foam::temperatureCoupledBase::KMethodType
>
Foam::temperatureCoupledBase::KMethodTypeNames_
({
    { KMethodType::mtFluidThermo, "fluidThermo" },
    { KMethodType::mtSolidThermo, "solidThermo" },
    { KMethodType::mtDirectionalSolidThermo, "directionalSolidThermo" },
    { KMethodType::mtLookup, "lookup" },
});


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName),
    alphaAniName_(alphaAniName)
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.get("kappaMethod", dict)),
    kappaName_(dict.lookupOrDefault<word>("kappa", "none")),
    alphaAniName_(dict.lookupOrDefault<word>("alphaAni","none"))
{
    switch (method_)
    {
        case mtDirectionalSolidThermo:
        {
            if (!dict.found("alphaAni"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'alphaAni'"
                       " required for 'kappaMethod' "
                    << KMethodTypeNames_[method_]
                    << exit(FatalIOError);
            }

            break;
        }

        case mtLookup:
        {
            if (!dict.found("kappa"))
            {
                FatalIOErrorInFunction(dict)
                    << "Did not find entry 'kappa'"
                       " required for 'kappaMethod' "
                    <<  KMethodTypeNames_[method_] << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                       " or volSymmTensorField"
                    << exit(FatalIOError);
            }

            break;
        }

        default:
        {
            break;
        }
    }
}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const temperatureCoupledBase& base
)
:
    patch_(patch),
    method_(base.method_),
    kappaName_(base.kappaName_),
    alphaAniName_(base.alphaAniName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase::kappa
(
    const scalarField& Tp
) const
{
    const fvMesh& mesh = patch_.boundaryMesh().mesh();
    const label patchi = patch_.index();

    switch (method_)
    {
        case mtFluidThermo:
        {
            typedef compressible::turbulenceModel turbulenceModel;

            const word turbName(turbulenceModel::propertiesName);

            if
            (
                mesh.foundObject<turbulenceModel>(turbName)
            )
            {
                const turbulenceModel& turbModel =
                    mesh.lookupObject<turbulenceModel>(turbName);

                return turbModel.kappaEff(patchi);
            }
            else if (mesh.foundObject<fluidThermo>(basicThermo::dictName))
            {
                const fluidThermo& thermo =
                    mesh.lookupObject<fluidThermo>(basicThermo::dictName);

                return thermo.kappa(patchi);
            }
            else if (mesh.foundObject<basicThermo>(basicThermo::dictName))
            {
                const basicThermo& thermo =
                    mesh.lookupObject<basicThermo>(basicThermo::dictName);

                return thermo.kappa(patchi);
            }
            else if (mesh.foundObject<basicThermo>("phaseProperties"))
            {
                const basicThermo& thermo =
                    mesh.lookupObject<basicThermo>("phaseProperties");

                return thermo.kappa(patchi);
            }
            else
            {
                FatalErrorInFunction
                    << "kappaMethod defined to employ "
                    << KMethodTypeNames_[method_]
                    << " method, but thermo package not available"
                    << exit(FatalError);
            }

            break;
        }

        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            if (!thermo.isotropic())
            {
                word regionName = "";
                if (mesh.name() != polyMesh::defaultRegion)
                {
                    regionName = " for region " + mesh.name();
                }

                const word& patchName = mesh.boundaryMesh()[patchi].name();

                WarningInFunction
                    << "Applying isotropic thermal conductivity assumption to "
                    << "anisotropic model" << regionName << " at patch "
                    << patchName << nl
                    << "Consider using an isotropic conductivity model or "
                    << "set 'kappaMethod' to "
                    << KMethodTypeNames_[mtDirectionalSolidThermo]
                    << nl << endl;
            }

            return thermo.kappa(patchi);
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            const symmTensorField& alphaAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    alphaAniName_
                );

            const scalarField& pp = thermo.p().boundaryField()[patchi];

            const symmTensorField kappa(alphaAni*thermo.Cp(pp, Tp, patchi));

            const vectorField n(patch_.nf());

            return n & kappa & n;
        }

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(kappaName_))
            {
                return patch_.lookupPatchField<volScalarField, scalar>
                (
                    kappaName_
                );
            }
            else if (mesh.foundObject<volSymmTensorField>(kappaName_))
            {
                const symmTensorField& KWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        kappaName_
                    );

                const vectorField n(patch_.nf());

                return n & KWall & n;
            }
            else
            {
                FatalErrorInFunction
                    << "Did not find field " << kappaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "    Please set 'kappa' to the name of a volScalarField"
                    << " or volSymmTensorField."
                    << exit(FatalError);
            }

            break;
        }

        default:
        {
            FatalErrorInFunction
                << "Unimplemented method " << KMethodTypeNames_[method_] << nl
                << "Please set 'kappaMethod' to one of "
                << flatOutput(KMethodTypeNames_.sortedToc()) << nl
                << "and 'kappa' to the name of the volScalar"
                << " or volSymmTensor field (if kappaMethod=lookup)"
                << exit(FatalError);
        }
    }

    return scalarField(0);
}


void Foam::temperatureCoupledBase::write(Ostream& os) const
{
    os.writeEntry("kappaMethod", KMethodTypeNames_[method_]);
    os.writeEntry("kappa", kappaName_);
    os.writeEntry("alphaAni", alphaAniName_);
}


// ************************************************************************* //
