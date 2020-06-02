/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017-2020 OpenCFD Ltd.
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
{
    { KMethodType::mtFluidThermo, "fluidThermo" },
    { KMethodType::mtSolidThermo, "solidThermo" },
    { KMethodType::mtDirectionalSolidThermo, "directionalSolidThermo" },
    { KMethodType::mtLookup, "lookup" }
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const word& calculationType,
    const word& kappaName,
    const word& alphaAniName,
    const word& alphaName
)
:
    patch_(patch),
    method_(KMethodTypeNames_[calculationType]),
    kappaName_(kappaName),
    alphaAniName_(alphaAniName),
    alphaName_(alphaName)
{}


Foam::temperatureCoupledBase::temperatureCoupledBase
(
    const fvPatch& patch,
    const dictionary& dict
)
:
    patch_(patch),
    method_(KMethodTypeNames_.get("kappaMethod", dict)),
    kappaName_(dict.getOrDefault<word>("kappa", "none")),
    alphaAniName_(dict.getOrDefault<word>("alphaAni", "none")),
    alphaName_(dict.getOrDefault<word>("alpha", "none"))
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
    alphaAniName_(base.alphaAniName_),
    alphaName_(base.alphaName_)
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

            {
                const auto* ptr =
                    mesh.cfindObject<turbulenceModel>
                    (
                        turbulenceModel::propertiesName
                    );

                if (ptr)
                {
                    return ptr->kappaEff(patchi);
                }
            }

            {
                const auto* ptr =
                    mesh.cfindObject<fluidThermo>(basicThermo::dictName);

                if (ptr)
                {
                    return ptr->kappa(patchi);
                }
            }

            {
                const auto* ptr =
                    mesh.cfindObject<basicThermo>(basicThermo::dictName);

                if (ptr)
                {
                    return ptr->kappa(patchi);
                }
            }

            {
                const auto* ptr =
                    mesh.cfindObject<basicThermo>("phaseProperties");

                if (ptr)
                {
                    return ptr->kappa(patchi);
                }
            }

            FatalErrorInFunction
                << "Using kappaMethod " << KMethodTypeNames_[method_]
                << ", but thermo package not available\n"
                << exit(FatalError);

            break;
        }

        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

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

            break;
        }
    }

    return scalarField();
}


Foam::tmp<Foam::scalarField> Foam::temperatureCoupledBase::alpha
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

            {
                const auto* ptr =
                    mesh.cfindObject<turbulenceModel>
                    (
                        turbulenceModel::propertiesName
                    );

                if (ptr)
                {
                    return ptr->alphaEff(patchi);
                }
            }

            {
                const auto* ptr =
                    mesh.cfindObject<fluidThermo>(basicThermo::dictName);

                if (ptr)
                {
                    return ptr->alpha(patchi);
                }
            }

            {
                const auto* ptr =
                    mesh.cfindObject<basicThermo>(basicThermo::dictName);

                if (ptr)
                {
                    return ptr->alpha(patchi);
                }
            }

            {
                const auto* ptr =
                    mesh.cfindObject<basicThermo>("phaseProperties");

                if (ptr)
                {
                    return ptr->alpha(patchi);
                }
            }

            FatalErrorInFunction
                << "Using kappaMethod " << KMethodTypeNames_[method_]
                << ", but thermo package not available\n"
                << exit(FatalError);

            break;
        }

        case mtSolidThermo:
        {
            const solidThermo& thermo =
                mesh.lookupObject<solidThermo>(basicThermo::dictName);

            return thermo.alpha(patchi);
            break;
        }

        case mtDirectionalSolidThermo:
        {
            const symmTensorField& alphaAni =
                patch_.lookupPatchField<volSymmTensorField, scalar>
                (
                    alphaAniName_
                );

            const vectorField n(patch_.nf());

            return n & alphaAni & n;
        }

        case mtLookup:
        {
            if (mesh.foundObject<volScalarField>(alphaName_))
            {
                return
                    patch_.lookupPatchField<volScalarField, scalar>
                    (
                        alphaName_
                    );
            }
            else if (mesh.foundObject<volSymmTensorField>(alphaName_))
            {
                const symmTensorField& alphaWall =
                    patch_.lookupPatchField<volSymmTensorField, scalar>
                    (
                        alphaName_
                    );

                const vectorField n(patch_.nf());

                return n & alphaWall & n;
            }
            else
            {
                FatalErrorInFunction
                    << "Did not find field " << alphaName_
                    << " on mesh " << mesh.name() << " patch " << patch_.name()
                    << nl
                    << "Please set 'kappaMethod' to one of "
                    << flatOutput(KMethodTypeNames_.sortedToc()) << nl
                    << "and 'alpha' to the name of the volScalar"
                    << " or volSymmTensor field (if kappaMethod=lookup)"
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
                << "and 'alpha' to the name of the volScalar"
                << " or volSymmTensor field (if kappaMethod=lookup)"
                << exit(FatalError);

            break;
        }
    }

    return scalarField();
}


void Foam::temperatureCoupledBase::write(Ostream& os) const
{
    os.writeEntry("kappaMethod", KMethodTypeNames_[method_]);
    os.writeEntry("kappa", kappaName_);
    os.writeEntry("alphaAni", alphaAniName_);
    os.writeEntry("alpha", alphaName_);
}


// ************************************************************************* //
