/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2017 OpenCFD Ltd.
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

#include "setFlow.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvcFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcSurfaceIntegrate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(setFlow, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        setFlow,
        dictionary
    );
}
}


const Foam::Enum
<
    Foam::functionObjects::setFlow::modeType
>
Foam::functionObjects::setFlow::modeTypeNames
{
    { functionObjects::setFlow::modeType::FUNCTION, "function" },
    { functionObjects::setFlow::modeType::ROTATION, "rotation" },
    { functionObjects::setFlow::modeType::VORTEX2D, "vortex2D" },
    { functionObjects::setFlow::modeType::VORTEX3D, "vortex3D" }
};


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::setFlow::setPhi(const volVectorField& U)
{
    surfaceScalarField* phiptr =
        mesh_.lookupObjectRefPtr<surfaceScalarField>(phiName_);

    if (!phiptr)
    {
        return;
    }

    if (rhoName_ != "none")
    {
        const volScalarField* rhoptr =
            mesh_.lookupObjectPtr<volScalarField>(rhoName_);

        if (rhoptr)
        {
            const volScalarField& rho = *rhoptr;
            *phiptr = fvc::flux(rho*U);
        }
        else
        {
            FatalErrorInFunction
                << "Unable to find rho field'" << rhoName_
                << "' in the mesh database.  Available fields are:"
                << mesh_.names<volScalarField>()
                << exit(FatalError);
        }
    }
    else
    {
        *phiptr = fvc::flux(U);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::setFlow::setFlow
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    UName_("U"),
    rhoName_("none"),
    phiName_("phi"),
    mode_(modeType::FUNCTION),
    reverseTime_(VGREAT),
    scalePtr_(nullptr),
    origin_(vector::zero),
    R_(tensor::I),
    omegaPtr_(nullptr),
    velocityPtr_(nullptr)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::setFlow::~setFlow()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::setFlow::read(const dictionary& dict)
{
    if (fvMeshFunctionObject::read(dict))
    {
        Info<< name() << ":" << endl;
        mode_ = modeTypeNames.read(dict.lookup("mode"));

        Info<< "    operating mode: " << modeTypeNames[mode_] << endl;

        if (dict.readIfPresent("U", UName_))
        {
            Info<< "    U field name: " << UName_ << endl;
        }

        if (dict.readIfPresent("rho", rhoName_))
        {
            Info<< "    rho field name: " << rhoName_ << endl;
        }

        if (dict.readIfPresent("phi", phiName_))
        {
            Info<< "    phi field name: " << phiName_ << endl;
        }

        if (dict.readIfPresent("reverseTime", reverseTime_))
        {
            Info<< "    reverse flow direction at time: " << reverseTime_
                << endl;
            reverseTime_ = mesh_.time().userTimeToTime(reverseTime_);
        }

        // Scaling is applied across all modes
        scalePtr_ = Function1<scalar>::New("scale", dict);

        switch (mode_)
        {
            case modeType::FUNCTION:
            {
                velocityPtr_ = Function1<vector>::New("velocity", dict);
                break;
            }
            case modeType::ROTATION:
            {
                omegaPtr_ = Function1<scalar>::New("omega", dict);
                dict.lookup("origin") >> origin_;
                vector refDir(dict.lookup("refDir"));
                refDir /= mag(refDir) + ROOTVSMALL;
                vector axis(dict.lookup("axis"));
                axis /= mag(axis) + ROOTVSMALL;
                R_ = tensor(refDir, axis, refDir^axis);
                break;
            }
            case modeType::VORTEX2D:
            case modeType::VORTEX3D:
            {
                dict.lookup("origin") >> origin_;
                vector refDir(dict.lookup("refDir"));
                refDir /= mag(refDir) + ROOTVSMALL;
                vector axis(dict.lookup("axis"));
                axis /= mag(axis) + ROOTVSMALL;
                R_ = tensor(refDir, axis, refDir^axis);
                break;
            }
        }

        Info<< endl;

        return true;
    }

    return false;
}


bool Foam::functionObjects::setFlow::execute()
{
    volVectorField* Uptr = mesh_.lookupObjectRefPtr<volVectorField>(UName_);

    surfaceScalarField* phiptr =
        mesh_.lookupObjectRefPtr<surfaceScalarField>(phiName_);

    Log << nl << name() << ":" << nl;

    if (!Uptr || !phiptr)
    {
        Info<< "    Either field " << UName_ << " or " << phiName_
            << " not found in the mesh database" << nl;

        return true;
    }

    const scalar t = mesh_.time().timeOutputValue();

    Log << "    setting " << UName_ << " and " << phiName_ << nl;

    volVectorField& U = *Uptr;
    surfaceScalarField& phi = *phiptr;

    switch (mode_)
    {
        case modeType::FUNCTION:
        {
            const vector Uc = velocityPtr_->value(t);
            U == dimensionedVector("Uc", dimVelocity, Uc);
            U.correctBoundaryConditions();
            setPhi(U);

            break;
        }
        case modeType::ROTATION:
        {
            const volVectorField& C = mesh_.C();
            const volVectorField d
            (
                typeName + ":d",
                C - dimensionedVector("origin", dimLength, origin_)
            );
            const scalarField x(d.component(vector::X));
            const scalarField z(d.component(vector::Z));

            const scalar omega = omegaPtr_->value(t);
            vectorField& Uc = U.primitiveFieldRef();
            Uc.replace(vector::X, -omega*z);
            Uc.replace(vector::Y, scalar(0));
            Uc.replace(vector::Z, omega*x);

            volVectorField::Boundary& Ubf = U.boundaryFieldRef();
            forAll(Ubf, patchi)
            {
                vectorField& Uf = Ubf[patchi];
                if (Uf.size())
                {
                    const vectorField& Cp = C.boundaryField()[patchi];
                    const vectorField dp(Cp - origin_);
                    const scalarField xp(dp.component(vector::X));
                    const scalarField zp(dp.component(vector::Z));
                    Uf.replace(vector::X, -omega*zp);
                    Uf.replace(vector::Y, scalar(0));
                    Uf.replace(vector::Z, omega*xp);
                }
            }

            U = U & R_;
            U.correctBoundaryConditions();
            setPhi(U);

            break;
        }
        case modeType::VORTEX2D:
        {
            const scalar pi = Foam::constant::mathematical::pi;

            const volVectorField& C = mesh_.C();

            const volVectorField d
            (
                typeName + ":d",
                C - dimensionedVector("origin", dimLength, origin_)
            );
            const scalarField x(d.component(vector::X));
            const scalarField z(d.component(vector::Z));

            vectorField& Uc = U.primitiveFieldRef();
            Uc.replace(vector::X, -sin(2*pi*z)*sqr(sin(pi*x)));
            Uc.replace(vector::Y, scalar(0));
            Uc.replace(vector::Z, sin(2*pi*x)*sqr(sin(pi*z)));

            U = U & R_;

            // Calculating incompressible flux based on stream function
            // Note: R_ rotation not implemented in phi calculation
            const scalarField xp(mesh_.points().component(0) - origin_[0]);
            const scalarField yp(mesh_.points().component(1) - origin_[1]);
            const scalarField zp(mesh_.points().component(2) - origin_[2]);
            const scalarField psi((1.0/pi)*sqr(sin(pi*xp))*sqr(sin(pi*zp)));

            scalarField& phic = phi.primitiveFieldRef();
            forAll(phic, fi)
            {
                phic[fi] = 0;
                const face& f = mesh_.faces()[fi];
                const label nPoints = f.size();

                forAll(f, fpi)
                {
                    const label p1 = f[fpi];
                    const label p2 = f[(fpi + 1) % nPoints];
                    phic[fi] += 0.5*(psi[p1] + psi[p2])*(yp[p2] - yp[p1]);
                }
            }

            surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();
            forAll(phibf, patchi)
            {
                scalarField& phif = phibf[patchi];
                const label start = mesh_.boundaryMesh()[patchi].start();

                forAll(phif, fi)
                {
                    phif[fi] = 0;
                    const face& f = mesh_.faces()[start + fi];
                    const label nPoints = f.size();

                    forAll(f, fpi)
                    {
                        const label p1 = f[fpi];
                        const label p2 = f[(fpi + 1) % nPoints];
                        phif[fi] += 0.5*(psi[p1] + psi[p2])*(yp[p2] - yp[p1]);
                    }
                }
            }

            break;
        }
        case modeType::VORTEX3D:
        {
            const scalar pi = Foam::constant::mathematical::pi;
            const volVectorField& C = mesh_.C();

            const volVectorField d
            (
                typeName + ":d",
                C - dimensionedVector("origin", dimLength, origin_)
            );
            const scalarField x(d.component(vector::X));
            const scalarField y(d.component(vector::Y));
            const scalarField z(d.component(vector::Z));

            vectorField& Uc = U.primitiveFieldRef();
            Uc.replace(vector::X, 2*sqr(sin(pi*x))*sin(2*pi*y)*sin(2*pi*z));
            Uc.replace(vector::Y, -sin(2*pi*x)*sqr(sin(pi*y))*sin(2*pi*z));
            Uc.replace(vector::Z, -sin(2*pi*x)*sin(2*pi*y)*sqr(sin(pi*z)));

            U = U & R_;
            U.correctBoundaryConditions();

            // Calculating phi
            // Note: R_ rotation not implemented in phi calculation
            const vectorField Cf(mesh_.Cf().primitiveField() - origin_);
            const scalarField Xf(Cf.component(vector::X));
            const scalarField Yf(Cf.component(vector::Y));
            const scalarField Zf(Cf.component(vector::Z));
            vectorField Uf(Xf.size());
            Uf.replace(0, 2*sqr(sin(pi*Xf))*sin(2*pi*Yf)*sin(2*pi*Zf));
            Uf.replace(1, -sin(2*pi*Xf)*sqr(sin(pi*Yf))*sin(2*pi*Zf));
            Uf.replace(2, -sin(2*pi*Xf)*sin(2*pi*Yf)*sqr(sin(pi*Zf)));

            scalarField& phic = phi.primitiveFieldRef();
            const vectorField& Sfc = mesh_.Sf().primitiveField();
            phic = Uf & Sfc;

            surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();
            const surfaceVectorField::Boundary& Sfbf =
                mesh_.Sf().boundaryField();
            const surfaceVectorField::Boundary& Cfbf =
                mesh_.Cf().boundaryField();

            forAll(phibf, patchi)
            {
                scalarField& phif = phibf[patchi];
                const vectorField& Sff = Sfbf[patchi];
                const vectorField& Cff = Cfbf[patchi];
                const scalarField xf(Cff.component(vector::X));
                const scalarField yf(Cff.component(vector::Y));
                const scalarField zf(Cff.component(vector::Z));
                vectorField Uf(xf.size());
                Uf.replace(0, 2*sqr(sin(pi*xf))*sin(2*pi*yf)*sin(2*pi*zf));
                Uf.replace(1, -sin(2*pi*xf)*sqr(sin(pi*yf))*sin(2*pi*zf));
                Uf.replace(2, -sin(2*pi*xf)*sin(2*pi*yf)*sqr(sin(pi*zf)));
                phif = Uf & Sff;
            }

            break;
        }
    }

    if (t > reverseTime_)
    {
        Log << "    flow direction: reverse" << nl;
        U.negate();
        phi.negate();
    }

    // Apply scaling
    const scalar s = scalePtr_->value(t);
    U *= s;
    phi *= s;

    U.correctBoundaryConditions();

    const scalarField sumPhi(fvc::surfaceIntegrate(phi));
    Log << "    Continuity error: max(mag(sum(phi))) = "
        << gMax(mag(sumPhi)) << nl << endl;

    return true;
}


bool Foam::functionObjects::setFlow::write()
{
    const auto& Uptr = mesh_.lookupObjectRefPtr<volVectorField>(UName_);
    if (Uptr)
    {
        Uptr->write();
    }

    const auto& phiptr = mesh_.lookupObjectRefPtr<surfaceScalarField>(phiName_);
    if (phiptr)
    {
        phiptr->write();
    }

    return true;
}


// ************************************************************************* //
