/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
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

#include "directionalSolidThermo.H"
#include "addToRunTimeSelectionTable.H"
#include "interpolateXY.H"
#include "transform.H"
#include "transformField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(directionalSolidThermo, 0);
    addToRunTimeSelectionTable
    (
        basicSolidThermo,
        directionalSolidThermo,
        mesh
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::directionalSolidThermo::directionalSolidThermo(const fvMesh& mesh)
:
    basicSolidThermo(mesh),
    ccTransforms_
    (
        IOobject
        (
            "ccTransforms",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimLength
    )
{
    read();

    // Determine transforms for cell centres
    forAll(mesh.C(), cellI)
    {
        vector dir = mesh.C()[cellI] - coordSys_.origin();
        dir /= mag(dir);

        // Define local coordinate system with
        // - e1 : axis from cc to centre
        // - e3 : rotation axis
        coordinateSystem cs
        (
            "cc",
            coordSys_.origin(),
            coordSys_.e3(),     //z',e3
            dir                 //x',e1
        );

        ccTransforms_[cellI] = cs.R();
    }

    forAll(mesh.C().boundaryField(), patchI)
    {
        const fvPatchVectorField& patchC = mesh.C().boundaryField()[patchI];
        fvPatchTensorField& patchT = ccTransforms_.boundaryField()[patchI];

        tensorField tc(patchT.size());
        forAll(tc, i)
        {
            vector dir = patchC[i] - coordSys_.origin();
            dir /= mag(dir);

            coordinateSystem cs
            (
                "cc",
                coordSys_.origin(),
                coordSys_.e3(),     //z',e3
                dir                 //x',e1
            );

            tc[i] = cs.R();
        }
        patchT = tc;
    }

    if (debug)
    {
        Info<< "directionalSolidThermo : dumping converted Kxx, Kyy, Kzz"
            << endl;
        {
            volVectorField Kxx
            (
                IOobject
                (
                    "Kxx",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimless
            );
            Kxx.internalField() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(1, 0, 0)
                )
            );
            forAll(Kxx.boundaryField(), patchI)
            {
                Kxx.boundaryField()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(1, 0, 0)
                    )
                );
            }
            Kxx.write();
        }
        {
            volVectorField Kyy
            (
                IOobject
                (
                    "Kyy",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimless
            );
            Kyy.internalField() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(0, 1, 0)
                )
            );
            forAll(Kyy.boundaryField(), patchI)
            {
                Kyy.boundaryField()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(0, 1, 0)
                    )
                );
            }
            Kyy.write();
        }
        {
            volVectorField Kzz
            (
                IOobject
                (
                    "Kzz",
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    false
                ),
                mesh,
                dimless
            );
            Kzz.internalField() = transform
            (
                ccTransforms_.internalField(),
                vectorField
                (
                    ccTransforms_.internalField().size(),
                    point(0, 0, 1)
                )
            );
            forAll(Kzz.boundaryField(), patchI)
            {
                Kzz.boundaryField()[patchI] = transform
                (
                    ccTransforms_.boundaryField()[patchI],
                    vectorField
                    (
                        ccTransforms_.boundaryField()[patchI].size(),
                        point(0, 0, 1)
                    )
                );
            }
            Kzz.write();
        }
    }



    correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::directionalSolidThermo::~directionalSolidThermo()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::symmTensor Foam::directionalSolidThermo::transformPrincipal
(
    const tensor& tt,
    const vector& st
) const
{
    return symmTensor
    (
        tt.xx()*st.x()*tt.xx()
      + tt.xy()*st.y()*tt.xy()
      + tt.xz()*st.z()*tt.xz(),

        tt.xx()*st.x()*tt.yx()
      + tt.xy()*st.y()*tt.yy()
      + tt.xz()*st.z()*tt.yz(),

        tt.xx()*st.x()*tt.zx()
      + tt.xy()*st.y()*tt.zy()
      + tt.xz()*st.z()*tt.zz(),

        tt.yx()*st.x()*tt.yx()
      + tt.yy()*st.y()*tt.yy()
      + tt.yz()*st.z()*tt.yz(),

        tt.yx()*st.x()*tt.zx()
      + tt.yy()*st.y()*tt.zy()
      + tt.yz()*st.z()*tt.zz(),
        
        tt.zx()*st.x()*tt.zx()
      + tt.zy()*st.y()*tt.zy()
      + tt.zz()*st.z()*tt.zz()
    );
}


void Foam::directionalSolidThermo::transformField
(
    symmTensorField& fld,
    const tensorField& tt,
    const vectorField& st
) const
{
    fld.setSize(tt.size());
    forAll(fld, i)
    {
        fld[i] = transformPrincipal(tt[i], st[i]);
    }
}


void Foam::directionalSolidThermo::correct()
{}


Foam::tmp<Foam::volScalarField> Foam::directionalSolidThermo::rho() const
{
    tmp<volScalarField> trho
    (
        new volScalarField
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimDensity
        )
    );
    volScalarField& rho = trho();

    rho.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        rhoValues_
    );

    forAll(rho.boundaryField(), patchI)
    {
        rho.boundaryField()[patchI] == this->rho(patchI)();
    }

    return trho;
}


Foam::tmp<Foam::volScalarField> Foam::directionalSolidThermo::cp() const
{
    tmp<volScalarField> tcp
    (
        new volScalarField
        (
            IOobject
            (
                "cp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/(dimMass*dimTemperature)
        )
    );
    volScalarField& cp = tcp();

    cp.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        cpValues_
    );

    forAll(cp.boundaryField(), patchI)
    {
        cp.boundaryField()[patchI] == this->cp(patchI)();
    }

    return tcp;
}


Foam::tmp<Foam::volSymmTensorField> Foam::directionalSolidThermo::directionalK()
const
{
    tmp<volSymmTensorField> tK
    (
        new volSymmTensorField
        (
            IOobject
            (
                "K",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimTime/(dimLength*dimTemperature)
        )
    );
    volSymmTensorField& K = tK();

    // Get temperature interpolated properties (principal directions)
    Field<vector> localK
    (
        interpolateXY
        (
            T_.internalField(),
            TValues_,
            KValues_
        )
    );

    // Transform into global coordinate system
    transformField(K.internalField(), ccTransforms_.internalField(), localK);

    forAll(K.boundaryField(), patchI)
    {
        K.boundaryField()[patchI] == this->directionalK(patchI)();
    }

    return tK;
}


Foam::tmp<Foam::volScalarField> Foam::directionalSolidThermo::K() const
{
    forAll(KValues_, i)
    {
        const vector& v = KValues_[i];
        if
        (
            v.x() != v.y() 
         || v.x() != v.z() 
         || v.y() != v.z() 
        )
        {
            FatalErrorIn("directionalSolidThermo::K() const")
                << "Supplied K values " << KValues_
                << " are not isotropic." << exit(FatalError);
        }
    }

    // Get temperature interpolated properties (principal directions)
    Field<vector> localK
    (
        interpolateXY
        (
            T_.internalField(),
            TValues_,
            KValues_
        )
    );

    tmp<volScalarField> tK
    (
        new volScalarField
        (
            IOobject
            (
                "K",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimTime/(dimLength*dimTemperature)
        )
    );
    volScalarField& K = tK();

    K.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        KValues_.component(0)()
    );

    forAll(K.boundaryField(), patchI)
    {
        K.boundaryField()[patchI] == this->K(patchI)();
    }

    return tK;
}


Foam::tmp<Foam::volScalarField> Foam::directionalSolidThermo::Hf() const
{
    tmp<volScalarField> tHf
    (
        new volScalarField
        (
            IOobject
            (
                "Hf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimEnergy/dimMass
        )
    );
    volScalarField& Hf = tHf();

    Hf.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        HfValues_
    );

    forAll(Hf.boundaryField(), patchI)
    {
        Hf.boundaryField()[patchI] == this->Hf(patchI)();
    }

    return tHf;
}


Foam::tmp<Foam::volScalarField> Foam::directionalSolidThermo::emissivity() const
{
    tmp<volScalarField> temissivity
    (
        new volScalarField
        (
            IOobject
            (
                "emissivity",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimless
        )
    );
    volScalarField& emissivity = temissivity();

    emissivity.internalField() = interpolateXY
    (
        T_.internalField(),
        TValues_,
        emissivityValues_
    );

    forAll(emissivity.boundaryField(), patchI)
    {
        emissivity.boundaryField()[patchI] == this->emissivity(patchI)();
    }

    return temissivity;
}


Foam::tmp<Foam::scalarField> Foam::directionalSolidThermo::rho
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                rhoValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::directionalSolidThermo::cp
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                cpValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::directionalSolidThermo::K
(
    const label patchI
) const
{
    forAll(KValues_, i)
    {
        const vector& v = KValues_[i];
        if
        (
            v.x() != v.y() 
         || v.x() != v.z() 
         || v.y() != v.z() 
        )
        {
            FatalErrorIn("directionalSolidThermo::K() const")
                << "Supplied K values " << KValues_
                << " are not isotropic." << exit(FatalError);
        }
    }

    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                KValues_.component(0)()
            )
        )
    );
}


Foam::tmp<Foam::symmTensorField> Foam::directionalSolidThermo::directionalK
(
    const label patchI
) const
{
    const fvPatchScalarField& patchT = T_.boundaryField()[patchI];

    Field<vector> localK(interpolateXY(patchT, TValues_, KValues_));

    tmp<symmTensorField> tglobalK(new symmTensorField(localK.size()));
    transformField(tglobalK(), ccTransforms_.boundaryField()[patchI], localK);
    
    return tglobalK;
}


Foam::tmp<Foam::scalarField> Foam::directionalSolidThermo::Hf
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                HfValues_
            )
        )
    );
}


Foam::tmp<Foam::scalarField> Foam::directionalSolidThermo::emissivity
(
    const label patchI
) const
{
    return tmp<scalarField>
    (
        new scalarField
        (
            interpolateXY
            (
                T_.boundaryField()[patchI],
                TValues_,
                emissivityValues_
            )
        )
    );
}


bool Foam::directionalSolidThermo::read()
{
    return read(subDict(typeName + "Coeffs"));
}


bool Foam::directionalSolidThermo::read(const dictionary& dict)
{
    TValues_ = Field<scalar>(dict.lookup("TValues"));
    rhoValues_ = Field<scalar>(dict.lookup("rhoValues"));
    cpValues_ = Field<scalar>(dict.lookup("cpValues"));
    KValues_ = Field<vector>(dict.lookup("KValues"));
    HfValues_ = Field<scalar>(dict.lookup("HfValues"));
    emissivityValues_ = Field<scalar>(dict.lookup("emissivityValues"));
    coordSys_ = coordinateSystem(dict, mesh_);

    Info<< "Constructed directionalSolidThermo with samples" << nl
        << "    T          : " << TValues_ << nl
        << "    rho        : " << rhoValues_ << nl
        << "    cp         : " << cpValues_ << nl
        << "    K          : " << KValues_ << nl
        << "    in coordinates system" << nl
        << "        type   : " << coordSys_.type() << nl
        << "        e3     : " << coordSys_.e3() << nl
        << "        e1     : " << coordSys_.e1() << nl
        << "    Hf         : " << HfValues_ << nl
        << "    emissivity : " << emissivityValues_ << nl
        << endl;


    if
    (
        (TValues_.size() != rhoValues_.size())
     && (TValues_.size() != cpValues_.size())
     && (TValues_.size() != rhoValues_.size())
     && (TValues_.size() != KValues_.size())
     && (TValues_.size() != HfValues_.size())
     && (TValues_.size() != emissivityValues_.size())
    )
    {
        FatalIOErrorIn("directionalSolidThermo::read()", dict)
            << "Size of property tables should be equal to size of Temperature"
            << " values " << TValues_.size()
            << exit(FatalIOError);
    }

    for (label i = 1; i < TValues_.size(); i++)
    {
        if (TValues_[i] <= TValues_[i-1])
        {
            FatalIOErrorIn("directionalSolidThermo::read()", dict)
                << "Temperature values are not in increasing order "
                << TValues_ << exit(FatalIOError);
        }
    }
    return true;
}


bool Foam::directionalSolidThermo::writeData(Ostream& os) const
{
    bool ok = basicSolidThermo::writeData(os);
    os.writeKeyword("TValues") << TValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("rhoValues") << rhoValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("cpValues") << cpValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("KValues") << KValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("HfValues") << HfValues_ << token::END_STATEMENT << nl;
    os.writeKeyword("emissivityValues") << emissivityValues_
        << token::END_STATEMENT << nl;

    return ok && os.good();
}


// * * * * * * * * * * * * * * IOStream operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const directionalSolidThermo& s)
{
    s.writeData(os);
    return os;
}


// ************************************************************************* //
