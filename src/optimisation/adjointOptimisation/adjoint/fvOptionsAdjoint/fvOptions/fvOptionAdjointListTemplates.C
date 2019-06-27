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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fv::optionAdjointList::correct
(
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    const word& fieldName = fld.name();

    forAll(*this, i)
    {
        optionAdjoint& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Correcting source " << source.name()
                        << " for field " << fieldName << endl;
                }

                source.correct(fld);
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionAdjointList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    return this->operator()(fld, fld.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionAdjointList::operator()
(
    GeometricField<Type, fvPatchField, volMesh>& fld,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds = fld.dimensions()/dimTime*dimVolume;

    tmp<fvMatrix<Type>> tmtx(new fvMatrix<Type>(fld, ds));
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        optionAdjoint& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(mtx, fieldI);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionAdjointList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    return this->operator()(rho, fld, fld.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionAdjointList::operator()
(
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& fld,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds = rho.dimensions()*fld.dimensions()/dimTime*dimVolume;

    tmp<fvMatrix<Type>> tmtx(new fvMatrix<Type>(fld, ds));
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        optionAdjoint& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(rho, mtx, fieldI);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionAdjointList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& fld
)
{
    return this->operator()(alpha, rho, fld, fld.name());
}


template<class Type>
Foam::tmp<Foam::fvMatrix<Type>> Foam::fv::optionAdjointList::operator()
(
    const volScalarField& alpha,
    const volScalarField& rho,
    GeometricField<Type, fvPatchField, volMesh>& fld,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds =
        alpha.dimensions()*rho.dimensions()*fld.dimensions()/dimTime*dimVolume;

    tmp<fvMatrix<Type>> tmtx(new fvMatrix<Type>(fld, ds));
    fvMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        optionAdjoint& source = this->operator[](i);

        label fieldI = source.applyToField(fieldName);

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(alpha, rho, mtx, fieldI);
            }
        }
    }

    return tmtx;
}


template<class Type>
void Foam::fv::optionAdjointList::constrain(fvMatrix<Type>& eqn)
{
    checkApplied();

    forAll(*this, i)
    {
        optionAdjoint& source = this->operator[](i);

        label fieldI = source.applyToField(eqn.psi().name());

        if (fieldI != -1)
        {
            source.setApplied(fieldI);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying constraint " << source.name()
                        << " to field " << eqn.psi().name() << endl;
                }

                source.constrain(eqn, fieldI);
            }
        }
    }
}


// ************************************************************************* //
