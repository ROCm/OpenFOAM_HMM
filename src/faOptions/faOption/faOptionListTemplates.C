/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2020 OpenCFD Ltd.
------------------------------------------------------------------------------
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

#include "profiling.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::source
(
    GeometricField<Type, faPatchField, areaMesh>& field,
    const areaScalarField& h,
    const word& fieldName,
    const dimensionSet& ds
)
{
    checkApplied();

    tmp<faMatrix<Type>> tmtx(new faMatrix<Type>(field, ds));
    faMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(faopt, "faOption()." + source.name());

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(h, mtx, fieldi);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::operator()
(
    const areaScalarField& h,
    GeometricField<Type, faPatchField, areaMesh>& field
)
{
    return this->operator()(h, field, field.name());
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::operator()
(
    const areaScalarField& h,
    GeometricField<Type, faPatchField, areaMesh>& field,
    const word& fieldName
)
{
    return source(field, h, fieldName, field.dimensions()/dimTime*dimArea);
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::operator()
(
    const areaScalarField& h,
    const areaScalarField& rho,
    GeometricField<Type, faPatchField, areaMesh>& field
)
{
    return this->operator()(h, rho, field, field.name());
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::operator()
(
    const areaScalarField& h,
    const areaScalarField& rho,
    GeometricField<Type, faPatchField, areaMesh>& field,
    const word& fieldName
)
{
    checkApplied();

    const dimensionSet ds
    (
        rho.dimensions()*field.dimensions()/dimTime*dimArea
    );

    tmp<faMatrix<Type>> tmtx(new faMatrix<Type>(field, ds));
    faMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(faopt, "faOption()." + source.name());

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << fieldName << endl;
                }

                source.addSup(h, rho, mtx, fieldi);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::operator()
(
    const areaScalarField& rho,
    GeometricField<Type, faPatchField, areaMesh>& field,
    const dimensionSet& ds
)
{
    checkApplied();

    const dimensionSet dsMat(ds*dimArea);

    tmp<faMatrix<Type>> tmtx(new faMatrix<Type>(field, dsMat));
    faMatrix<Type>& mtx = tmtx.ref();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(field.name());

        if (fieldi != -1)
        {
            addProfiling(faopt, "faOption()." + source.name());

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying source " << source.name() << " to field "
                        << field.name() << endl;
                }

                source.addSup(rho, mtx, fieldi);
            }
        }
    }

    return tmtx;
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::d2dt2
(
    GeometricField<Type, faPatchField, areaMesh>& field
)
{
    return this->d2dt2(field, field.name());
}


template<class Type>
Foam::tmp<Foam::faMatrix<Type>> Foam::fa::optionList::d2dt2
(
    GeometricField<Type, faPatchField, areaMesh>& field,
    const word& fieldName
)
{
    return source(field, fieldName, field.dimensions()/sqr(dimTime)*dimArea);
}


template<class Type>
void Foam::fa::optionList::constrain(faMatrix<Type>& eqn)
{
    checkApplied();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(eqn.psi().name());

        if (fieldi != -1)
        {
            addProfiling(faopt, "faOption::constrain." + eqn.psi().name());

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Applying constraint " << source.name()
                        << " to field " << eqn.psi().name() << endl;
                }

                source.constrain(eqn, fieldi);
            }
        }
    }
}


template<class Type>
void Foam::fa::optionList::correct
(
    GeometricField<Type, faPatchField, areaMesh>& field
)
{
    const word& fieldName = field.name();

    forAll(*this, i)
    {
        option& source = this->operator[](i);

        label fieldi = source.applyToField(fieldName);

        if (fieldi != -1)
        {
            addProfiling(faopt, "faOption::correct." + source.name());

            source.setApplied(fieldi);

            if (source.isActive())
            {
                if (debug)
                {
                    Info<< "Correcting source " << source.name()
                        << " for field " << fieldName << endl;
                }

                source.correct(field);
            }
        }
    }
}


// ************************************************************************* //
