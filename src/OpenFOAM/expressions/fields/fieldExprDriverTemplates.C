/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019 OpenCFD Ltd.
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
Foam::tmp<Foam::Field<Type>>
Foam::expressions::fieldExpr::parseDriver::getField
(
    const word& name
) const
{
    bool isPointVal = false;
    bool isUniformVal = false;

    tmp<Field<Type>> tfield;

    if (hasVariable(name) && variable(name).isType<Type>())
    {
        const expressions::exprResult& var = variable(name);

        isPointVal = var.isPointValue();
        isUniformVal = var.isUniform();

        tfield = var.cref<Type>().clone();
    }

    if (tfield.valid())
    {
        const label fldLen = tfield().size();
        const label len = (isPointVal ? this->pointSize() : this->size());

        if (returnReduce((fldLen == len), andOp<bool>()))
        {
            return tfield;
        }

        if (!isUniformVal)
        {
            WarningInFunction
                << "Variable " << name
                << " does not fit the size and is not a uniform value." << nl
                << "Using average value" << endl;
        }

        return tmp<Field<Type>>::New(this->size(), gAverage(tfield));
    }

    return tfield;
}


// ************************************************************************* //
