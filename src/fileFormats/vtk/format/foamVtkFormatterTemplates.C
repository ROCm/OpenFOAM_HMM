/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "foamVtkPTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::vtk::formatter&
Foam::vtk::formatter::writeAttribute
(
    const word& k,
    const Type& v,
    const char quote
)
{
    if (!inTag_)
    {
        WarningInFunction
            << "xml attribute '" << k << "' but not within a tag!"
            << endl;
    }

    os_ << ' ' << k << '=' << quote << v << quote;

    return *this;
}


template<class Type, int nComp>
Foam::vtk::formatter&
Foam::vtk::formatter::openDataArray
(
    const word& dataName
)
{
    openTag("DataArray");
    xmlAttr("type", vtkPTraits<Type>::typeName);
    xmlAttr("Name", dataName);
    if (nComp > 1)
    {
        xmlAttr(fileAttr::NUMBER_OF_COMPONENTS, nComp);
    }
    xmlAttr("format", name());

    return *this;
}


template<class Type, int nComp>
Foam::vtk::formatter&
Foam::vtk::formatter::openDataArray
(
    const vtk::dataArrayAttr& attrEnum
)
{
    return openDataArray<Type, nComp>(vtk::dataArrayAttrNames[attrEnum]);
}


template<class Type, int nComp>
Foam::vtk::formatter&
Foam::vtk::formatter::PDataArray
(
    const word& dataName
)
{
    openTag("PDataArray");
    xmlAttr("type", vtkPTraits<Type>::typeName);
    if (dataName.size())
    {
        xmlAttr("Name", dataName);
    }
    if (nComp > 1)
    {
        xmlAttr(fileAttr::NUMBER_OF_COMPONENTS, nComp);
    }

    closeTag(true);

    return *this;
}


// ************************************************************************* //
