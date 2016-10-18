/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "vtkTools.H"
#include "endian.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::vtkTools::swapWord(label& word32)
{
    char* mem = reinterpret_cast<char*>(&word32);

    char a = mem[0];
    mem[0] = mem[3];
    mem[3] = a;

    a = mem[1];
    mem[1] = mem[2];
    mem[2] = a;
}


void Foam::vtkTools::swapWords(const label nWords, label* words32)
{
    for (label i = 0; i < nWords; i++)
    {
        swapWord(words32[i]);
    }
}


void Foam::vtkTools::write
(
    std::ostream& os,
    const bool binary,
    List<floatScalar>& fField
)
{
    if (binary)
    {
        #ifdef WM_LITTLE_ENDIAN
        swapWords(fField.size(), reinterpret_cast<label*>(fField.begin()));
        #endif
        os.write
        (
            reinterpret_cast<char*>(fField.begin()),
            fField.size()*sizeof(float)
        );

        os  << std::endl;
    }
    else
    {
        forAll(fField, i)
        {
            os  << fField[i];

            if (i > 0 && (i % 10) == 0)
            {
                os  << std::endl;
            }
            else
            {
                os  << ' ';
            }
        }
        os  << std::endl;
    }
}


void Foam::vtkTools::write
(
    std::ostream& os,
    const bool binary,
    DynamicList<floatScalar>& fField
)
{
    List<floatScalar>& fld = fField.shrink();

    write(os, binary, fld);
}


void Foam::vtkTools::write
(
    std::ostream& os,
    const bool binary,
    labelList& elems
)
{
    if (binary)
    {
        #ifdef WM_LITTLE_ENDIAN
        swapWords(elems.size(), reinterpret_cast<label*>(elems.begin()));
        #endif
        os.write
        (
            reinterpret_cast<char*>(elems.begin()),
            elems.size()*sizeof(label)
        );

        os  << std::endl;
    }
    else
    {
        forAll(elems, i)
        {
            os  << elems[i];

            if (i > 0 && (i % 10) == 0)
            {
                os  << std::endl;
            }
            else
            {
                os  << ' ';
            }
        }
        os  << std::endl;
    }
}


void Foam::vtkTools::write
(
    std::ostream& os,
    const bool binary,
    DynamicList<label>& elems
)
{
    labelList& fld = elems.shrink();

    write(os, binary, fld);
}


void Foam::vtkTools::writeHeader
(
    std::ostream& os,
    const bool binary,
    const std::string& title
)
{
    os  << "# vtk DataFile Version 2.0" << std::endl
        << title << std::endl;

    if (binary)
    {
        os  << "BINARY" << std::endl;
    }
    else
    {
        os  << "ASCII" << std::endl;
    }
}


void Foam::vtkTools::writeCellDataHeader
(
    std::ostream& os,
    const label nCells,
    const label nFields
)
{
    os  << "CELL_DATA " << nCells << std::endl
        << "FIELD attributes " << nFields << std::endl;
}


void Foam::vtkTools::writePointDataHeader
(
    std::ostream& os,
    const label nPoints,
    const label nFields
)
{
    os  << "POINT_DATA  " << nPoints << std::endl
        << "FIELD attributes " << nFields << std::endl;
}


void Foam::vtkTools::insert(const scalar src, DynamicList<floatScalar>& dest)
{
    dest.append(float(src));
}


void Foam::vtkTools::insert(const vector& src, DynamicList<floatScalar>& dest)
{
    for (direction cmpt = 0; cmpt < vector::nComponents; ++cmpt)
    {
        dest.append(float(src[cmpt]));
    }
}


void Foam::vtkTools::insert
(
    const sphericalTensor& src,
    DynamicList<floatScalar>& dest
)
{
    for (direction cmpt = 0; cmpt < sphericalTensor::nComponents; ++cmpt)
    {
        dest.append(float(src[cmpt]));
    }
}


void Foam::vtkTools::insert
(
    const symmTensor& src,
    DynamicList<floatScalar>& dest
)
{
    dest.append(float(src.xx()));
    dest.append(float(src.yy()));
    dest.append(float(src.zz()));
    dest.append(float(src.xy()));
    dest.append(float(src.yz()));
    dest.append(float(src.xz()));
}


void Foam::vtkTools::insert(const tensor& src, DynamicList<floatScalar>& dest)
{
    for (direction cmpt = 0; cmpt < tensor::nComponents; ++cmpt)
    {
        dest.append(float(src[cmpt]));
    }
}


void Foam::vtkTools::insert(const labelList& src, DynamicList<label>& dest)
{
    dest.append(src);
}


// ************************************************************************* //
