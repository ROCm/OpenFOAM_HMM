/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "Matrix.H"
#include "Istream.H"
#include "Ostream.H"
#include "token.H"
#include "contiguous.H"
#include "ListPolicy.H"
#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(Istream& is)
:
    mRows_(0),
    nCols_(0),
    v_(nullptr)
{
    operator>>(is, *this);
}


template<class Form, class Type>
bool Foam::Matrix<Form, Type>::readMatrix(Istream& is)
{
    // Anull matrix
    clear();

    is.fatalCheck(FUNCTION_NAME);

    token firstToken(is);

    is.fatalCheck("readMatrix : reading first token");

    if (firstToken.isLabel())
    {
        mRows_ = firstToken.labelToken();
        nCols_ = readLabel(is);
        doAlloc();

        // The total size
        const label len = size();

        // Read list contents depending on data format
        if (is.format() == IOstream::ASCII || !is_contiguous<Type>::value)
        {
            // Read beginning of contents
            char listDelimiter = is.readBeginList("Matrix");

            if (len)
            {
                if (listDelimiter == token::BEGIN_LIST)
                {
                    label idx = 0;

                    // Loop over rows
                    for (label i = 0; i < mRows_; ++i)
                    {
                        listDelimiter = is.readBeginList("MatrixRow");

                        for (label j = 0; j < nCols_; ++j)
                        {
                            is >> v_[idx++];

                            is.fatalCheck("readMatrix : reading reading entry");
                        }

                        is.readEndList("MatrixRow");
                    }
                }
                else
                {
                    Type element;
                    is >> element;

                    is.fatalCheck("readMatrix : reading the single entry");

                    std::fill(begin(), end(), element);
                }
            }

            // Read end of contents
            is.readEndList("Matrix");
        }
        else
        {
            if (len)
            {
                Detail::readContiguous<Type>
                (
                    is,
                    reinterpret_cast<char*>(v_),
                    len*sizeof(Type)
                );

                is.fatalCheck("readMatrix : reading the binary block");
            }
        }

        return len;
    }


    FatalIOErrorInFunction(is)
        << "incorrect first token, expected <int>, found "
        << firstToken.info()
        << exit(FatalIOError);

    return 0;
}


template<class Form, class Type>
Foam::Ostream& Foam::Matrix<Form, Type>::writeMatrix
(
    Ostream& os,
    const label shortLen
) const
{
    const Matrix<Form, Type>& mat = *this;

    // The total size
    const label len = mat.size();

    // Rows, columns size
    os  << mat.m() << token::SPACE << mat.n();

    // Write list contents depending on data format
    if (os.format() == IOstream::ASCII || !is_contiguous<Type>::value)
    {
        if (len)
        {
            const Type* v = mat.cdata();

            // Can the contents be considered 'uniform' (ie, identical)
            if (len > 1 && is_contiguous<Type>::value && mat.uniform())
            {
                // Two or more entries, and all entries have identical values.
                os  << token::BEGIN_BLOCK << v[0] << token::END_BLOCK;
            }
            else if (len < shortLen && is_contiguous<Type>::value)
            {
                // Write start contents delimiter
                os  << token::BEGIN_LIST;

                label idx = 0;

                // Loop over rows
                for (label i = 0; i < mat.m(); ++i)
                {
                    os  << token::BEGIN_LIST;

                    // Write row
                    for (label j = 0; j < mat.n(); ++j)
                    {
                        if (j) os << token::SPACE;
                        os << v[idx++];
                    }

                    os << token::END_LIST;
                }

                // Write end of contents delimiter
                os << token::END_LIST;
            }
            else
            {
                // Write start contents delimiter
                os  << nl << token::BEGIN_LIST;

                label idx = 0;

                // Loop over rows
                for (label i=0; i < mat.m(); ++i)
                {
                    os  << nl << token::BEGIN_LIST;

                    // Write row
                    for (label j=0; j < mat.n(); ++j)
                    {
                        os << nl << v[idx++];
                    }

                    os << nl << token::END_LIST;
                }

                // Write end of contents delimiter
                os << nl << token::END_LIST << nl;
            }
        }
        else
        {
            os  << token::BEGIN_LIST << token::END_LIST << nl;
        }
    }
    else
    {
        // Contents are binary and contiguous

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write
            (
                reinterpret_cast<const char*>(mat.cdata()),
                len*sizeof(Type)
            );
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}


template<class Form, class Type>
Foam::Istream& Foam::operator>>(Istream& is, Matrix<Form, Type>& mat)
{
    mat.readMatrix(is);
    return is;
}


template<class Form, class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const Matrix<Form, Type>& mat)
{
    return mat.writeMatrix(os, Detail::ListPolicy::short_length<Type>::value);
}


// ************************************************************************* //
