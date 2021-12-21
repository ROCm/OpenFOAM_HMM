/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(Istream& is)
:
    mRows_(0),
    nCols_(0),
    v_(nullptr)
{
    this->readMatrix(is);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

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

        if (is.format() == IOstream::BINARY && is_contiguous<Type>::value)
        {
            // Binary and contiguous

            if (len)
            {
                Detail::readContiguous<Type>
                (
                    is,
                    this->data_bytes(),
                    this->size_bytes()
                );

                is.fatalCheck("readMatrix : reading the binary block");
            }

        }
        else
        {
            // Begin of contents marker
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
                            is.fatalCheck("readMatrix : reading entry");
                        }

                        is.readEndList("MatrixRow");
                    }
                }
                else  // BEGIN_BLOCK
                {
                    Type element;
                    is >> element;

                    is.fatalCheck("readMatrix : reading the single entry");

                    std::fill(begin(), end(), element);
                }
            }

            // End of contents marker
            is.readEndList("Matrix");
        }

        return len;
    }

    FatalIOErrorInFunction(is)
        << "incorrect first token, expected <int>, found "
        << firstToken.info() << nl
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
    os  << mat.nRows() << token::SPACE << mat.nCols();

    if (os.format() == IOstream::BINARY && is_contiguous<Type>::value)
    {
        // Binary and contiguous

        if (len)
        {
            // write(...) includes surrounding start/end delimiters
            os.write(mat.cdata_bytes(), mat.size_bytes());
        }
    }
    else
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
                for (label i = 0; i < mat.nRows(); ++i)
                {
                    os  << token::BEGIN_LIST;

                    // Write row
                    for (label j = 0; j < mat.nCols(); ++j)
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
                for (label i=0; i < mat.nRows(); ++i)
                {
                    os  << nl << token::BEGIN_LIST;

                    // Write row
                    for (label j = 0; j < mat.nCols(); ++j)
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
            // Empty matrix
            os  << token::BEGIN_LIST << token::END_LIST << nl;
        }
    }

    os.check(FUNCTION_NAME);
    return os;
}


// ************************************************************************* //
