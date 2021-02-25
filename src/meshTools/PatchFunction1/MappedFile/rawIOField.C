/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016-2021 OpenCFD Ltd.
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

#include "rawIOField.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rawIOField<Type>::rawIOField(const IOobject& io, const bool readAverage)
:
    regIOobject(io),
    average_(Zero)
{
    // Check for MUST_READ_IF_MODIFIED
    warnNoRereading<rawIOField<Type>>();

    if
    (
        io.readOpt() == IOobject::MUST_READ
     || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
     || io.readOpt() == IOobject::READ_IF_PRESENT
    )
    {
        bool haveFile = false;
        bool headerOk = false;

        // Replacement of regIOobject::headerok() since that one complains
        // if there is no header. TBD - Move up to headerOk()/fileHandler.
        {
            const fileName fName(filePath());

            // Try to open raw first
            autoPtr<ISstream> isPtr(fileHandler().NewIFstream(fName));

            if (isPtr && isPtr->good())
            {
                haveFile = true;

                ISstream& is = isPtr();

                const token firstToken(is);

                headerOk = is.good() && firstToken.isWord("FoamFile");
            }

            isPtr.clear();

            if (debug)
            {
                Pout<< "rawIOField : object:" << io.name()
                    << " haveFile:" << haveFile
                    << " headerOk:" << headerOk << endl;
            }
        }


        if (headerOk)
        {
            // Read but don't fail upon wrong class. Could extend by providing
            // wanted typeName. Tbd.
            Istream& is = readStream(word::null);

            if (is.good())
            {
                is  >> static_cast<Field<Type>&>(*this);
                if (readAverage)
                {
                    average_ = pTraits<Type>(is);
                }
                close();
            }
        }
        else if (haveFile)
        {
            // Failed reading - fall back to IFstream
            autoPtr<ISstream> isPtr(fileHandler().NewIFstream(io.objectPath()));

            if (!isPtr || !isPtr->good())
            {
                if (io.readOpt() != IOobject::READ_IF_PRESENT)
                {
                    FatalIOErrorInFunction(*isPtr)
                        << "Trying to read raw field" << endl
                        << exit(FatalIOError);
                }
            }
            else
            {
                ISstream& is = isPtr();

                is  >> static_cast<Field<Type>&>(*this);
                if (readAverage)
                {
                    average_ = pTraits<Type>(is);
                }
            }
        }

        if (debug)
        {
            Pout<< "rawIOField : object:" << io.name()
                << " size:" << this->size() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::rawIOField<Type>::writeData(Ostream& os) const
{
    os  << static_cast<const Field<Type>&>(*this);
    if (average_ != pTraits<Type>::zero)
    {
        os << token::NL << average_;
    }
    return os.good();
}


// ************************************************************************* //
