/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016-2017 OpenCFD Ltd.
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

#include "CSV.H"
#include "DynamicList.H"
//#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<>
Foam::label Foam::Function1Types::CSV<Foam::label>::readValue
(
    const List<string>& splitted
)
{
    if (componentColumns_[0] >= splitted.size())
    {
        FatalErrorInFunction
            << "No column " << componentColumns_[0] << " in "
            << splitted << endl
            << exit(FatalError);
    }

    return readLabel(splitted[componentColumns_[0]]);
}


template<>
Foam::scalar Foam::Function1Types::CSV<Foam::scalar>::readValue
(
    const List<string>& splitted
)
{
    if (componentColumns_[0] >= splitted.size())
    {
        FatalErrorInFunction
            << "No column " << componentColumns_[0] << " in "
            << splitted << endl
            << exit(FatalError);
    }

    return readScalar(splitted[componentColumns_[0]]);
}


template<class Type>
Type Foam::Function1Types::CSV<Type>::readValue(const List<string>& splitted)
{
    Type result;

    for (label i = 0; i < pTraits<Type>::nComponents; ++i)
    {
        if (componentColumns_[i] >= splitted.size())
        {
            FatalErrorInFunction
                << "No column " << componentColumns_[i] << " in "
                << splitted << endl
                << exit(FatalError);
        }

        result[i] = readScalar(splitted[componentColumns_[i]]);
    }

    return result;
}


template<class Type>
void Foam::Function1Types::CSV<Type>::read()
{
    fileName expandedFile(fName_);
    //IFstream is(expandedFile.expand());
    autoPtr<ISstream> isPtr(fileHandler().NewIFstream(expandedFile.expand()));
    ISstream& is = isPtr();

    if (!is.good())
    {
        FatalIOErrorInFunction(is)
            << "Cannot open CSV file for reading."
            << exit(FatalIOError);
    }

    DynamicList<Tuple2<scalar, Type>> values;

    // skip header
    for (label i = 0; i < nHeaderLine_; i++)
    {
        string line;
        is.getLine(line);
    }

    label nEntries = max(componentColumns_);

    // read data
    while (is.good())
    {
        string line;
        is.getLine(line);


        label n = 0;
        std::size_t pos = 0;
        DynamicList<string> splitted;

        if (mergeSeparators_)
        {
            std::size_t nPos = 0;

            while ((pos != std::string::npos) && (n <= nEntries))
            {
                bool found = false;
                while (!found)
                {
                    nPos = line.find(separator_, pos);

                    if ((nPos != std::string::npos) && (nPos - pos == 0))
                    {
                        pos = nPos + 1;
                    }
                    else
                    {
                        found = true;
                    }
                }

                nPos = line.find(separator_, pos);

                if (nPos == std::string::npos)
                {
                    splitted.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    splitted.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }
        else
        {
            while ((pos != std::string::npos) && (n <= nEntries))
            {
                std::size_t nPos = line.find(separator_, pos);

                if (nPos == std::string::npos)
                {
                    splitted.append(line.substr(pos));
                    pos = nPos;
                    n++;
                }
                else
                {
                    splitted.append(line.substr(pos, nPos - pos));
                    pos = nPos + 1;
                    n++;
                }
            }
        }


        if (splitted.size() <= 1)
        {
            break;
        }

        scalar x = readScalar(splitted[refColumn_]);
        Type value = readValue(splitted);

        values.append(Tuple2<scalar,Type>(x, value));
    }

    this->table_.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::CSV<Type>::CSV
(
    const word& entryName,
    const dictionary& dict,
    const fileName& fName
)
:
    TableBase<Type>(entryName, dict),
    nHeaderLine_(readLabel(dict.lookup("nHeaderLine"))),
    refColumn_(readLabel(dict.lookup("refColumn"))),
    componentColumns_(dict.lookup("componentColumns")),
    separator_(dict.lookupOrDefault<string>("separator", string(","))[0]),
    mergeSeparators_(readBool(dict.lookup("mergeSeparators"))),
    fName_(fName != fileName::null ? fName : dict.lookup("file"))
{
    if (componentColumns_.size() != pTraits<Type>::nComponents)
    {
        FatalErrorInFunction
            << componentColumns_ << " does not have the expected length of "
            << pTraits<Type>::nComponents << endl
            << exit(FatalError);
    }

    read();

    TableBase<Type>::check();
}


template<class Type>
Foam::Function1Types::CSV<Type>::CSV(const CSV<Type>& csv)
:
    TableBase<Type>(csv),
    nHeaderLine_(csv.nHeaderLine_),
    refColumn_(csv.refColumn_),
    componentColumns_(csv.componentColumns_),
    separator_(csv.separator_),
    mergeSeparators_(csv.mergeSeparators_),
    fName_(csv.fName_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::CSV<Type>::~CSV()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
const Foam::fileName& Foam::Function1Types::CSV<Type>::fName() const
{
    return fName_;
}


template<class Type>
void Foam::Function1Types::CSV<Type>::writeData(Ostream& os) const
{
    Function1<Type>::writeData(os);
    os.endEntry();

    os.beginBlock(word(this->name() + "Coeffs"));

    // Note: for TableBase write the dictionary entries it needs but not
    // the values themselves
    TableBase<Type>::writeEntries(os);

    os.writeEntry("nHeaderLine", nHeaderLine_);
    os.writeEntry("refColumn", refColumn_);

    // Force writing labelList in ascii
    const enum IOstream::streamFormat fmt = os.format();
    os.format(IOstream::ASCII);
    os.writeEntry("componentColumns", componentColumns_);
    os.format(fmt);

    os.writeEntry("separator", string(separator_));
    os.writeEntry("mergeSeparators", mergeSeparators_);
    os.writeEntry("file", fName_);

    os.endBlock() << flush;
}


// ************************************************************************* //
