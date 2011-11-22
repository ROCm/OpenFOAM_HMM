/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "CSV.H"
#include "DynamicList.H"
#include "IFstream.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

namespace Foam
{
    // doesn't recognize specialization otherwise
    template<>
    scalar CSV<scalar>::readValue(const List<string>& splitted)
    {
        if (componentColumns_[0] >= splitted.size())
        {
            FatalErrorIn("CSV<scalar>::readValue(const List<string>&)")
                << "No column " << componentColumns_[0] << " in "
                << splitted << endl
                << exit(FatalError);
        }

        return readScalar(IStringStream(splitted[componentColumns_[0]])());
    }


    template<class Type>
    Type CSV<Type>::readValue(const List<string>& splitted)
    {
        Type result;

        for (label i = 0; i < pTraits<Type>::nComponents; i++)
        {
            if (componentColumns_[i] >= splitted.size())
            {
                FatalErrorIn("CSV<Type>::readValue(const List<string>&)")
                    << "No column " << componentColumns_[i] << " in "
                    << splitted << endl
                    << exit(FatalError);
            }

            result[i] =
                readScalar(IStringStream(splitted[componentColumns_[i]])());
        }

        return result;
    }
}


template<class Type>
void Foam::CSV<Type>::read()
{
    IFstream is(fName_.expand());

    DynamicList<Tuple2<scalar, Type> > values;

    // skip header
    if (headerLine_)
    {
        string line;
        is.getLine(line);
    }

    // read data
    while (is.good())
    {
        string line;
        is.getLine(line);

        DynamicList<string> splitted;

        std::size_t pos = 0;
        while (pos != std::string::npos)
        {
            std::size_t nPos = line.find(separator_, pos);

            if (nPos == std::string::npos)
            {
                splitted.append(line.substr(pos));
                pos = nPos;
            }
            else
            {
                splitted.append(line.substr(pos, nPos - pos));
                pos = nPos + 1;
            }
        }

        if (splitted.size() <= 1)
        {
            break;
        }

        scalar x = readScalar(IStringStream(splitted[refColumn_])());
        Type value = readValue(splitted);

        values.append(Tuple2<scalar,Type>(x, value));
    }

    table_.transfer(values);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::CSV<Type>::CSV(const word& entryName, const dictionary& dict)
:
    DataEntry<Type>(entryName),
    headerLine_(false),
    refColumn_(0),
    componentColumns_(),
    separator_(string(",")[0]),
    fName_("none"),
    table_()
{
    const dictionary coeffs(dict.subDict(type() + "Coeffs"));
    coeffs.lookup("hasHeaderLine") >> headerLine_;
    coeffs.lookup("refColumn") >> refColumn_;
    coeffs.lookup("componentColumns") >> componentColumns_;
    coeffs.readIfPresent("separator", string(separator_)[0]);
    coeffs.lookup("fileName") >> fName_;

    if (componentColumns_.size() != pTraits<Type>::nComponents)
    {
        FatalErrorIn("Foam::CSV<Type>::CSV(const word&, Istream&)")
            << componentColumns_ << " does not have the expected length "
            << pTraits<Type>::nComponents << endl
            << exit(FatalError);
    }

    read();

    if (!table_.size())
    {
        FatalErrorIn("Foam::CSV<Type>::CSV(const Istream&)")
            << "CSV for entry " << this->name_ << " is invalid (empty)"
            << nl << exit(FatalError);
    }
}


template<class Type>
Foam::CSV<Type>::CSV(const CSV<Type>& tbl)
:
    DataEntry<Type>(tbl),
    headerLine_(tbl.headerLine_),
    refColumn_(tbl.refColumn_),
    componentColumns_(tbl.componentColumns_),
    separator_(tbl.separator_),
    fName_(tbl.fName_),
    table_(tbl.table_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::CSV<Type>::~CSV()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::CSV<Type>::value(const scalar x) const
{
    // Return zero if out of bounds
    if (x < table_[0].first() || x > table_.last().first())
    {
        return pTraits<Type>::zero;
    }

    // Find i such that x(i) < x < x(i+1)
    label i = 0;
    while ((table_[i+1].first() < x) && (i+1 < table_.size()))
    {
        i++;
    }

    // Linear interpolation to find value. Note constructor needed for
    // CSV<label> to convert intermediate scalar back to label.
    return Type
    (
        (x - table_[i].first())/(table_[i+1].first() - table_[i].first())
      * (table_[i+1].second() - table_[i].second())
      + table_[i].second()
    );
}


template<class Type>
Type Foam::CSV<Type>::integrate(const scalar x1, const scalar x2) const
{
    // Initialise return value
    Type sum = pTraits<Type>::zero;

    // Return zero if out of bounds
    if ((x1 > table_.last().first()) || (x2 < table_[0].first()))
    {
        return sum;
    }

    // Find next index greater than x1
    label id1 = 0;
    while ((table_[id1].first() < x1) && (id1 < table_.size()))
    {
        id1++;
    }

    // Find next index less than x2
    label id2 = table_.size() - 1;
    while ((table_[id2].first() > x2) && (id2 >= 1))
    {
        id2--;
    }

    if ((id1 - id2) == 1)
    {
        // x1 and x2 lie within 1 interval
        sum = 0.5*(value(x1) + value(x2))*(x2 - x1);
    }
    else
    {
        // x1 and x2 cross multiple intervals

        // Integrate table body
        for (label i=id1; i<id2; i++)
        {
            sum +=
                (table_[i].second() + table_[i+1].second())
              * (table_[i+1].first() - table_[i].first());
        }
        sum *= 0.5;

        // Add table ends (partial segments)
        if (id1 > 0)
        {
            sum += 0.5
              * (value(x1) + table_[id1].second())
              * (table_[id1].first() - x1);
        }
        if (id2 < table_.size() - 1)
        {
            sum += 0.5
              * (table_[id2].second() + value(x2))
              * (x2 - table_[id2].first());
        }
    }

    return sum;
}


// * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * * * //

#include "CSVIO.C"


// ************************************************************************* //
