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

Class
    Foam::CSV

Description
    Templated CSV container data entry.  Reference column is always a scalar,
    e.g. time

    \verbatim
        <entryName>   CSV
        {
            hasHeaderLine   true;
            refColumn       0;              // reference column index
            componentColumns (0 1 2);       // component column indices
            separator       ",";            // optional (defaults to ",")
            fileName        fileXYZ;        // name of csv data file
        }
    \endverbatim

SourceFiles
    CSV.C

\*---------------------------------------------------------------------------*/

#ifndef CSV_H
#define CSV_H

#include "DataEntry.H"
#include "Tuple2.H"
#include "labelList.H"
#include "ISstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class CSV;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const CSV<Type>&
);

/*---------------------------------------------------------------------------*\
                           Class CSV Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class CSV
:
    public DataEntry<Type>
{
    // Private data

        //- Does the file have a header line?
        bool headerLine_;

        //- Column of the time
        label refColumn_;

        //- Labels of the components
        labelList componentColumns_;

        //- Separator character
        char separator_;

        //- File name for csv table (optional)
        fileName fName_;

        //- CSV data
        List<Tuple2<scalar, Type> > table_;


    // Private Member Functions

        //- Read csv data table
        void read();

        //- Read the next value from the splitted string
        Type readValue(const List<string>&);

        //- Disallow default bitwise assignment
        void operator=(const CSV<Type>&);


public:

    //- Runtime type information
    TypeName("csv");


    // Constructors

        //- Construct from entry name and Istream
        CSV(const word& entryName, const dictionary& dict);

        //- Copy constructor
        CSV(const CSV<Type>& tbl);

        //- Construct and return a clone
        virtual tmp<DataEntry<Type> > clone() const
        {
            return tmp<DataEntry<Type> >(new CSV<Type>(*this));
        }


    //- Destructor
    virtual ~CSV();


    // Member Functions

        //- Return CSV value
        Type value(const scalar x) const;

        //- Integrate between two (scalar) values
        Type integrate(const scalar x1, const scalar x2) const;


    // I/O

        //- Ostream Operator
        friend Ostream& operator<< <Type>
        (
            Ostream& os,
            const CSV<Type>& cnst
        );

        //- Write in dictionary format
        virtual void writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "CSV.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
