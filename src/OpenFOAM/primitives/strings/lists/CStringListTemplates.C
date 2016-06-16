/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenCFD Ltd.
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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class StringType>
Foam::CStringList::CStringList
(
    const UList<StringType>& input
)
:
    argc_(0),
    len_(0),
    argv_(0),
    data_(0)
{
    reset(input);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class StringType>
void Foam::CStringList::reset
(
    const UList<StringType>& input
)
{
    clear();

    argc_ = input.size();
    forAll(input, argI)
    {
        len_ += input[argI].size();
        ++len_;  // nul terminator for C-strings
    }

    argv_ = new char*[argc_+1];
    argv_[argc_] = NULL; // extra terminator

    if (argc_ > 0)
    {
        // allocation includes final nul terminator,
        // but overall count does not
        data_ = new char[len_--];

        char* ptr = data_;
        forAll(input, argI)
        {
            argv_[argI] = ptr;

            const std::string& str =
                static_cast<const std::string&>(input[argI]);

            for
            (
                std::string::const_iterator iter = str.begin();
                iter != str.end();
                ++iter
            )
            {
                *(ptr++) = *iter;
            }
            *(ptr++) = '\0';
        }
    }
}


template<class StringType>
Foam::List<StringType>
Foam::CStringList::asList(int argc, const char * const argv[])
{
    List<StringType> lst(argc);

    for (int i=0; i < argc; ++i)
    {
        lst[i] = argv[i];
    }

    return lst;
}


template<class StringType>
Foam::List<StringType>
Foam::CStringList::asList(const char * const argv[])
{
    return asList<StringType>(count(argv), argv);
}


// ************************************************************************* //
