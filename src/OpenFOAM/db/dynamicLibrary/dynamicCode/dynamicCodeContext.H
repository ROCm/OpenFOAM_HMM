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

Class
    Foam::dynamicCodeContext

Description
    Encapsulation of dynamic code dictionaries

SourceFiles
    dynamicCodeContext.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCodeContext_H
#define dynamicCodeContext_H

#include <functional>
#include "dictionary.H"
#include "SHA1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicCodeContext Declaration
\*---------------------------------------------------------------------------*/

class dynamicCodeContext
{
    // Private Data

        //- The parent dictionary context
        std::reference_wrapper<const dictionary> dict_;

        //- The SHA1 of the contents
        SHA1 sha1_;

        //- Optional "codeOptions" entry
        string options_;

        //- Optional "codeLibs" entry
        string libs_;

        //- Optional "codeInclude" entry
        string include_;

        //- Optional "code" entry
        string code_;

        //- Optional "localCode" entry
        string localCode_;


public:

    // Constructors

        //- Construct null
        dynamicCodeContext();

        //- Construct from a dictionary
        explicit dynamicCodeContext(const dictionary& dict);


    // Static Member Functions

        //- Cleanup string and expand with dictionary parameters
        static void inplaceExpand(string& code, const dictionary& dict);

        //- Prefix a \#line directive to code.
        //  The input lineNum is 0-based.
        //  Is a no-op if any of the arguments are invalid
        //  (lineNum is negative, code or file are empty)
        //
        //  \return The change in string length caused by the directive.
        //  This can potentially be used to recover the substring portions.
        static unsigned addLineDirective
        (
            string& code,
            label lineNum,
            const fileName& file
        );

        //- Prefix a \#line directive to code.
        //  The name of the dictionary is used for the 'file' name.
        static unsigned addLineDirective
        (
            string& code,
            label lineNum,
            const dictionary& dict
        );


    // Member Functions

        //- Considered valid if not using dictionary::null as the context
        bool valid() const;

        //- Set code context from a dictionary
        void setCodeContext(const dictionary& dict);

        //- Return the parent dictionary context
        const dictionary& dict() const
        {
            return dict_.get();
        }

        //- Return the code-includes
        const string& include() const
        {
            return include_;
        }

        //- Return the code-options
        const string& options() const
        {
            return options_;
        }

        //- Return the code-libs
        const string& libs() const
        {
            return libs_;
        }

        //- Return the code
        const string& code() const
        {
            return code_;
        }

        //- Return the local (file-scope) code
        const string& localCode() const
        {
            return localCode_;
        }

        //- Return SHA1 calculated from options, libs, include, code
        const SHA1& sha1() const
        {
            return sha1_;
        }

        //- Add content to SHA1 hashing
        void append(const std::string& str)
        {
            sha1_.append(str);
        }


    // Member Operators

        //- Cast to dictionary
        operator const dictionary&() const
        {
            return dict_.get();
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
