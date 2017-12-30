/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

%option prefix="yySTL"
%option yyclass="yySTLFlexLexer"

%{

 /* ------------------------------------------------------------------------ *\
   ------ local definitions
 \* ------------------------------------------------------------------------ */

#include "STLReader.H"
#include "OSspecific.H"

using namespace Foam;
// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! \cond dummy
#if YY_FLEX_MAJOR_VERSION <= 2 && YY_FLEX_MINOR_VERSION <= 5 && YY_FLEX_SUBMINOR_VERSION < 34
int yyFlexLexer::yylex()
{
    FatalErrorInFunction
        << "Should not have called this function"
        << abort(FatalError);
    return 0;
}
#endif
//! \endcond

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! \cond dummy
#if YY_FLEX_MAJOR_VERSION <= 2 && YY_FLEX_MINOR_VERSION <= 5 && YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yySTLFlexLexer::yywrap()
#endif
{
    return 1;
}
//! \endcond

//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (zoneIds).
//  The facets are within a solid/endsolid grouping
class STLASCIILexer
:
    public yySTLFlexLexer
{
    // Private data

        bool  sorted_;
        label groupID_;      // current solid group
        label lineNo_;
        word  startError_;

        DynamicList<STLpoint> points_;
        DynamicList<label> facets_;
        DynamicList<word>  names_;
        DynamicList<label> sizes_;
        HashTable<label>   lookup_;

public:

    // Constructors

        //- From input stream and the approximate number of vertices in the STL
        STLASCIILexer(istream* is, const label approxNpoints);


    // Member Functions

        //- The lexer function itself
        int lex();

      // Access

        //- Do all the solid groups appear in order?
        inline bool sorted() const
        {
            return sorted_;
        }

        //- A list of unstitched triangle points
        inline DynamicList<STLpoint>& points()
        {
            return points_;
        }

        //- A list of facet IDs (group IDs)
        //  corresponds to the number of triangles
        inline DynamicList<label>& facets()
        {
            return facets_;
        }

        //- Solid names in the order of their appearance.
        inline DynamicList<word>& names()
        {
            return names_;
        }

        //- Solid sizes in the order of their appearance.
        inline DynamicList<label>& sizes()
        {
            return sizes_;
        }
};


STLASCIILexer::STLASCIILexer(istream* is, const label approxNpoints)
:
    yySTLFlexLexer(is),
    sorted_(true),
    groupID_(-1),
    lineNo_(1),
    points_(approxNpoints),
    facets_(approxNpoints)
{}


 /* ------------------------------------------------------------------------ *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------ */

#undef YY_DECL
#define YY_DECL int STLASCIILexer::lex()

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+

alpha                 [_A-Za-z]
digit                 [0-9]

integer               {digit}+
signedInteger         [-+]?{integer}

word                  ([[:alnum:]]|[[:punct:]])*
string                {word}({some_space}{word})*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

floatNum              (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))

x                     {floatNum}
y                     {floatNum}
z                     {floatNum}

solid                 {space}("solid"|"SOLID"){space}
color                 {space}("color"|"COLOR"){some_space}{floatNum}{some_space}{floatNum}{some_space}{floatNum}{space}
facet                 {space}("facet"|"FACET"){space}
normal                {space}("normal"|"NORMAL"){space}
point                 {space}{x}{some_space}{y}{some_space}{z}{space}
outerloop             {space}("outer"{some_space}"loop")|("OUTER"{some_space}"LOOP"){space}
vertex                {space}("vertex"|"VERTEX"){space}
endloop               {space}("endloop"|"ENDLOOP"){space}
endfacet              {space}("endfacet"|"ENDFACET"){space}
endsolid              {space}("endsolid"|"ENDSOLID")({some_space}{word})*


 /* ------------------------------------------------------------------------ *\
                      -----  Exclusive start states -----
 \* ------------------------------------------------------------------------ */

%option stack

%x readSolidName
%x readFacet
%x readNormal
%x readVertices
%x readVertex
%x stlError

%%

%{
    // End of read character pointer returned by strtof
    // char* endPtr;

    label cmpt = 0;   // Component index when reading vertex
    STLpoint vertex;
    // STLpoint normal;

    static const char* stateNames[7] =
    {
        "reading solid",
        "reading solid name",
        "reading facet",
        "reading normal",
        "reading vertices",
        "reading vertex",
        "error"
    };

    static const char* stateExpects[7] =
    {
        "'solid', 'color', 'facet' or 'endsolid'",
        "<string>",
        "'normal', 'outer loop' or 'endfacet'",
        "<float> <float> <float>",
        "'vertex' or 'endloop'",
        "<float> <float> <float>",
        ""
    };
%}


 /* ------------------------------------------------------------------------ *\
                            ------ Start Lexing ------
 \* ------------------------------------------------------------------------ */

 /*                      ------ Reading control header ------                */

{solid} {
    BEGIN(readSolidName);
}

<readSolidName>{string} {
    const word solidName(word::validate(YYText()));

    auto iter = lookup_.cfind(solidName);
    if (iter.found())
    {
        if (groupID_ != iter.object())
        {
            sorted_ = false;  // Group appeared out of order
            groupID_ = iter.object();
        }
    }
    else
    {
        groupID_ = sizes_.size();
        if (lookup_.insert(solidName, groupID_))
        {
            names_.append(solidName);
            sizes_.append(0);
        }
        else
        {
            FatalErrorInFunction<< "Duplicate solid-name: " << solidName
                << exit(FatalError);
        }
    }
    BEGIN(INITIAL);
}

<readSolidName>{space}\n {
    const word solidName("solid");  // Could also use solid0, solid1, ...

    auto iter = lookup_.cfind(solidName);
    if (iter.found())
    {
        if (groupID_ != iter.object())
        {
            sorted_ = false;  // Group appeared out of order
            groupID_ = iter.object();
        }
    }
    else
    {
        groupID_ = sizes_.size();
        if (lookup_.insert(solidName, groupID_))
        {
            names_.append(solidName);
            sizes_.append(0);
        }
        else
        {
            FatalErrorInFunction<< "Duplicate solid-name: " << solidName
                << exit(FatalError);
        }
    }

    ++lineNo_;
    BEGIN(INITIAL);
}

{color} {
    /* ignore 'color' */
}

{facet} {
    BEGIN(readFacet);
}

<readFacet>{normal} {
    BEGIN(readNormal);
}

<readNormal>{point} {
    /*
     skip reading normals:
     normal.x() = strtof(YYText(), &endPtr);
     normal.y() = strtof(endPtr, &endPtr);
     normal.z() = strtof(endPtr, &endPtr);
     normals_.append(normal);
     */
    BEGIN(readFacet);
}

<readFacet>{outerloop} {
    BEGIN(readVertices);
}

<readVertices>{vertex} {
    BEGIN(readVertex);
}

<readVertex>{space}{signedInteger}{space} {
    vertex[cmpt++] = atol(YYText());

    if (cmpt == 3)
    {
        cmpt = 0;
        points_.append(vertex);
        BEGIN(readVertices);
    }
}

<readVertex>{space}{floatNum}{space} {
    vertex[cmpt++] = atof(YYText());

    if (cmpt == 3)
    {
        cmpt = 0;
        points_.append(vertex);
        BEGIN(readVertices);
    }
}

<readVertices>{endloop} {
    BEGIN(readFacet);
}

<readFacet>{endfacet} {
    facets_.append(groupID_);
    sizes_[groupID_]++;
    BEGIN(INITIAL);
}

{endsolid} {
}


 /* ---------------- Ignore remaining spaces and newlines ------------------ */

<*>{space} {}
<*>\n      { ++lineNo_; }


 /* ------------------- Any other characters are errors -------------------- */

<*>. {
    startError_ = YYText();
    yy_push_state(stlError);
}


 /* ---------------------------- Error handler ----------------------------- */

<stlError>.* {
    yy_pop_state();
    FatalErrorInFunction
        << "while " << stateNames[YY_START] << " on line " << lineNo_ << nl
        << "    expected " << stateExpects[YY_START]
        << " but found '" << startError_.c_str() << YYText() << "'"
        << exit(FatalError);
}


 /*  ------------------------ On EOF terminate ----------------------------  */

<<EOF>> {
    yyterminate();
}
%%


//
// member function
//
bool Foam::fileFormats::STLReader::readASCII
(
    const fileName& filename
)
{
    format_ = STLFormat::UNKNOWN;

    IFstream is(filename);
    if (!is)
    {
        FatalErrorInFunction
            << "file " << filename << " not found"
            << exit(FatalError);
    }

    // Create the lexer with the approximate number of vertices in the STL
    // from the file size
    STLASCIILexer lexer(&(is.stdStream()), Foam::fileSize(filename)/400);
    while (lexer.lex() != 0) {}

    sorted_ = lexer.sorted();

    // Transfer to normal lists
    points_.transfer(lexer.points());
    zoneIds_.transfer(lexer.facets());
    names_.transfer(lexer.names());
    sizes_.transfer(lexer.sizes());

    format_ = STLFormat::ASCII;
    return true;
}

 /* ------------------------------------------------------------------------ *\
    ------ End of STLReaderASCII.L
 \* ------------------------------------------------------------------------ */
