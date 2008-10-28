/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2008 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

%{

#undef yyFlexLexer

 /* ------------------------------------------------------------------------ *\
   ------ local definitions
 \* ------------------------------------------------------------------------ */

#include "STLfileFormat.H"

using namespace Foam;

// Dummy yyFlexLexer::yylex() to keep the linker happy. It is not called
//! @cond dummy
int yyFlexLexer::yylex()
{
    FatalErrorIn("yyFlexLexer::yylex()")
        << "Should not have called this function"
        << abort(FatalError);
    return 0;
}
//! @endcond dummy

// Dummy yywrap to keep yylex happy at compile time.
// It is called by yylex but is not used as the mechanism to change file.
// See <<EOF>>
//! @cond dummy
#if YY_FLEX_SUBMINOR_VERSION < 34
extern "C" int yywrap()
#else
int yyFlexLexer::yywrap()
#endif
{
    return 1;
}
//! @endcond dummy


//- A lexer for parsing STL ASCII files.
//  Returns DynamicList(s) of points and facets (regionIds).
//  Since the facets appear within a solid/endsolid grouping,
//  they are always within a region
class STLASCIILexer
:
    public yyFlexLexer
{
    // Private data

        short groupID_;      // current region
        short maxGroupID_;   // max region
        label lineNo_;
        word  startError_;

        DynamicList<point> points_;
        DynamicList<label> facets_;
        HashTable<label> groupToPatch_;

public:

    // Constructors

        //- From input stream and the approximate number of vertices in the STL
        STLASCIILexer(istream* is, const label approxNpoints);


    // Member Functions

        //- The lexer function itself
        int lex();

    // Access

        //- A list of points corresponding to a pointField
        DynamicList<point>& points()
        {
            return points_;
        }

        //- A list of facet IDs (region IDs)
        //  corresponds to the number of triangles
        DynamicList<label>& facets()
        {
            return facets_;
        }

        //- region names
        const HashTable<label>& groupToPatch() const
        {
            return groupToPatch_;
        }

};


STLASCIILexer::STLASCIILexer(istream* is, const label approxNpoints)
:
    yyFlexLexer(is),
    groupID_(-1),
    maxGroupID_(-1),
    lineNo_(1),
    points_(approxNpoints),
    facets_(approxNpoints)
{}


 /* ------------------------------------------------------------------------ *\
   ------ cppLexer::yylex()
 \* ------------------------------------------------------------------------ */

#define YY_DECL int STLASCIILexer::lex()

%}

one_space             [ \t\f\r]
space                 {one_space}*
some_space            {one_space}+
cspace                ","{space}

alpha                 [_A-Za-z]
digit                 [0-9]
dec_digit             [0-9]
octal_digit           [0-7]
hex_digit             [0-9a-fA-F]

identifier            {alpha}({alpha}|{digit})*
integer               {dec_digit}+
label                 [1-9]{dec_digit}*
zeroLabel             {digit}*
signedInteger         [-+]?{integer}

word                  ([[:alnum:]]|[[:punct:]])*
string                {word}({some_space}{word})*

exponent_part         [eE][-+]?{digit}+
fractional_constant   [-+]?(({digit}*"."{digit}+)|({digit}+"."?))

double                (({fractional_constant}{exponent_part}?)|({digit}+{exponent_part}))
float                 {double}

x                     {float}
y                     {float}
z                     {float}

solid                 {space}("solid"|"SOLID"){space}
color                 {space}("color"|"COLOR"){some_space}{float}{some_space}{float}{some_space}{float}{space}
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

    STLpoint normal;
    STLpoint vertex;
    label cmpt = 0;   // component index used for reading vertex

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
        word group(Foam::string::validate<word>(YYText()));

        HashTable<label>::const_iterator findGroup =
            groupToPatch_.find(group);

        if (findGroup != groupToPatch_.end())
        {
            groupID_ = findGroup();
        }
        else
        {
            groupID_ = ++maxGroupID_;
            groupToPatch_.insert(group, groupID_);
        }
        BEGIN(INITIAL);
    }

<readSolidName>{space}\n {
        word group("solid");

        HashTable<label>::const_iterator findGroup =
            groupToPatch_.find(group);

        if (findGroup != groupToPatch_.end())
        {
            groupID_ = findGroup();
        }
        else
        {
            groupID_ = ++maxGroupID_;
            groupToPatch_.insert(group, groupID_);
        }

        lineNo_++;
        BEGIN(INITIAL);
    }

{color} {
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

<readVertex>{space}{float}{space} {
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
        BEGIN(INITIAL);
    }

{endsolid} {
    }


 /* ------------------ Ignore remaining space and \n s. -------------------- */

<*>{space} {}
<*>\n      { lineNo_++; }


 /* ------------------- Any other characters are errors -------------------- */

<*>. {
        startError_ = YYText();
        yy_push_state(stlError);
    }


 /* ---------------------------- Error handler ----------------------------- */

<stlError>.* {
        yy_pop_state();
        FatalErrorIn
        (
            "fileFormats::STLfileFormat::readASCII(const fileName& STLfileName)"
        )   << "while " << stateNames[YY_START] << " on line " << lineNo_ << nl
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
bool
Foam::fileFormats::STLfileFormat::readASCII
(
    IFstream& ifs,
    const off_t fileSize
)
{
    // Create the lexer with the approximate number of vertices in the STL
    // from the file size
    STLASCIILexer lexer(&ifs.stdStream(), fileSize/400);
    while (lexer.lex() != 0) {}

    DynamicList<point>& pointLst = lexer.points();
    DynamicList<label>& facetsLst = lexer.facets();

    // transfer to normal list
    points().transfer(pointLst);

    // make our triangles
    facetsLst.shrink();
    List<keyedFace> faceLst(facetsLst.size());

    face fTri(3);
    label ptI = 0;
    forAll(facetsLst, faceI)
    {
        fTri[0] = ptI++;
        fTri[1] = ptI++;
        fTri[2] = ptI++;

        faceLst[faceI] = keyedFace(fTri, facetsLst[faceI]);
    }
    facetsLst.clearStorage();

    // transfer to normal list
    faces().transfer(faceLst);
    setPatches(lexer.groupToPatch());
    stitchFaces(SMALL);

    return true;
}


 /* ------------------------------------------------------------------------ *\
    ------ End of STLfileFormatASCII.L
 \* ------------------------------------------------------------------------ */
