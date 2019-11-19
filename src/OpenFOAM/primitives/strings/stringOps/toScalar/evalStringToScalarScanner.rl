/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
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

Description
    Ragel lexer interface for lemon grammar of a simple string to
    scalar evaluation

\*---------------------------------------------------------------------------*/

#include "evalStringToScalarScanner.H"
#include "evalStringToScalarDriver.H"
#include "evalStringToScalarLemonParser.h"
#include "evalStringToScalarParser.H"
#include "error.H"
#include "macros.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wold-style-cast"


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel lexer with lemon parser integration

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

%%{
    machine evalScanner;
    write  data;
}%%


#undef DebugScannerInfo
#define DebugScannerInfo if (debug & 4) InfoErr


#define TOKEN_OF(T)         TOK_##T
#define EMIT_TOKEN(T)                                                         \
    driver.parsePosition() = (ts-buf);                                        \
    DebugScannerInfo << STRINGIFY(T) << ": "<< driver.parsePosition() << nl;  \
    parser_->parse(TOKEN_OF(T), 0);                                           \
    driver.parsePosition() = (p-buf);


%%{
    machine evalScanner;

    action emit_number {
        driver.parsePosition() = (ts-buf);

        DebugScannerInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver.parsePosition() << nl;

        scalar val(0);

        if (readScalar(std::string(ts, te-ts), val))
        {
            // Emit number
            parser_->parse(TOKEN_OF(NUMBER), val);
        }
        else
        {
            // Range error
            driver.reportFatal("Error parsing number");
        }

        driver.parsePosition() = (p-buf);
    }

    action emit_ident {
        driver.parsePosition() = (ts-buf);
        const word ident = word::validate(ts, te);

        driver.reportFatal("Unknown function/type: " + ident);
        driver.parsePosition() = (p-buf);
    }

    decimal = ((digit* '.' digit+) | (digit+ '.'?)) ;
    number  = (digit+ | decimal) ([Ee][\-+]? digit+)? ;
    ident   = ((alpha|'_') . ((alnum|'_')**)) ;


    ## The scanner
    main := |*
        space*;

        number => emit_number;

    ## operators
    '!'  => { EMIT_TOKEN(NOT); };
    '%'  => { EMIT_TOKEN(PERCENT); };
    '('  => { EMIT_TOKEN(LPAREN); };
    ')'  => { EMIT_TOKEN(RPAREN); };
    '*'  => { EMIT_TOKEN(TIMES); };
    '+'  => { EMIT_TOKEN(PLUS); };
    '-'  => { EMIT_TOKEN(MINUS); };
    ','  => { EMIT_TOKEN(COMMA); };
    '/'  => { EMIT_TOKEN(DIVIDE); };
    '!'  => { EMIT_TOKEN(NOT); };
    '?'  => { EMIT_TOKEN(QUESTION); };
    ':'  => { EMIT_TOKEN(COLON); };
    '<'  => { EMIT_TOKEN(LESS); };
    '<=' => { EMIT_TOKEN(LESS_EQ); };
    '>'  => { EMIT_TOKEN(GREATER); };
    '>=' => { EMIT_TOKEN(GREATER_EQ); };
    '==' => { EMIT_TOKEN(EQUAL); };
    '!=' => { EMIT_TOKEN(NOT_EQUAL); };
    '&&' => { EMIT_TOKEN(LAND); };
    '||' => { EMIT_TOKEN(LOR); };
## Not needed  '&'  => { EMIT_TOKEN(BIT_AND); };
## Not needed  '|'  => { EMIT_TOKEN(BIT_OR); };
## Not needed  '^'  => { EMIT_TOKEN(BIT_XOR); };

    ## Regular functions
    'pi'        => { EMIT_TOKEN(PI); };
    'degToRad'  => { EMIT_TOKEN(DEG_TO_RAD); };
    'radToDeg'  => { EMIT_TOKEN(RAD_TO_DEG); };
    'exp'       => { EMIT_TOKEN(EXP); };
    'log'       => { EMIT_TOKEN(LOG); };
    'log10'     => { EMIT_TOKEN(LOG10); };
    'pow'       => { EMIT_TOKEN(POW); };
    'sqr'       => { EMIT_TOKEN(SQR); };
    'sqrt'      => { EMIT_TOKEN(SQRT); };
    'cbrt'      => { EMIT_TOKEN(CBRT); };
    'sin'       => { EMIT_TOKEN(SIN); };
    'cos'       => { EMIT_TOKEN(COS); };
    'tan'       => { EMIT_TOKEN(TAN); };
    'asin'      => { EMIT_TOKEN(ASIN); };
    'acos'      => { EMIT_TOKEN(ACOS); };
    'atan'      => { EMIT_TOKEN(ATAN); };
    'atan2'     => { EMIT_TOKEN(ATAN2); };
    'hypot'     => { EMIT_TOKEN(HYPOT); };
    'sinh'      => { EMIT_TOKEN(SINH); };
    'cosh'      => { EMIT_TOKEN(COSH); };
    'tanh'      => { EMIT_TOKEN(TANH); };
    'min'       => { EMIT_TOKEN(MIN); };
    'max'       => { EMIT_TOKEN(MAX); };
    'mag'       => { EMIT_TOKEN(MAG); };
    'magSqr'    => { EMIT_TOKEN(MAGSQR); };
    'floor'     => { EMIT_TOKEN(FLOOR); };
    'ceil'      => { EMIT_TOKEN(CEIL); };
    'round'     => { EMIT_TOKEN(ROUND); };
    'rand'      => { EMIT_TOKEN(RAND); };
    'bool'      => { EMIT_TOKEN(BOOL); };

    ## Constants
    'false'     =>{ EMIT_TOKEN(BOOL_FALSE); };
    'true'      =>{ EMIT_TOKEN(BOOL_TRUE); };

    ## Catch-all for identifiers/errors
    ident       => emit_ident;

        space*;
    *|;
}%%


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parsing::evalStringToScalar::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::parsing::evalStringToScalar::scanner::process
(
    const std::string& str,
    size_t strBeg,
    size_t strLen,
    parseDriver& driver
)
{
    if (!parser_)
    {
        parser_ = new parser();
    }

    driver.content(str, strBeg, strLen);

    size_t strEnd = str.length();

    if (strBeg > str.length())
    {
        strBeg = str.length();
    }
    else if (strLen != std::string::npos)
    {
        strLen += strBeg;

        if (strLen < str.length())
        {
            strEnd = strLen;
        }
    }


    parser_->start(driver);

    // Ragel token start/end (required naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are required Ragel naming
    // - buf is our own naming

    const char* buf = &(str[strBeg]);
    const char* eof = &(str[strEnd]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    %%{write init;}%%   /* ^^^ FSM initialization here ^^^ */;

    %%{write exec;}%%  /* ^^^ FSM execution here ^^^ */;

    if (%%{write error;}%% == cs)
    {
        driver.reportFatal("Parse error while scanning", (p-buf));
    }

    if (p != eof)
    {
        driver.reportFatal("Parsing failed with remaining content", (p-buf));
    }

    // Terminate parser execution
    parser_->parse(0, 0);
    parser_->stop();

    return true;
}


// ************************************************************************* //
