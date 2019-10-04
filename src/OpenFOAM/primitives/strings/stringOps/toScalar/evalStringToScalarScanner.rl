/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2019 OpenCFD Ltd.
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

Description
    Ragel lexer interface for lemon grammar of a simple string to
    scalar evaluation

\*---------------------------------------------------------------------------*/

#include "evalStringToScalarScanner.H"
#include "evalStringToScalarDriver.H"
#include "evalStringToScalarLemonParser.h"
#include "evalStringToScalarParser.H"
#include "error.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wold-style-cast"

#ifndef FULLDEBUG
#define NDEBUG
#endif

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

int Foam::parsing::evalStringToScalar::scanner::debug = 0;


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

#define TOKEN_OF(T)         TOK_##T
#define EMIT_TOKEN(T)                                                         \
    driver.parsePosition() = (ts-buf);                                        \
    DebugInfo<< "TOKEN_" #T << " at " << driver.parsePosition() << nl;        \
    parser_->parse(TOKEN_OF(T), 0);                                           \
    driver.parsePosition() = (p-buf);


%%{
    machine evalScanner;

    action emit_number {
        driver.parsePosition() = (ts-buf);

        DebugInfo
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
            // Catch range errors
            driver.reportFatal("Error parsing number");
        }

        driver.parsePosition() = (p-buf);
    }

    decimal = ((digit* '.' digit+) | (digit+ '.'?)) ;
    number  = (digit+ | decimal) ([Ee][\-+]? digit+)? ;
    lfunc   = space* '(';           # Require functions to have '('

    operators = (
        '('  @{ EMIT_TOKEN(LPAREN); }
      | ')'  @{ EMIT_TOKEN(RPAREN); }
      | '+'  @{ EMIT_TOKEN(PLUS); }
      | '-'  @{ EMIT_TOKEN(MINUS); }
      | '*'  @{ EMIT_TOKEN(TIMES); }
      | '/'  @{ EMIT_TOKEN(DIVIDE); }
      | ','  @{ EMIT_TOKEN(COMMA); }
    );

    functions = (
        'pi'         lfunc  @{ fhold; EMIT_TOKEN(PI); }
      | 'degToRad'   lfunc  @{ fhold; EMIT_TOKEN(DEG_TO_RAD); }
      | 'radToDeg'   lfunc  @{ fhold; EMIT_TOKEN(RAD_TO_DEG); }
      | 'exp'        lfunc  @{ fhold; EMIT_TOKEN(EXP); }
      | 'log'        lfunc  @{ fhold; EMIT_TOKEN(LOG); }
      | 'log10'      lfunc  @{ fhold; EMIT_TOKEN(LOG10); }
      | 'pow'        lfunc  @{ fhold; EMIT_TOKEN(POW); }
      | 'sqr'        lfunc  @{ fhold; EMIT_TOKEN(SQR); }
      | 'sqrt'       lfunc  @{ fhold; EMIT_TOKEN(SQRT); }
      | 'cbrt'       lfunc  @{ fhold; EMIT_TOKEN(CBRT); }
      | 'sin'        lfunc  @{ fhold; EMIT_TOKEN(SIN); }
      | 'cos'        lfunc  @{ fhold; EMIT_TOKEN(COS); }
      | 'tan'        lfunc  @{ fhold; EMIT_TOKEN(TAN); }
      | 'asin'       lfunc  @{ fhold; EMIT_TOKEN(ASIN); }
      | 'acos'       lfunc  @{ fhold; EMIT_TOKEN(ACOS); }
      | 'atan'       lfunc  @{ fhold; EMIT_TOKEN(ATAN); }
      | 'atan2'      lfunc  @{ fhold; EMIT_TOKEN(ATAN2); }
      | 'hypot'      lfunc  @{ fhold; EMIT_TOKEN(HYPOT); }
      | 'sinh'       lfunc  @{ fhold; EMIT_TOKEN(SINH); }
      | 'cosh'       lfunc  @{ fhold; EMIT_TOKEN(COSH); }
      | 'tanh'       lfunc  @{ fhold; EMIT_TOKEN(TANH); }
      | 'min'        lfunc  @{ fhold; EMIT_TOKEN(MIN); }
      | 'max'        lfunc  @{ fhold; EMIT_TOKEN(MAX); }
      | 'mag'        lfunc  @{ fhold; EMIT_TOKEN(MAG); }
      | 'magSqr'     lfunc  @{ fhold; EMIT_TOKEN(MAGSQR); }
      | 'floor'      lfunc  @{ fhold; EMIT_TOKEN(FLOOR); }
      | 'ceil'       lfunc  @{ fhold; EMIT_TOKEN(CEIL); }
      | 'round'      lfunc  @{ fhold; EMIT_TOKEN(ROUND); }
      | 'rand'       lfunc  @{ fhold; EMIT_TOKEN(RAND); }
    );

    main := |*
        space*;
        number => emit_number;
        operators;
        functions;
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
