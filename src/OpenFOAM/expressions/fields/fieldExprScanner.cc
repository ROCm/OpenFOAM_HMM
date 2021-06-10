
#line 1 "fieldExprScanner.rl"
/*--------------------------------*- C++ -*----------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
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
    Ragel lexer interface for lemon grammar for patch expressions

\*---------------------------------------------------------------------------*/

#include "fieldExprScanner.H"
#include "fieldExprDriver.H"
#include "fieldExprLemonParser.h"
#include "fieldExprParser.H"
#include "Enum.H"
#include "macros.H"

#pragma GCC diagnostic ignored "-Wimplicit-fallthrough"
#pragma GCC diagnostic ignored "-Wunused-const-variable"

// Debugging to stderr
#undef  DebugInfo
#define DebugInfo if (debug & 0x2) InfoErr


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

//- Paste token prefix
#define TOKEN_OF(T)         TOK_##T

//- An {int, c_str} enum pairing
#define TOKEN_PAIR(Name,T)  { TOKEN_OF(T), Name }

//- An {int, c_str} enum pairing for field types
#define FIELD_PAIR(Fld,T)  { TOKEN_OF(T), Fld::typeName.c_str() }

#undef HAS_LOOKBEHIND_TOKENS

// Special handling of predefined method types. Eg, .x(), .y(), ...
static const Enum<int> fieldMethodEnums
({
    TOKEN_PAIR("x", CMPT_X),
    TOKEN_PAIR("y", CMPT_Y),
    TOKEN_PAIR("z", CMPT_Z),
    TOKEN_PAIR("xx", CMPT_XX),
    TOKEN_PAIR("xy", CMPT_XY),
    TOKEN_PAIR("xz", CMPT_XZ),
    TOKEN_PAIR("yx", CMPT_YX),
    TOKEN_PAIR("yy", CMPT_YY),
    TOKEN_PAIR("yz", CMPT_YZ),
    TOKEN_PAIR("zx", CMPT_ZX),
    TOKEN_PAIR("zy", CMPT_ZY),
    TOKEN_PAIR("zz", CMPT_ZZ),
    TOKEN_PAIR("ii", CMPT_II),
    TOKEN_PAIR("diag", DIAG),   /* tensors only */
    TOKEN_PAIR("T", TRANSPOSE), /* tensors only */
});


// Simple compile-time function name declarations.
// Useful for handling driver-specific dispatching, or functions that
// are not universally available.
static const Enum<int> funcTokenEnums
({
#ifdef TOK_FLOOR
    TOKEN_PAIR("floor", FLOOR),
    TOKEN_PAIR("ceil", CEIL),
    TOKEN_PAIR("round", ROUND),
#endif
#ifdef TOK_HYPOT
    TOKEN_PAIR("hypot", HYPOT),
#endif
});

} // End namespace Foam


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace Foam
{

// Classifying token type based on an identifier name is indeed ugly.
//
// 1)
//   Handle special cases (eg, cellSet,...) first that have been tagged
//   as expected content with the stashed "look-behind" token.
//   Handle not-found errors here directly.
//
// 2)
//   Fallback to determining which field-type (volScalarField etc) the name
//   corresponds to.
//   Handle not-found errors by return -1.
//
static int driverTokenType
(
    const expressions::fieldExpr::parseDriver& driver_,
    const word& ident
)
{
    // Field variables
    #ifdef TOK_SCALAR_ID
    {
        #undef checkFieldToken
        #define checkFieldToken(TokType, Type)                                \
        if (driver_.isLocalVariable<Type>(ident, false))                      \
        {                                                                     \
            return TokType;                                                   \
        }

        checkFieldToken(TOK_SCALAR_ID, scalar);
        checkFieldToken(TOK_VECTOR_ID, vector);
        checkFieldToken(TOK_SYM_TENSOR_ID, symmTensor);
        checkFieldToken(TOK_SPH_TENSOR_ID, sphericalTensor);
        checkFieldToken(TOK_TENSOR_ID, tensor);
        // Not tested: checkFieldToken(TOK_BOOL_ID, bool);
    }
    #endif

    return -1;
}

} // End namespace Foam


/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

// Ragel machine definition
// Ragel variables (p, pe, eof, cs, top, stack, ts, te, act) defined later...
//
// Can use 'variable p xxx;' etc to change these names

#define EMIT_TOKEN(T)                                                         \
    driver_.parsePosition() = (ts-buf);                                       \
    DebugInfo<< STRINGIFY(T) << " at " << driver_.parsePosition() << nl;      \
    parser_->parse(TOKEN_OF(T), nullptr);                                     \
    driver_.parsePosition() = (p-buf);



#line 167 "fieldExprScanner.cc"
static const int fieldExpr_start = 11;
static const int fieldExpr_first_final = 11;
static const int fieldExpr_error = 0;

static const int fieldExpr_en_main = 11;


#line 305 "fieldExprScanner.rl"



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::expressions::fieldExpr::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::expressions::fieldExpr::scanner::dispatch_method
(
    const parseDriver& driver_,
    scanToken& scanTok,
    word&& ident
) const
{
    if (ident[0] == '.')
    {
        ident.erase(0, 1);
    }

    DebugInfo
        << "Method:" << ident
        << " at " << driver_.parsePosition() << nl;

    const int methType = fieldMethodEnums.lookup(ident, -1);

    if (methType > 0)
    {
        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT, nullptr);
        parser_->parse(methType, nullptr);

        return true;
    }

    driver_.reportFatal("Unknown method: " + ident);
    return false;
}


bool Foam::expressions::fieldExpr::scanner::dispatch_ident
(
    const parseDriver& driver_,
    scanToken& scanTok,
    word&& ident
) const
{
    int tokType = -1;

    const bool quoted =
    (
        (ident.front() == '"' || ident.front() == '\'')
     && (ident.front() == ident.back())
    );

    if (quoted)
    {
        ident.erase(ident.size()-1);
        ident.erase(0, 1);
    }
    else
    {
        // Check for function name
        tokType = funcTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " function:"
                << parser_->tokenName(tokType) << nl;

            parser_->parse(tokType, nullptr);
            return true;
        }

        #ifdef HAS_LOOKBEHIND_TOKENS
        // Specials such "cset" also reset the look-behind
        tokType = lookBehindTokenEnums.lookup(ident, -1);

        if (tokType > 0)
        {
            DebugInfo
                << "Emit:" << ident << " as look-behind:"
                << parser_->tokenName(tokType) << nl;

            driver_.resetStashedTokenId(tokType);
            parser_->parse(tokType, nullptr);
            return true;
        }
        #endif
    }


    // Can also peek at stashed "look-behind"
    // const int lookBehind = driver_.stashedTokenId();

    tokType = driverTokenType(driver_, ident);

    if (tokType > 0)
    {
        DebugInfo
            << "Emit:" << ident << " token:"
            << parser_->tokenName(tokType) << nl;

        scanTok.name = new Foam::word(std::move(ident));
        parser_->parse(tokType, &scanTok);

        return true;
    }


    // Not found? Attempt to strip off '.x' endings etc,
    // but not when quoted

    const auto dot = ident.rfind('.');
    const int methType =
    (
        quoted || dot == std::string::npos
      ? -1
      : fieldMethodEnums.lookup(ident.substr(dot+1), -1)
    );

    if
    (
        methType > 0
     && (tokType = driverTokenType(driver_, ident.substr(0, dot))) > 0
    )
    {
        DebugInfo
            << "Emit:" << ident.substr(0, dot).c_str() << " token:"
            << parser_->tokenName(tokType) << " with "
            << ident.substr(dot).c_str() << " token:"
            << parser_->tokenName(methType) << nl;

        // The field (before the ".")
        ident.erase(dot);

        scanTok.name = new Foam::word(std::move(ident));
        parser_->parse(tokType, &scanTok);

        // Dispatch '.' and "method" separately
        parser_->parse(TOK_DOT, nullptr);
        parser_->parse(methType, nullptr);

        return true;
    }

    driver_.reportFatal
    (
        "Object " + ident + " does not exist or wrong type"
    );

    return false;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::expressions::fieldExpr::scanner::process
(
    const std::string& str,
    size_t strBeg,
    size_t strLen,
    parseDriver& driver_
)
{
    // Save debug value
    const int oldDebug = debug;

    if (driver_.debugScanner()) { debug |= 0x2; }
    if (driver_.debugParser())  { debug |= 0x4; }

    if (debug & 0x6)
    {
        InfoErr
            << "Begin parse {"
            << str.substr(strBeg, strLen).c_str() << '}' << nl;
    }

    if (!parser_)
    {
        parser_ = new parser();
    }

    driver_.content(str, strBeg, strLen);

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


    parser_->start(driver_);

    // Scan token type
    scanToken scanTok;

    // Token start/end (Ragel naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are Ragel naming
    // - buf is our own naming

    const char* buf = &(str[strBeg]);
    const char* eof = &(str[strEnd]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    
#line 407 "fieldExprScanner.cc"
	{
	cs = fieldExpr_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 535 "fieldExprScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 419 "fieldExprScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 189 "fieldExprScanner.rl"
	{te = p+1;{
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr4:
#line 189 "fieldExprScanner.rl"
	{te = p+1;{
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr5:
#line 167 "fieldExprScanner.rl"
	{{p = ((te))-1;}{
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        if (readScalar(std::string(ts, te-ts), scanTok.svalue))
        {
            parser_->parse(TOKEN_OF(NUMBER), &scanTok);
        }
        else
        {
            driver_.reportFatal
            (
                "Error parsing number: " + std::string(ts, te-ts)
            );
        }

        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr8:
#line 232 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(EQUAL); }}
	goto st11;
tr9:
#line 284 "fieldExprScanner.rl"
	{{p = ((te))-1;}{ EMIT_TOKEN(TENSOR); }}
	goto st11;
tr11:
#line 292 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(IDENTITY_TENSOR); }}
	goto st11;
tr12:
#line 235 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LOR); }}
	goto st11;
tr16:
#line 217 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PERCENT); }}
	goto st11;
tr19:
#line 218 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LPAREN); }}
	goto st11;
tr20:
#line 219 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(RPAREN); }}
	goto st11;
tr21:
#line 220 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(TIMES); }}
	goto st11;
tr22:
#line 221 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(PLUS); }}
	goto st11;
tr23:
#line 223 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COMMA); }}
	goto st11;
tr24:
#line 222 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(MINUS); }}
	goto st11;
tr26:
#line 225 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(DIVIDE); }}
	goto st11;
tr28:
#line 227 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(COLON); }}
	goto st11;
tr32:
#line 226 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(QUESTION); }}
	goto st11;
tr35:
#line 238 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(BIT_XOR); }}
	goto st11;
tr51:
#line 211 "fieldExprScanner.rl"
	{te = p;p--;}
	goto st11;
tr52:
#line 216 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NOT); }}
	goto st11;
tr53:
#line 233 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(NOT_EQUAL); }}
	goto st11;
tr54:
#line 236 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(BIT_AND); }}
	goto st11;
tr55:
#line 234 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LAND); }}
	goto st11;
tr56:
#line 224 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(DOT); }}
	goto st11;
tr59:
#line 167 "fieldExprScanner.rl"
	{te = p;p--;{
        driver_.parsePosition() = (ts-buf);

        DebugInfo
            << "Number:" << std::string(ts, te-ts).c_str()
            << " at " << driver_.parsePosition() << nl;

        if (readScalar(std::string(ts, te-ts), scanTok.svalue))
        {
            parser_->parse(TOKEN_OF(NUMBER), &scanTok);
        }
        else
        {
            driver_.reportFatal
            (
                "Error parsing number: " + std::string(ts, te-ts)
            );
        }

        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr61:
#line 195 "fieldExprScanner.rl"
	{te = p;p--;{
        // Tokenized ".method" - dispatch '.' and "method" separately
        driver_.parsePosition() = (ts-buf);
        dispatch_method(driver_, scanTok, word(ts+1, te-ts-1, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr62:
#line 228 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LESS); }}
	goto st11;
tr63:
#line 229 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(LESS_EQ); }}
	goto st11;
tr64:
#line 230 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(GREATER); }}
	goto st11;
tr65:
#line 231 "fieldExprScanner.rl"
	{te = p+1;{ EMIT_TOKEN(GREATER_EQ); }}
	goto st11;
tr66:
#line 189 "fieldExprScanner.rl"
	{te = p;p--;{
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }}
	goto st11;
tr68:
#line 1 "NONE"
	{	switch( act ) {
	case 26:
	{{p = ((te))-1;} EMIT_TOKEN(PI); }
	break;
	case 27:
	{{p = ((te))-1;} EMIT_TOKEN(DEG_TO_RAD); }
	break;
	case 28:
	{{p = ((te))-1;} EMIT_TOKEN(RAD_TO_DEG); }
	break;
	case 29:
	{{p = ((te))-1;} EMIT_TOKEN(EXP); }
	break;
	case 31:
	{{p = ((te))-1;} EMIT_TOKEN(LOG10); }
	break;
	case 32:
	{{p = ((te))-1;} EMIT_TOKEN(POW); }
	break;
	case 34:
	{{p = ((te))-1;} EMIT_TOKEN(SQRT); }
	break;
	case 35:
	{{p = ((te))-1;} EMIT_TOKEN(CBRT); }
	break;
	case 39:
	{{p = ((te))-1;} EMIT_TOKEN(ASIN); }
	break;
	case 40:
	{{p = ((te))-1;} EMIT_TOKEN(ACOS); }
	break;
	case 42:
	{{p = ((te))-1;} EMIT_TOKEN(ATAN2); }
	break;
	case 43:
	{{p = ((te))-1;} EMIT_TOKEN(SINH); }
	break;
	case 44:
	{{p = ((te))-1;} EMIT_TOKEN(COSH); }
	break;
	case 45:
	{{p = ((te))-1;} EMIT_TOKEN(TANH); }
	break;
	case 47:
	{{p = ((te))-1;} EMIT_TOKEN(MAGSQR); }
	break;
	case 50:
	{{p = ((te))-1;} EMIT_TOKEN(POS0); }
	break;
	case 51:
	{{p = ((te))-1;} EMIT_TOKEN(NEG0); }
	break;
	case 52:
	{{p = ((te))-1;} EMIT_TOKEN(SIGN); }
	break;
	case 53:
	{{p = ((te))-1;} EMIT_TOKEN(MIN); }
	break;
	case 54:
	{{p = ((te))-1;} EMIT_TOKEN(MAX); }
	break;
	case 55:
	{{p = ((te))-1;} EMIT_TOKEN(AVERAGE); }
	break;
	case 56:
	{{p = ((te))-1;} EMIT_TOKEN(SUM); }
	break;
	case 57:
	{{p = ((te))-1;} EMIT_TOKEN(RAND); }
	break;
	case 58:
	{{p = ((te))-1;} EMIT_TOKEN(BOOL); }
	break;
	case 59:
	{{p = ((te))-1;} EMIT_TOKEN(VECTOR); }
	break;
	case 61:
	{{p = ((te))-1;} EMIT_TOKEN(SYM_TENSOR); }
	break;
	case 62:
	{{p = ((te))-1;} EMIT_TOKEN(SPH_TENSOR); }
	break;
	case 63:
	{{p = ((te))-1;} EMIT_TOKEN(ZERO); }
	break;
	case 64:
	{{p = ((te))-1;} EMIT_TOKEN(LTRUE); }
	break;
	case 65:
	{{p = ((te))-1;} EMIT_TOKEN(LFALSE); }
	break;
	case 67:
	{{p = ((te))-1;} EMIT_TOKEN(ARG); }
	break;
	case 68:
	{{p = ((te))-1;}
        driver_.parsePosition() = (ts-buf);
        dispatch_ident(driver_, scanTok, word(ts, te-ts, false));
        driver_.parsePosition() = (p-buf);
    }
	break;
	}
	}
	goto st11;
tr84:
#line 260 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(ATAN); }}
	goto st11;
tr99:
#line 256 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(COS); }}
	goto st11;
tr116:
#line 249 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(LOG); }}
	goto st11;
tr123:
#line 265 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(MAG); }}
	goto st11;
tr130:
#line 269 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(NEG); }}
	goto st11;
tr136:
#line 268 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(POS); }}
	goto st11;
tr155:
#line 255 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SIN); }}
	goto st11;
tr171:
#line 252 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(SQR); }}
	goto st11;
tr186:
#line 257 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TAN); }}
	goto st11;
tr192:
#line 284 "fieldExprScanner.rl"
	{te = p;p--;{ EMIT_TOKEN(TENSOR); }}
	goto st11;
st11:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof11;
case 11:
#line 1 "NONE"
	{ts = p;}
#line 760 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 32: goto st12;
		case 33: goto st13;
		case 34: goto st1;
		case 37: goto tr16;
		case 38: goto st14;
		case 39: goto st3;
		case 40: goto tr19;
		case 41: goto tr20;
		case 42: goto tr21;
		case 43: goto tr22;
		case 44: goto tr23;
		case 45: goto tr24;
		case 46: goto st15;
		case 47: goto tr26;
		case 58: goto tr28;
		case 60: goto st20;
		case 61: goto st7;
		case 62: goto st21;
		case 63: goto tr32;
		case 90: goto st24;
		case 94: goto tr35;
		case 95: goto st22;
		case 97: goto st27;
		case 98: goto st41;
		case 99: goto st44;
		case 100: goto st49;
		case 101: goto st56;
		case 102: goto st58;
		case 108: goto st62;
		case 109: goto st66;
		case 110: goto st72;
		case 112: goto st75;
		case 114: goto st78;
		case 115: goto st86;
		case 116: goto st114;
		case 118: goto st124;
		case 124: goto st10;
	}
	if ( (*p) < 48 ) {
		if ( 9 <= (*p) && (*p) <= 13 )
			goto st12;
	} else if ( (*p) > 57 ) {
		if ( (*p) > 89 ) {
			if ( 103 <= (*p) && (*p) <= 122 )
				goto st22;
		} else if ( (*p) >= 65 )
			goto st22;
	} else
		goto tr27;
	goto st0;
st0:
cs = 0;
	goto _out;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	if ( (*p) == 32 )
		goto st12;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st12;
	goto tr51;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 61 )
		goto tr53;
	goto tr52;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( (*p) == 34 )
		goto st0;
	goto st2;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	if ( (*p) == 34 )
		goto tr2;
	goto st2;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 38 )
		goto tr55;
	goto tr54;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( (*p) == 39 )
		goto st0;
	goto st4;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 39 )
		goto tr4;
	goto st4;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr57;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st18;
	} else
		goto st18;
	goto tr56;
tr57:
#line 1 "NONE"
	{te = p+1;}
	goto st16;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
#line 887 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr57;
	goto tr59;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 43: goto st6;
		case 45: goto st6;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr5;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr5;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st17;
	goto tr59;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto st18;
	} else if ( (*p) >= 65 )
		goto st18;
	goto tr61;
tr27:
#line 1 "NONE"
	{te = p+1;}
	goto st19;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
#line 938 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr57;
		case 69: goto st5;
		case 101: goto st5;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr27;
	goto tr59;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 61 )
		goto tr63;
	goto tr62;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 61 )
		goto tr8;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	if ( (*p) == 61 )
		goto tr65;
	goto tr64;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
tr67:
#line 1 "NONE"
	{te = p+1;}
#line 189 "fieldExprScanner.rl"
	{act = 68;}
	goto st23;
tr71:
#line 1 "NONE"
	{te = p+1;}
#line 289 "fieldExprScanner.rl"
	{act = 63;}
	goto st23;
tr78:
#line 1 "NONE"
	{te = p+1;}
#line 259 "fieldExprScanner.rl"
	{act = 40;}
	goto st23;
tr79:
#line 1 "NONE"
	{te = p+1;}
#line 293 "fieldExprScanner.rl"
	{act = 67;}
	goto st23;
tr81:
#line 1 "NONE"
	{te = p+1;}
#line 258 "fieldExprScanner.rl"
	{act = 39;}
	goto st23;
tr85:
#line 1 "NONE"
	{te = p+1;}
#line 261 "fieldExprScanner.rl"
	{act = 42;}
	goto st23;
tr90:
#line 1 "NONE"
	{te = p+1;}
#line 277 "fieldExprScanner.rl"
	{act = 55;}
	goto st23;
tr93:
#line 1 "NONE"
	{te = p+1;}
#line 282 "fieldExprScanner.rl"
	{act = 58;}
	goto st23;
tr97:
#line 1 "NONE"
	{te = p+1;}
#line 254 "fieldExprScanner.rl"
	{act = 35;}
	goto st23;
tr100:
#line 1 "NONE"
	{te = p+1;}
#line 263 "fieldExprScanner.rl"
	{act = 44;}
	goto st23;
tr107:
#line 1 "NONE"
	{te = p+1;}
#line 246 "fieldExprScanner.rl"
	{act = 27;}
	goto st23;
tr109:
#line 1 "NONE"
	{te = p+1;}
#line 248 "fieldExprScanner.rl"
	{act = 29;}
	goto st23;
tr113:
#line 1 "NONE"
	{te = p+1;}
#line 291 "fieldExprScanner.rl"
	{act = 65;}
	goto st23;
tr118:
#line 1 "NONE"
	{te = p+1;}
#line 250 "fieldExprScanner.rl"
	{act = 31;}
	goto st23;
tr122:
#line 1 "NONE"
	{te = p+1;}
#line 276 "fieldExprScanner.rl"
	{act = 54;}
	goto st23;
tr126:
#line 1 "NONE"
	{te = p+1;}
#line 266 "fieldExprScanner.rl"
	{act = 47;}
	goto st23;
tr127:
#line 1 "NONE"
	{te = p+1;}
#line 275 "fieldExprScanner.rl"
	{act = 53;}
	goto st23;
tr131:
#line 1 "NONE"
	{te = p+1;}
#line 271 "fieldExprScanner.rl"
	{act = 51;}
	goto st23;
tr132:
#line 1 "NONE"
	{te = p+1;}
#line 245 "fieldExprScanner.rl"
	{act = 26;}
	goto st23;
tr135:
#line 1 "NONE"
	{te = p+1;}
#line 251 "fieldExprScanner.rl"
	{act = 32;}
	goto st23;
tr137:
#line 1 "NONE"
	{te = p+1;}
#line 270 "fieldExprScanner.rl"
	{act = 50;}
	goto st23;
tr145:
#line 1 "NONE"
	{te = p+1;}
#line 247 "fieldExprScanner.rl"
	{act = 28;}
	goto st23;
tr146:
#line 1 "NONE"
	{te = p+1;}
#line 279 "fieldExprScanner.rl"
	{act = 57;}
	goto st23;
tr154:
#line 1 "NONE"
	{te = p+1;}
#line 272 "fieldExprScanner.rl"
	{act = 52;}
	goto st23;
tr156:
#line 1 "NONE"
	{te = p+1;}
#line 262 "fieldExprScanner.rl"
	{act = 43;}
	goto st23;
tr169:
#line 1 "NONE"
	{te = p+1;}
#line 286 "fieldExprScanner.rl"
	{act = 62;}
	goto st23;
tr172:
#line 1 "NONE"
	{te = p+1;}
#line 253 "fieldExprScanner.rl"
	{act = 34;}
	goto st23;
tr173:
#line 1 "NONE"
	{te = p+1;}
#line 278 "fieldExprScanner.rl"
	{act = 56;}
	goto st23;
tr181:
#line 1 "NONE"
	{te = p+1;}
#line 285 "fieldExprScanner.rl"
	{act = 61;}
	goto st23;
tr187:
#line 1 "NONE"
	{te = p+1;}
#line 264 "fieldExprScanner.rl"
	{act = 45;}
	goto st23;
tr195:
#line 1 "NONE"
	{te = p+1;}
#line 290 "fieldExprScanner.rl"
	{act = 64;}
	goto st23;
tr200:
#line 1 "NONE"
	{te = p+1;}
#line 283 "fieldExprScanner.rl"
	{act = 59;}
	goto st23;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
#line 1181 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr68;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st25;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto st26;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto tr71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 99: goto st28;
		case 114: goto st30;
		case 115: goto st31;
		case 116: goto st33;
		case 118: goto st36;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st29;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto tr78;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto tr79;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 105: goto st32;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto tr81;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st34;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto st35;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	switch( (*p) ) {
		case 46: goto tr67;
		case 50: goto tr85;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr84;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st37;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto st38;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st39;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto st40;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto tr90;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st42;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st43;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 108: goto tr93;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 98: goto st45;
		case 111: goto st47;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto st46;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 116: goto tr97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto st48;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 104: goto tr100;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr99;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st50;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto st51;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 46: goto tr67;
		case 84: goto st52;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st53;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	switch( (*p) ) {
		case 46: goto tr67;
		case 82: goto st54;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st55;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 100: goto tr107;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 120: goto st57;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 112: goto tr109;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st59;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 108: goto st60;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto st61;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto tr113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st63;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto st64;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	switch( (*p) ) {
		case 46: goto tr67;
		case 49: goto st65;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr116;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	switch( (*p) ) {
		case 46: goto tr67;
		case 48: goto tr118;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st67;
		case 105: goto st71;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto st68;
		case 120: goto tr122;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	switch( (*p) ) {
		case 46: goto tr67;
		case 83: goto st69;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr123;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 113: goto st70;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto tr126;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto tr127;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st73;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto st74;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 46: goto tr67;
		case 48: goto tr131;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr130;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 105: goto tr132;
		case 111: goto st76;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto st77;
		case 119: goto tr135;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	switch( (*p) ) {
		case 46: goto tr67;
		case 48: goto tr137;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 49 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr136;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st79;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 100: goto st80;
		case 110: goto st85;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 46: goto tr67;
		case 84: goto st81;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st82;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	switch( (*p) ) {
		case 46: goto tr67;
		case 68: goto st83;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st84;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto tr145;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 100: goto tr146;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st86:
	if ( ++p == pe )
		goto _test_eof86;
case 86:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 105: goto st87;
		case 112: goto st90;
		case 113: goto st103;
		case 117: goto st105;
		case 121: goto st106;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 103: goto st88;
		case 110: goto st89;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto tr154;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 104: goto tr156;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr155;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 104: goto st91;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st92;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto st93;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 105: goto st94;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st94:
	if ( ++p == pe )
		goto _test_eof94;
case 94:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 99: goto st95;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st95:
	if ( ++p == pe )
		goto _test_eof95;
case 95:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st96;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st96:
	if ( ++p == pe )
		goto _test_eof96;
case 96:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 108: goto st97;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st97:
	if ( ++p == pe )
		goto _test_eof97;
case 97:
	switch( (*p) ) {
		case 46: goto tr67;
		case 84: goto st98;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st98:
	if ( ++p == pe )
		goto _test_eof98;
case 98:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st99;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st99:
	if ( ++p == pe )
		goto _test_eof99;
case 99:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto st100;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st100:
	if ( ++p == pe )
		goto _test_eof100;
case 100:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto st101;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st101:
	if ( ++p == pe )
		goto _test_eof101;
case 101:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st102;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st102:
	if ( ++p == pe )
		goto _test_eof102;
case 102:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto tr169;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st103:
	if ( ++p == pe )
		goto _test_eof103;
case 103:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto st104;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st104:
	if ( ++p == pe )
		goto _test_eof104;
case 104:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 116: goto tr172;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr171;
st105:
	if ( ++p == pe )
		goto _test_eof105;
case 105:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 109: goto tr173;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st106:
	if ( ++p == pe )
		goto _test_eof106;
case 106:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 109: goto st107;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st107:
	if ( ++p == pe )
		goto _test_eof107;
case 107:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 109: goto st108;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st108:
	if ( ++p == pe )
		goto _test_eof108;
case 108:
	switch( (*p) ) {
		case 46: goto tr67;
		case 84: goto st109;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st109:
	if ( ++p == pe )
		goto _test_eof109;
case 109:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st110;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st110:
	if ( ++p == pe )
		goto _test_eof110;
case 110:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto st111;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st111:
	if ( ++p == pe )
		goto _test_eof111;
case 111:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto st112;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st112:
	if ( ++p == pe )
		goto _test_eof112;
case 112:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st113;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st113:
	if ( ++p == pe )
		goto _test_eof113;
case 113:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto tr181;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st114:
	if ( ++p == pe )
		goto _test_eof114;
case 114:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 97: goto st115;
		case 101: goto st117;
		case 114: goto st122;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 98 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st115:
	if ( ++p == pe )
		goto _test_eof115;
case 115:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto st116;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st116:
	if ( ++p == pe )
		goto _test_eof116;
case 116:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 104: goto tr187;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr186;
st117:
	if ( ++p == pe )
		goto _test_eof117;
case 117:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 110: goto st118;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st118:
	if ( ++p == pe )
		goto _test_eof118;
case 118:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 115: goto st119;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st119:
	if ( ++p == pe )
		goto _test_eof119;
case 119:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st120;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st120:
	if ( ++p == pe )
		goto _test_eof120;
case 120:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto tr191;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
tr191:
#line 1 "NONE"
	{te = p+1;}
	goto st121;
st121:
	if ( ++p == pe )
		goto _test_eof121;
case 121:
#line 2966 "fieldExprScanner.cc"
	switch( (*p) ) {
		case 46: goto tr67;
		case 58: goto st8;
		case 95: goto tr67;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr192;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 58 )
		goto st9;
	goto tr9;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	if ( (*p) == 73 )
		goto tr11;
	goto tr9;
st122:
	if ( ++p == pe )
		goto _test_eof122;
case 122:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 117: goto st123;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st123:
	if ( ++p == pe )
		goto _test_eof123;
case 123:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto tr195;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st124:
	if ( ++p == pe )
		goto _test_eof124;
case 124:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 101: goto st125;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st125:
	if ( ++p == pe )
		goto _test_eof125;
case 125:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 99: goto st126;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st126:
	if ( ++p == pe )
		goto _test_eof126;
case 126:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 116: goto st127;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st127:
	if ( ++p == pe )
		goto _test_eof127;
case 127:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 111: goto st128;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st128:
	if ( ++p == pe )
		goto _test_eof128;
case 128:
	switch( (*p) ) {
		case 46: goto tr67;
		case 95: goto tr67;
		case 114: goto tr200;
	}
	if ( (*p) < 65 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr67;
	} else if ( (*p) > 90 ) {
		if ( 97 <= (*p) && (*p) <= 122 )
			goto tr67;
	} else
		goto tr67;
	goto tr66;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 124 )
		goto tr12;
	goto st0;
	}
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof21: cs = 21; goto _test_eof; 
	_test_eof22: cs = 22; goto _test_eof; 
	_test_eof23: cs = 23; goto _test_eof; 
	_test_eof24: cs = 24; goto _test_eof; 
	_test_eof25: cs = 25; goto _test_eof; 
	_test_eof26: cs = 26; goto _test_eof; 
	_test_eof27: cs = 27; goto _test_eof; 
	_test_eof28: cs = 28; goto _test_eof; 
	_test_eof29: cs = 29; goto _test_eof; 
	_test_eof30: cs = 30; goto _test_eof; 
	_test_eof31: cs = 31; goto _test_eof; 
	_test_eof32: cs = 32; goto _test_eof; 
	_test_eof33: cs = 33; goto _test_eof; 
	_test_eof34: cs = 34; goto _test_eof; 
	_test_eof35: cs = 35; goto _test_eof; 
	_test_eof36: cs = 36; goto _test_eof; 
	_test_eof37: cs = 37; goto _test_eof; 
	_test_eof38: cs = 38; goto _test_eof; 
	_test_eof39: cs = 39; goto _test_eof; 
	_test_eof40: cs = 40; goto _test_eof; 
	_test_eof41: cs = 41; goto _test_eof; 
	_test_eof42: cs = 42; goto _test_eof; 
	_test_eof43: cs = 43; goto _test_eof; 
	_test_eof44: cs = 44; goto _test_eof; 
	_test_eof45: cs = 45; goto _test_eof; 
	_test_eof46: cs = 46; goto _test_eof; 
	_test_eof47: cs = 47; goto _test_eof; 
	_test_eof48: cs = 48; goto _test_eof; 
	_test_eof49: cs = 49; goto _test_eof; 
	_test_eof50: cs = 50; goto _test_eof; 
	_test_eof51: cs = 51; goto _test_eof; 
	_test_eof52: cs = 52; goto _test_eof; 
	_test_eof53: cs = 53; goto _test_eof; 
	_test_eof54: cs = 54; goto _test_eof; 
	_test_eof55: cs = 55; goto _test_eof; 
	_test_eof56: cs = 56; goto _test_eof; 
	_test_eof57: cs = 57; goto _test_eof; 
	_test_eof58: cs = 58; goto _test_eof; 
	_test_eof59: cs = 59; goto _test_eof; 
	_test_eof60: cs = 60; goto _test_eof; 
	_test_eof61: cs = 61; goto _test_eof; 
	_test_eof62: cs = 62; goto _test_eof; 
	_test_eof63: cs = 63; goto _test_eof; 
	_test_eof64: cs = 64; goto _test_eof; 
	_test_eof65: cs = 65; goto _test_eof; 
	_test_eof66: cs = 66; goto _test_eof; 
	_test_eof67: cs = 67; goto _test_eof; 
	_test_eof68: cs = 68; goto _test_eof; 
	_test_eof69: cs = 69; goto _test_eof; 
	_test_eof70: cs = 70; goto _test_eof; 
	_test_eof71: cs = 71; goto _test_eof; 
	_test_eof72: cs = 72; goto _test_eof; 
	_test_eof73: cs = 73; goto _test_eof; 
	_test_eof74: cs = 74; goto _test_eof; 
	_test_eof75: cs = 75; goto _test_eof; 
	_test_eof76: cs = 76; goto _test_eof; 
	_test_eof77: cs = 77; goto _test_eof; 
	_test_eof78: cs = 78; goto _test_eof; 
	_test_eof79: cs = 79; goto _test_eof; 
	_test_eof80: cs = 80; goto _test_eof; 
	_test_eof81: cs = 81; goto _test_eof; 
	_test_eof82: cs = 82; goto _test_eof; 
	_test_eof83: cs = 83; goto _test_eof; 
	_test_eof84: cs = 84; goto _test_eof; 
	_test_eof85: cs = 85; goto _test_eof; 
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof94: cs = 94; goto _test_eof; 
	_test_eof95: cs = 95; goto _test_eof; 
	_test_eof96: cs = 96; goto _test_eof; 
	_test_eof97: cs = 97; goto _test_eof; 
	_test_eof98: cs = 98; goto _test_eof; 
	_test_eof99: cs = 99; goto _test_eof; 
	_test_eof100: cs = 100; goto _test_eof; 
	_test_eof101: cs = 101; goto _test_eof; 
	_test_eof102: cs = 102; goto _test_eof; 
	_test_eof103: cs = 103; goto _test_eof; 
	_test_eof104: cs = 104; goto _test_eof; 
	_test_eof105: cs = 105; goto _test_eof; 
	_test_eof106: cs = 106; goto _test_eof; 
	_test_eof107: cs = 107; goto _test_eof; 
	_test_eof108: cs = 108; goto _test_eof; 
	_test_eof109: cs = 109; goto _test_eof; 
	_test_eof110: cs = 110; goto _test_eof; 
	_test_eof111: cs = 111; goto _test_eof; 
	_test_eof112: cs = 112; goto _test_eof; 
	_test_eof113: cs = 113; goto _test_eof; 
	_test_eof114: cs = 114; goto _test_eof; 
	_test_eof115: cs = 115; goto _test_eof; 
	_test_eof116: cs = 116; goto _test_eof; 
	_test_eof117: cs = 117; goto _test_eof; 
	_test_eof118: cs = 118; goto _test_eof; 
	_test_eof119: cs = 119; goto _test_eof; 
	_test_eof120: cs = 120; goto _test_eof; 
	_test_eof121: cs = 121; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof122: cs = 122; goto _test_eof; 
	_test_eof123: cs = 123; goto _test_eof; 
	_test_eof124: cs = 124; goto _test_eof; 
	_test_eof125: cs = 125; goto _test_eof; 
	_test_eof126: cs = 126; goto _test_eof; 
	_test_eof127: cs = 127; goto _test_eof; 
	_test_eof128: cs = 128; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 12: goto tr51;
	case 13: goto tr52;
	case 14: goto tr54;
	case 15: goto tr56;
	case 16: goto tr59;
	case 5: goto tr5;
	case 6: goto tr5;
	case 17: goto tr59;
	case 18: goto tr61;
	case 19: goto tr59;
	case 20: goto tr62;
	case 21: goto tr64;
	case 22: goto tr66;
	case 23: goto tr68;
	case 24: goto tr66;
	case 25: goto tr66;
	case 26: goto tr66;
	case 27: goto tr66;
	case 28: goto tr66;
	case 29: goto tr66;
	case 30: goto tr66;
	case 31: goto tr66;
	case 32: goto tr66;
	case 33: goto tr66;
	case 34: goto tr66;
	case 35: goto tr84;
	case 36: goto tr66;
	case 37: goto tr66;
	case 38: goto tr66;
	case 39: goto tr66;
	case 40: goto tr66;
	case 41: goto tr66;
	case 42: goto tr66;
	case 43: goto tr66;
	case 44: goto tr66;
	case 45: goto tr66;
	case 46: goto tr66;
	case 47: goto tr66;
	case 48: goto tr99;
	case 49: goto tr66;
	case 50: goto tr66;
	case 51: goto tr66;
	case 52: goto tr66;
	case 53: goto tr66;
	case 54: goto tr66;
	case 55: goto tr66;
	case 56: goto tr66;
	case 57: goto tr66;
	case 58: goto tr66;
	case 59: goto tr66;
	case 60: goto tr66;
	case 61: goto tr66;
	case 62: goto tr66;
	case 63: goto tr66;
	case 64: goto tr116;
	case 65: goto tr66;
	case 66: goto tr66;
	case 67: goto tr66;
	case 68: goto tr123;
	case 69: goto tr66;
	case 70: goto tr66;
	case 71: goto tr66;
	case 72: goto tr66;
	case 73: goto tr66;
	case 74: goto tr130;
	case 75: goto tr66;
	case 76: goto tr66;
	case 77: goto tr136;
	case 78: goto tr66;
	case 79: goto tr66;
	case 80: goto tr66;
	case 81: goto tr66;
	case 82: goto tr66;
	case 83: goto tr66;
	case 84: goto tr66;
	case 85: goto tr66;
	case 86: goto tr66;
	case 87: goto tr66;
	case 88: goto tr66;
	case 89: goto tr155;
	case 90: goto tr66;
	case 91: goto tr66;
	case 92: goto tr66;
	case 93: goto tr66;
	case 94: goto tr66;
	case 95: goto tr66;
	case 96: goto tr66;
	case 97: goto tr66;
	case 98: goto tr66;
	case 99: goto tr66;
	case 100: goto tr66;
	case 101: goto tr66;
	case 102: goto tr66;
	case 103: goto tr66;
	case 104: goto tr171;
	case 105: goto tr66;
	case 106: goto tr66;
	case 107: goto tr66;
	case 108: goto tr66;
	case 109: goto tr66;
	case 110: goto tr66;
	case 111: goto tr66;
	case 112: goto tr66;
	case 113: goto tr66;
	case 114: goto tr66;
	case 115: goto tr66;
	case 116: goto tr186;
	case 117: goto tr66;
	case 118: goto tr66;
	case 119: goto tr66;
	case 120: goto tr66;
	case 121: goto tr192;
	case 8: goto tr9;
	case 9: goto tr9;
	case 122: goto tr66;
	case 123: goto tr66;
	case 124: goto tr66;
	case 125: goto tr66;
	case 126: goto tr66;
	case 127: goto tr66;
	case 128: goto tr66;
	}
	}

	_out: {}
	}

#line 537 "fieldExprScanner.rl"
   /* ^^^ FSM execution here ^^^ */;

    if (0 == cs)
    {
        driver_.reportFatal("Parse error while scanning", (p-buf));
    }

    if (p != eof)
    {
        driver_.reportFatal("Parsing failed with remaining content", (p-buf));
    }

    // Terminate parser execution
    parser_->parse(0, nullptr);
    parser_->stop();

    if (debug & 0x6)
    {
        InfoErr<< "Done parse." << nl;
    }

    // Restore debug value
    debug = oldDebug;

    return true;
}


// ************************************************************************* //
