
#line 1 "evalStringToScalarScanner.rl"
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


#line 61 "evalStringToScalarScanner.cc"
static const int evalScanner_start = 86;
static const int evalScanner_first_final = 86;
static const int evalScanner_error = 0;

static const int evalScanner_en_main = 86;


#line 60 "evalStringToScalarScanner.rl"


#define TOKEN_OF(T)         TOK_##T
#define EMIT_TOKEN(T)                                                         \
    /* Inform driver of last position */                                      \
    driver.parsePosition() = (p-buf);                                         \
    DebugInfo<< "TOKEN_" #T << " at " << driver.parsePosition() << nl;        \
    parser_->parse(TOKEN_OF(T), 0)



#line 163 "evalStringToScalarScanner.rl"



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::parsing::evalStringToScalar::scanner::~scanner()
{
    if (parser_)
    {
        delete parser_;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::parsing::evalStringToScalar::scanner::process
(
    const std::string& str,
    parseDriver& driver
)
{
    if (!parser_)
    {
        parser_ = new parser();
    }

    driver.content(str);
    parser_->start(driver);

    // Ragel token start/end (required naming)
    const char* ts;
    const char* te;

    // Local buffer data.
    // - p, pe, eof are required Ragel naming
    // - buf is our own naming

    const char* buf = &(str[0]);
    const char* eof = &(str[str.length()]);
    const char* p = buf;
    const char* pe = eof;

    // Initialize FSM variables
    
#line 127 "evalStringToScalarScanner.cc"
	{
	cs = evalScanner_start;
	ts = 0;
	te = 0;
	act = 0;
	}

#line 207 "evalStringToScalarScanner.rl"
   /* ^^^ FSM initialization here ^^^ */;

    
#line 139 "evalStringToScalarScanner.cc"
	{
	if ( p == pe )
		goto _test_eof;
	switch ( cs )
	{
tr2:
#line 86 "evalStringToScalarScanner.rl"
	{{p = ((te))-1;}{
        // Inform driver of last position
        driver.parsePosition() = (p-buf);

        DebugInfo
            << "NUMBER:" << std::string(ts, te-ts).c_str()
            << " at " << driver.parsePosition() << nl;

        scalar val;

        if (readScalar(std::string(ts, te-ts), val))
        {
            // Emit number
            parser_->parse(TOKEN_OF(NUMBER), val);
        }
        else
        {
            driver.reportFatal("Error reading scalar value");
        }
    }}
	goto st86;
tr5:
#line 157 "evalStringToScalarScanner.rl"
	{{p = ((te))-1;}}
	goto st86;
tr8:
#line 74 "evalStringToScalarScanner.rl"
	{
        // Inform driver of last position
        driver.parsePosition() = (p-buf);
    }
#line 159 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr14:
#line 127 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(ACOS); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr17:
#line 126 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(ASIN); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr21:
#line 128 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(ATAN); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr23:
#line 129 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(ATAN2); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr28:
#line 122 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(CBRT); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr31:
#line 124 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(COS); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr33:
#line 132 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(COSH); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr41:
#line 114 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(DEG_TO_RAD); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr44:
#line 116 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(EXP); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr49:
#line 130 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(HYPOT); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr53:
#line 117 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(LOG); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr56:
#line 118 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(LOG10); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr62:
#line 136 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(MAG); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr66:
#line 137 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(MAGSQR); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr67:
#line 135 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(MAX); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr69:
#line 134 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(MIN); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr72:
#line 113 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(PI); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr74:
#line 119 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(POW); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr83:
#line 115 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(RAD_TO_DEG); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr85:
#line 138 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(RAND); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr90:
#line 123 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(SIN); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr92:
#line 131 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(SINH); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr95:
#line 120 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(SQR); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr97:
#line 121 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(SQRT); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr101:
#line 125 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(TAN); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr103:
#line 133 "evalStringToScalarScanner.rl"
	{ p--; EMIT_TOKEN(TANH); }
#line 156 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr105:
#line 142 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(LPAREN); }
#line 157 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr106:
#line 143 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(RPAREN); }
#line 157 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr107:
#line 146 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(TIMES); }
#line 157 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr108:
#line 144 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(PLUS); }
#line 157 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr109:
#line 148 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(COMMA); }
#line 157 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr110:
#line 145 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(MINUS); }
#line 157 "evalStringToScalarScanner.rl"
	{te = p+1;}
	goto st86;
tr125:
#line 153 "evalStringToScalarScanner.rl"
	{te = p;p--;}
	goto st86;
tr126:
#line 86 "evalStringToScalarScanner.rl"
	{te = p;p--;{
        // Inform driver of last position
        driver.parsePosition() = (p-buf);

        DebugInfo
            << "NUMBER:" << std::string(ts, te-ts).c_str()
            << " at " << driver.parsePosition() << nl;

        scalar val;

        if (readScalar(std::string(ts, te-ts), val))
        {
            // Emit number
            parser_->parse(TOKEN_OF(NUMBER), val);
        }
        else
        {
            driver.reportFatal("Error reading scalar value");
        }
    }}
	goto st86;
tr128:
#line 157 "evalStringToScalarScanner.rl"
	{te = p;p--;}
	goto st86;
tr130:
#line 160 "evalStringToScalarScanner.rl"
	{te = p;p--;}
	goto st86;
st86:
#line 1 "NONE"
	{ts = 0;}
	if ( ++p == pe )
		goto _test_eof86;
case 86:
#line 1 "NONE"
	{ts = p;}
#line 416 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 32: goto st87;
		case 40: goto tr105;
		case 41: goto tr106;
		case 42: goto tr107;
		case 43: goto tr108;
		case 44: goto tr109;
		case 45: goto tr110;
		case 46: goto st1;
		case 47: goto tr112;
		case 97: goto st6;
		case 99: goto st18;
		case 100: goto st26;
		case 101: goto st34;
		case 104: goto st37;
		case 108: goto st42;
		case 109: goto st48;
		case 112: goto st58;
		case 114: goto st62;
		case 115: goto st72;
		case 116: goto st81;
	}
	if ( (*p) > 13 ) {
		if ( 48 <= (*p) && (*p) <= 57 )
			goto tr113;
	} else if ( (*p) >= 9 )
		goto st87;
	goto st0;
st0:
cs = 0;
	goto _out;
st87:
	if ( ++p == pe )
		goto _test_eof87;
case 87:
	if ( (*p) == 32 )
		goto st87;
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st87;
	goto tr125;
st1:
	if ( ++p == pe )
		goto _test_eof1;
case 1:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr0;
	goto st0;
tr0:
#line 1 "NONE"
	{te = p+1;}
	goto st88;
st88:
	if ( ++p == pe )
		goto _test_eof88;
case 88:
#line 472 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 69: goto st2;
		case 101: goto st2;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr0;
	goto tr126;
st2:
	if ( ++p == pe )
		goto _test_eof2;
case 2:
	switch( (*p) ) {
		case 43: goto st3;
		case 45: goto st3;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st89;
	goto tr2;
st3:
	if ( ++p == pe )
		goto _test_eof3;
case 3:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st89;
	goto tr2;
st89:
	if ( ++p == pe )
		goto _test_eof89;
case 89:
	if ( 48 <= (*p) && (*p) <= 57 )
		goto st89;
	goto tr126;
tr112:
#line 1 "NONE"
	{te = p+1;}
#line 147 "evalStringToScalarScanner.rl"
	{ EMIT_TOKEN(DIVIDE); }
	goto st90;
st90:
	if ( ++p == pe )
		goto _test_eof90;
case 90:
#line 515 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 42: goto st4;
		case 47: goto st91;
	}
	goto tr128;
st4:
	if ( ++p == pe )
		goto _test_eof4;
case 4:
	if ( (*p) == 42 )
		goto st5;
	goto st4;
st5:
	if ( ++p == pe )
		goto _test_eof5;
case 5:
	switch( (*p) ) {
		case 42: goto st5;
		case 47: goto tr8;
	}
	goto st4;
st91:
	if ( ++p == pe )
		goto _test_eof91;
case 91:
	if ( (*p) == 10 )
		goto tr131;
	goto st91;
tr131:
#line 74 "evalStringToScalarScanner.rl"
	{
        // Inform driver of last position
        driver.parsePosition() = (p-buf);
    }
	goto st92;
st92:
	if ( ++p == pe )
		goto _test_eof92;
case 92:
#line 555 "evalStringToScalarScanner.cc"
	if ( (*p) == 10 )
		goto tr131;
	goto tr130;
tr113:
#line 1 "NONE"
	{te = p+1;}
	goto st93;
st93:
	if ( ++p == pe )
		goto _test_eof93;
case 93:
#line 567 "evalStringToScalarScanner.cc"
	switch( (*p) ) {
		case 46: goto tr0;
		case 69: goto st2;
		case 101: goto st2;
	}
	if ( 48 <= (*p) && (*p) <= 57 )
		goto tr113;
	goto tr126;
st6:
	if ( ++p == pe )
		goto _test_eof6;
case 6:
	switch( (*p) ) {
		case 99: goto st7;
		case 115: goto st10;
		case 116: goto st13;
	}
	goto st0;
st7:
	if ( ++p == pe )
		goto _test_eof7;
case 7:
	if ( (*p) == 111 )
		goto st8;
	goto st0;
st8:
	if ( ++p == pe )
		goto _test_eof8;
case 8:
	if ( (*p) == 115 )
		goto st9;
	goto st0;
st9:
	if ( ++p == pe )
		goto _test_eof9;
case 9:
	switch( (*p) ) {
		case 32: goto st9;
		case 40: goto tr14;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st9;
	goto st0;
st10:
	if ( ++p == pe )
		goto _test_eof10;
case 10:
	if ( (*p) == 105 )
		goto st11;
	goto st0;
st11:
	if ( ++p == pe )
		goto _test_eof11;
case 11:
	if ( (*p) == 110 )
		goto st12;
	goto st0;
st12:
	if ( ++p == pe )
		goto _test_eof12;
case 12:
	switch( (*p) ) {
		case 32: goto st12;
		case 40: goto tr17;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st12;
	goto st0;
st13:
	if ( ++p == pe )
		goto _test_eof13;
case 13:
	if ( (*p) == 97 )
		goto st14;
	goto st0;
st14:
	if ( ++p == pe )
		goto _test_eof14;
case 14:
	if ( (*p) == 110 )
		goto st15;
	goto st0;
st15:
	if ( ++p == pe )
		goto _test_eof15;
case 15:
	switch( (*p) ) {
		case 32: goto st16;
		case 40: goto tr21;
		case 50: goto st17;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st16;
	goto st0;
st16:
	if ( ++p == pe )
		goto _test_eof16;
case 16:
	switch( (*p) ) {
		case 32: goto st16;
		case 40: goto tr21;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st16;
	goto st0;
st17:
	if ( ++p == pe )
		goto _test_eof17;
case 17:
	switch( (*p) ) {
		case 32: goto st17;
		case 40: goto tr23;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st17;
	goto st0;
st18:
	if ( ++p == pe )
		goto _test_eof18;
case 18:
	switch( (*p) ) {
		case 98: goto st19;
		case 111: goto st22;
	}
	goto st0;
st19:
	if ( ++p == pe )
		goto _test_eof19;
case 19:
	if ( (*p) == 114 )
		goto st20;
	goto st0;
st20:
	if ( ++p == pe )
		goto _test_eof20;
case 20:
	if ( (*p) == 116 )
		goto st21;
	goto st0;
st21:
	if ( ++p == pe )
		goto _test_eof21;
case 21:
	switch( (*p) ) {
		case 32: goto st21;
		case 40: goto tr28;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st21;
	goto st0;
st22:
	if ( ++p == pe )
		goto _test_eof22;
case 22:
	if ( (*p) == 115 )
		goto st23;
	goto st0;
st23:
	if ( ++p == pe )
		goto _test_eof23;
case 23:
	switch( (*p) ) {
		case 32: goto st24;
		case 40: goto tr31;
		case 104: goto st25;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st24;
	goto st0;
st24:
	if ( ++p == pe )
		goto _test_eof24;
case 24:
	switch( (*p) ) {
		case 32: goto st24;
		case 40: goto tr31;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st24;
	goto st0;
st25:
	if ( ++p == pe )
		goto _test_eof25;
case 25:
	switch( (*p) ) {
		case 32: goto st25;
		case 40: goto tr33;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st25;
	goto st0;
st26:
	if ( ++p == pe )
		goto _test_eof26;
case 26:
	if ( (*p) == 101 )
		goto st27;
	goto st0;
st27:
	if ( ++p == pe )
		goto _test_eof27;
case 27:
	if ( (*p) == 103 )
		goto st28;
	goto st0;
st28:
	if ( ++p == pe )
		goto _test_eof28;
case 28:
	if ( (*p) == 84 )
		goto st29;
	goto st0;
st29:
	if ( ++p == pe )
		goto _test_eof29;
case 29:
	if ( (*p) == 111 )
		goto st30;
	goto st0;
st30:
	if ( ++p == pe )
		goto _test_eof30;
case 30:
	if ( (*p) == 82 )
		goto st31;
	goto st0;
st31:
	if ( ++p == pe )
		goto _test_eof31;
case 31:
	if ( (*p) == 97 )
		goto st32;
	goto st0;
st32:
	if ( ++p == pe )
		goto _test_eof32;
case 32:
	if ( (*p) == 100 )
		goto st33;
	goto st0;
st33:
	if ( ++p == pe )
		goto _test_eof33;
case 33:
	switch( (*p) ) {
		case 32: goto st33;
		case 40: goto tr41;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st33;
	goto st0;
st34:
	if ( ++p == pe )
		goto _test_eof34;
case 34:
	if ( (*p) == 120 )
		goto st35;
	goto st0;
st35:
	if ( ++p == pe )
		goto _test_eof35;
case 35:
	if ( (*p) == 112 )
		goto st36;
	goto st0;
st36:
	if ( ++p == pe )
		goto _test_eof36;
case 36:
	switch( (*p) ) {
		case 32: goto st36;
		case 40: goto tr44;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st36;
	goto st0;
st37:
	if ( ++p == pe )
		goto _test_eof37;
case 37:
	if ( (*p) == 121 )
		goto st38;
	goto st0;
st38:
	if ( ++p == pe )
		goto _test_eof38;
case 38:
	if ( (*p) == 112 )
		goto st39;
	goto st0;
st39:
	if ( ++p == pe )
		goto _test_eof39;
case 39:
	if ( (*p) == 111 )
		goto st40;
	goto st0;
st40:
	if ( ++p == pe )
		goto _test_eof40;
case 40:
	if ( (*p) == 116 )
		goto st41;
	goto st0;
st41:
	if ( ++p == pe )
		goto _test_eof41;
case 41:
	switch( (*p) ) {
		case 32: goto st41;
		case 40: goto tr49;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st41;
	goto st0;
st42:
	if ( ++p == pe )
		goto _test_eof42;
case 42:
	if ( (*p) == 111 )
		goto st43;
	goto st0;
st43:
	if ( ++p == pe )
		goto _test_eof43;
case 43:
	if ( (*p) == 103 )
		goto st44;
	goto st0;
st44:
	if ( ++p == pe )
		goto _test_eof44;
case 44:
	switch( (*p) ) {
		case 32: goto st45;
		case 40: goto tr53;
		case 49: goto st46;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st45;
	goto st0;
st45:
	if ( ++p == pe )
		goto _test_eof45;
case 45:
	switch( (*p) ) {
		case 32: goto st45;
		case 40: goto tr53;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st45;
	goto st0;
st46:
	if ( ++p == pe )
		goto _test_eof46;
case 46:
	if ( (*p) == 48 )
		goto st47;
	goto st0;
st47:
	if ( ++p == pe )
		goto _test_eof47;
case 47:
	switch( (*p) ) {
		case 32: goto st47;
		case 40: goto tr56;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st47;
	goto st0;
st48:
	if ( ++p == pe )
		goto _test_eof48;
case 48:
	switch( (*p) ) {
		case 97: goto st49;
		case 105: goto st56;
	}
	goto st0;
st49:
	if ( ++p == pe )
		goto _test_eof49;
case 49:
	switch( (*p) ) {
		case 103: goto st50;
		case 120: goto st55;
	}
	goto st0;
st50:
	if ( ++p == pe )
		goto _test_eof50;
case 50:
	switch( (*p) ) {
		case 32: goto st51;
		case 40: goto tr62;
		case 83: goto st52;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st51;
	goto st0;
st51:
	if ( ++p == pe )
		goto _test_eof51;
case 51:
	switch( (*p) ) {
		case 32: goto st51;
		case 40: goto tr62;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st51;
	goto st0;
st52:
	if ( ++p == pe )
		goto _test_eof52;
case 52:
	if ( (*p) == 113 )
		goto st53;
	goto st0;
st53:
	if ( ++p == pe )
		goto _test_eof53;
case 53:
	if ( (*p) == 114 )
		goto st54;
	goto st0;
st54:
	if ( ++p == pe )
		goto _test_eof54;
case 54:
	switch( (*p) ) {
		case 32: goto st54;
		case 40: goto tr66;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st54;
	goto st0;
st55:
	if ( ++p == pe )
		goto _test_eof55;
case 55:
	switch( (*p) ) {
		case 32: goto st55;
		case 40: goto tr67;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st55;
	goto st0;
st56:
	if ( ++p == pe )
		goto _test_eof56;
case 56:
	if ( (*p) == 110 )
		goto st57;
	goto st0;
st57:
	if ( ++p == pe )
		goto _test_eof57;
case 57:
	switch( (*p) ) {
		case 32: goto st57;
		case 40: goto tr69;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st57;
	goto st0;
st58:
	if ( ++p == pe )
		goto _test_eof58;
case 58:
	switch( (*p) ) {
		case 105: goto st59;
		case 111: goto st60;
	}
	goto st0;
st59:
	if ( ++p == pe )
		goto _test_eof59;
case 59:
	switch( (*p) ) {
		case 32: goto st59;
		case 40: goto tr72;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st59;
	goto st0;
st60:
	if ( ++p == pe )
		goto _test_eof60;
case 60:
	if ( (*p) == 119 )
		goto st61;
	goto st0;
st61:
	if ( ++p == pe )
		goto _test_eof61;
case 61:
	switch( (*p) ) {
		case 32: goto st61;
		case 40: goto tr74;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st61;
	goto st0;
st62:
	if ( ++p == pe )
		goto _test_eof62;
case 62:
	if ( (*p) == 97 )
		goto st63;
	goto st0;
st63:
	if ( ++p == pe )
		goto _test_eof63;
case 63:
	switch( (*p) ) {
		case 100: goto st64;
		case 110: goto st70;
	}
	goto st0;
st64:
	if ( ++p == pe )
		goto _test_eof64;
case 64:
	if ( (*p) == 84 )
		goto st65;
	goto st0;
st65:
	if ( ++p == pe )
		goto _test_eof65;
case 65:
	if ( (*p) == 111 )
		goto st66;
	goto st0;
st66:
	if ( ++p == pe )
		goto _test_eof66;
case 66:
	if ( (*p) == 68 )
		goto st67;
	goto st0;
st67:
	if ( ++p == pe )
		goto _test_eof67;
case 67:
	if ( (*p) == 101 )
		goto st68;
	goto st0;
st68:
	if ( ++p == pe )
		goto _test_eof68;
case 68:
	if ( (*p) == 103 )
		goto st69;
	goto st0;
st69:
	if ( ++p == pe )
		goto _test_eof69;
case 69:
	switch( (*p) ) {
		case 32: goto st69;
		case 40: goto tr83;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st69;
	goto st0;
st70:
	if ( ++p == pe )
		goto _test_eof70;
case 70:
	if ( (*p) == 100 )
		goto st71;
	goto st0;
st71:
	if ( ++p == pe )
		goto _test_eof71;
case 71:
	switch( (*p) ) {
		case 32: goto st71;
		case 40: goto tr85;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st71;
	goto st0;
st72:
	if ( ++p == pe )
		goto _test_eof72;
case 72:
	switch( (*p) ) {
		case 105: goto st73;
		case 113: goto st77;
	}
	goto st0;
st73:
	if ( ++p == pe )
		goto _test_eof73;
case 73:
	if ( (*p) == 110 )
		goto st74;
	goto st0;
st74:
	if ( ++p == pe )
		goto _test_eof74;
case 74:
	switch( (*p) ) {
		case 32: goto st75;
		case 40: goto tr90;
		case 104: goto st76;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st75;
	goto st0;
st75:
	if ( ++p == pe )
		goto _test_eof75;
case 75:
	switch( (*p) ) {
		case 32: goto st75;
		case 40: goto tr90;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st75;
	goto st0;
st76:
	if ( ++p == pe )
		goto _test_eof76;
case 76:
	switch( (*p) ) {
		case 32: goto st76;
		case 40: goto tr92;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st76;
	goto st0;
st77:
	if ( ++p == pe )
		goto _test_eof77;
case 77:
	if ( (*p) == 114 )
		goto st78;
	goto st0;
st78:
	if ( ++p == pe )
		goto _test_eof78;
case 78:
	switch( (*p) ) {
		case 32: goto st79;
		case 40: goto tr95;
		case 116: goto st80;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st79;
	goto st0;
st79:
	if ( ++p == pe )
		goto _test_eof79;
case 79:
	switch( (*p) ) {
		case 32: goto st79;
		case 40: goto tr95;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st79;
	goto st0;
st80:
	if ( ++p == pe )
		goto _test_eof80;
case 80:
	switch( (*p) ) {
		case 32: goto st80;
		case 40: goto tr97;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st80;
	goto st0;
st81:
	if ( ++p == pe )
		goto _test_eof81;
case 81:
	if ( (*p) == 97 )
		goto st82;
	goto st0;
st82:
	if ( ++p == pe )
		goto _test_eof82;
case 82:
	if ( (*p) == 110 )
		goto st83;
	goto st0;
st83:
	if ( ++p == pe )
		goto _test_eof83;
case 83:
	switch( (*p) ) {
		case 32: goto st84;
		case 40: goto tr101;
		case 104: goto st85;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st84;
	goto st0;
st84:
	if ( ++p == pe )
		goto _test_eof84;
case 84:
	switch( (*p) ) {
		case 32: goto st84;
		case 40: goto tr101;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st84;
	goto st0;
st85:
	if ( ++p == pe )
		goto _test_eof85;
case 85:
	switch( (*p) ) {
		case 32: goto st85;
		case 40: goto tr103;
	}
	if ( 9 <= (*p) && (*p) <= 13 )
		goto st85;
	goto st0;
	}
	_test_eof86: cs = 86; goto _test_eof; 
	_test_eof87: cs = 87; goto _test_eof; 
	_test_eof1: cs = 1; goto _test_eof; 
	_test_eof88: cs = 88; goto _test_eof; 
	_test_eof2: cs = 2; goto _test_eof; 
	_test_eof3: cs = 3; goto _test_eof; 
	_test_eof89: cs = 89; goto _test_eof; 
	_test_eof90: cs = 90; goto _test_eof; 
	_test_eof4: cs = 4; goto _test_eof; 
	_test_eof5: cs = 5; goto _test_eof; 
	_test_eof91: cs = 91; goto _test_eof; 
	_test_eof92: cs = 92; goto _test_eof; 
	_test_eof93: cs = 93; goto _test_eof; 
	_test_eof6: cs = 6; goto _test_eof; 
	_test_eof7: cs = 7; goto _test_eof; 
	_test_eof8: cs = 8; goto _test_eof; 
	_test_eof9: cs = 9; goto _test_eof; 
	_test_eof10: cs = 10; goto _test_eof; 
	_test_eof11: cs = 11; goto _test_eof; 
	_test_eof12: cs = 12; goto _test_eof; 
	_test_eof13: cs = 13; goto _test_eof; 
	_test_eof14: cs = 14; goto _test_eof; 
	_test_eof15: cs = 15; goto _test_eof; 
	_test_eof16: cs = 16; goto _test_eof; 
	_test_eof17: cs = 17; goto _test_eof; 
	_test_eof18: cs = 18; goto _test_eof; 
	_test_eof19: cs = 19; goto _test_eof; 
	_test_eof20: cs = 20; goto _test_eof; 
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

	_test_eof: {}
	if ( p == eof )
	{
	switch ( cs ) {
	case 87: goto tr125;
	case 88: goto tr126;
	case 2: goto tr2;
	case 3: goto tr2;
	case 89: goto tr126;
	case 90: goto tr128;
	case 4: goto tr5;
	case 5: goto tr5;
	case 91: goto tr130;
	case 92: goto tr130;
	case 93: goto tr126;
	}
	}

	_out: {}
	}

#line 209 "evalStringToScalarScanner.rl"
  /* ^^^ FSM execution here ^^^ */;

    if (0 == cs)
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
