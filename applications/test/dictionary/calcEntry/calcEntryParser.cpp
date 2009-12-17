

#include <wchar.h>
#include "calcEntryParser.h"


namespace Foam {
namespace functionEntries {
namespace calcEntryInternal {


// ----------------------------------------------------------------------------
// Parser Implementation
// ----------------------------------------------------------------------------

void Parser::SynErr(int n) {
	if (errDist >= minErrDist) errors->SynErr(la->line, la->col, n);
	errDist = 0;
}


void Parser::SemErr(const wchar_t* msg) {
	if (errDist >= minErrDist) errors->Error(t->line, t->col, msg);
	errDist = 0;
}


void Parser::Get() {
	for (;;) {
		t = la;
		la = scanner->Scan();
		if (la->kind <= maxT) {
			++errDist;
			break;
		}

		if (dummyToken != t) {
			dummyToken->kind = t->kind;
			dummyToken->pos = t->pos;
			dummyToken->col = t->col;
			dummyToken->line = t->line;
			dummyToken->next = NULL;
			coco_string_delete(dummyToken->val);
			dummyToken->val = coco_string_create(t->val);
			t = dummyToken;
		}
		la = t;
	}
}


void Parser::Expect(int n) {
	if (la->kind == n) {
		Get();
	}
	else {
		SynErr(n);
	}
}


void Parser::ExpectWeak(int n, int follow) {
	if (la->kind == n) {
		Get();
	}
	else {
		SynErr(n);
		while (!StartOf(follow)) {
			Get();
		}
	}
}


bool Parser::WeakSeparator(int n, int syFol, int repFol) {
	if (la->kind == n) {
		Get();
		return true;
	}
	else if (StartOf(repFol)) {
		return false;
	}
	else {
		SynErr(n);
		while (!(StartOf(syFol) || StartOf(repFol) || StartOf(0))) {
			Get();
		}
		return StartOf(syFol);
	}
}


void Parser::SimpleCalc() {
		val = 0;
		if (debug){Info<<"start val"<< nl;}
		
		if (la->kind == 5) {
			Get();
			Expr(val);
			Expect(6);
		} else if (StartOf(1)) {
			Expr(val);
		} else SynErr(14);
		Expect(0);
}

void Parser::Expr(scalar& val) {
		scalar val2 = 0;
		if (debug) {Info<<"Expr:"<< val<< nl;}
		
		Term(val);
		while (la->kind == 7 || la->kind == 8) {
			if (la->kind == 7) {
				Get();
				Term(val2);
				if (debug) {Info<<"+Term:"<<val2 <<nl;}
				val += val2;
				if (debug) {Info<<"="<< val << nl;}
				
			} else {
				Get();
				Term(val2);
				if (debug) {Info<<"-Term:"<<val2 <<nl;}
				val -= val2;
				if (debug) {Info<<"="<< val << nl;}
				
			}
		}
}

void Parser::Term(scalar& val) {
		scalar val2 = 0;
		if (debug) {Info<<"Term:"<< val<< nl;}
		
		Factor(val);
		while (la->kind == 9 || la->kind == 10) {
			if (la->kind == 9) {
				Get();
				Factor(val2);
				if (debug) {Info<<"*Factor:"<<val2 << nl;}
				val *= val2;
				if (debug) {Info<<"="<< val << nl; }
				
			} else {
				Get();
				Factor(val2);
				if (debug) {Info<<"/Factor:"<<val2 << nl;}
				val /= val2;
				if (debug) {Info<<"="<< val << nl; }
				
			}
		}
}

void Parser::Factor(scalar& val) {
		if (la->kind == 3) {
			Get();
			val = getDictLookup();
			if (debug) {Info<<"lookup:"<<val<<nl;}
			
		} else if (la->kind == 4) {
			Get();
			val = getScalar();
			if (debug) {Info<<"got num:"<<val<<nl;}
			
		} else if (la->kind == 8) {
			Get();
			Expect(11);
			Expr(val);
			Expect(12);
			val = -val;
			if (debug) {Info<<"inv:"<<val<<nl;}
			
		} else if (la->kind == 11) {
			Get();
			Expr(val);
			Expect(12);
			if (debug){Info<<"got Expr:"<<val<<nl;} 
		} else SynErr(15);
}



void Parser::Parse() {
	t = NULL;
	if (dummyToken) {    // safety: someone might call Parse() twice
		delete dummyToken;
	}
	la = dummyToken = new Token();
	la->val = coco_string_create(L"Dummy Token");
	Get();
	SimpleCalc();

	Expect(0);
}


Parser::Parser(Scanner* scan, Errors* err)
:
	dummyToken(NULL),
	deleteErrorsDestruct_(!err),
	minErrDist(2),
	errDist(minErrDist),
	scanner(scan),
	errors(err),
	t(NULL),
	la(NULL)
{
	maxT = 13;

	if (!errors) {   // add in default error handling
		errors = new Errors();
	}
}


bool Parser::StartOf(int s) {
	const bool T = true;
	const bool x = false;

	static bool set[2][15] = {
		{T,x,x,x, x,x,x,x, x,x,x,x, x,x,x},
		{x,x,x,T, T,x,x,x, T,x,x,T, x,x,x}
	};



	return set[s][la->kind];
}


Parser::~Parser() {
	if (deleteErrorsDestruct_) {    // delete default error handling
		delete errors;
	}
	delete dummyToken;
}


// ----------------------------------------------------------------------------
// Errors Implementation
// ----------------------------------------------------------------------------

Errors::Errors()
:
	count(0)
{}


Errors::~Errors()
{}


void Errors::clear() {
	count = 0;
}


wchar_t* Errors::strerror(int n)
{
	wchar_t* s;
	switch (n) {
			case 0: s = coco_string_create(L"EOF expected"); break;
			case 1: s = coco_string_create(L"ident expected"); break;
			case 2: s = coco_string_create(L"string expected"); break;
			case 3: s = coco_string_create(L"variable expected"); break;
			case 4: s = coco_string_create(L"number expected"); break;
			case 5: s = coco_string_create(L"\"{\" expected"); break;
			case 6: s = coco_string_create(L"\"}\" expected"); break;
			case 7: s = coco_string_create(L"\"+\" expected"); break;
			case 8: s = coco_string_create(L"\"-\" expected"); break;
			case 9: s = coco_string_create(L"\"*\" expected"); break;
			case 10: s = coco_string_create(L"\"/\" expected"); break;
			case 11: s = coco_string_create(L"\"(\" expected"); break;
			case 12: s = coco_string_create(L"\")\" expected"); break;
			case 13: s = coco_string_create(L"??? expected"); break;
			case 14: s = coco_string_create(L"invalid SimpleCalc"); break;
			case 15: s = coco_string_create(L"invalid Factor"); break;

		default:
		{
			wchar_t format[20];
			coco_swprintf(format, 20, L"error %d", n);
			s = coco_string_create(format);
		}
		break;
	}
	return s;
}


void Errors::Warning(const wchar_t* msg) {
	wprintf(L"%ls\n", msg);
}


void Errors::Warning(int line, int col, const wchar_t* msg) {
	wprintf(L"-- line %d col %d: %ls\n", line, col, msg);
}


void Errors::Error(int line, int col, const wchar_t* msg) {
	wprintf(L"-- line %d col %d: %ls\n", line, col, msg);
	count++;
}


void Errors::SynErr(int line, int col, int n) {
	wchar_t* msg = this->strerror(n);
	this->Error(line, col, msg);
	coco_string_delete(msg);
}


void Errors::Exception(const wchar_t* msg) {
	wprintf(L"%ls", msg);
	::exit(1);
}

} // namespace
} // namespace
} // namespace


