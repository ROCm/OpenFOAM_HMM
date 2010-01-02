

#ifndef COCO_calcEntryPARSER_H__
#define COCO_calcEntryPARSER_H__

#include "dictionary.H"
#include "wchar.H"
#include "calcEntryInternal.H"


#include "calcEntryScanner.h"

namespace Foam {
namespace functionEntries {
namespace calcEntryInternal {


//! Parser error handing
class Errors {
public:
	int count;      //!< The number of errors detected

	//! Allocate and return a string describing the given error code.
	/** It is the responsibility of the caller to free this string,
	 *  eg, with coco_string_delete()
	 */
	static wchar_t* strerror(int n);

	Errors();               //!< Construct null - start with no errors
	virtual ~Errors();      //!< Destructor
	virtual void clear();   //!< Clear the error count

	//! Handle a general warning 'msg'
	virtual void Warning(const wchar_t* msg);
	//! Handle a general warning 'msg'
	virtual void Warning(int line, int col, const wchar_t* msg);
	//! Handle general error 'msg' (eg, a semantic error)
	virtual void Error(int line, int col, const wchar_t* msg);
	//! Handle syntax error 'n', uses strerror for the message, calls Error()
	virtual void SynErr(int line, int col, int n);
	//! Handle a general exception 'msg'
	virtual void Exception(const wchar_t* msg);

}; // Errors


//! A Coco/R Parser
class Parser {
private:
	enum {
		_EOF=0,
		_ident=1,
		_string=2,
		_variable=3,
		_number=4,
	};
	static const int maxT = 14;

	static const int minErrDist = 2; //!< min. distance before reporting errors

	Token *dummyToken;
	bool deleteErrorsDestruct_; //!< delete the 'errors' member in destructor
	int  errDist;

	void SynErr(int n);         //!< Handle syntax error 'n'
	void Get();
	void Expect(int n);
	bool StartOf(int s);
	void ExpectWeak(int n, int follow);
	bool WeakSeparator(int n, int syFol, int repFol);

public:
	Scanner *scanner;
	Errors  *errors;

	Token *t;                   //!< last recognized token
	Token *la;                  //!< lookahead token

private:
    //- The parent dictionary
    dictionary* dict_;

    //- The calculation result
    scalar val;

    //- lookup dictionary entry
    scalar getDictLookup() const
    {
        scalar dictValue = 0;

        if (!dict_)
        {
            FatalErrorIn
            (
                "calcEntry::getDictEntry() const"
            )   << "No dictionary attached!"
                << exit(FatalError);

            return 0;
        }

        char* str = coco_string_create_char
        (
            t->val,
            1,
            (coco_string_length(t->val) - 1)
        );
        word keyword(str);
        coco_string_delete(str);

        entry* entryPtr = dict_->lookupEntryPtr(keyword, true, false);
        if (entryPtr && !entryPtr->isDict())
        {
            if (entryPtr->stream().size() != 1)
            {
                FatalErrorIn
                (
                    "calcEntry::getDictEntry() const"
                )   << "keyword " << keyword << " has "
                    << entryPtr->stream().size() << " values in dictionary "
                    << exit(FatalError);
            }
            entryPtr->stream() >> dictValue;
        }
        else
        {
            FatalErrorIn
            (
                "calcEntry::getDictEntry() const"
            )   << "keyword " << keyword << " is undefined in dictionary "
                << exit(FatalError);
        }

        return dictValue;
    }


public:

    //- attach a dictionary
    void dict(const dictionary& dict)
    {
        dict_ = const_cast<dictionary*>(&dict);
    }

    //- Return the calculated result
    scalar Result() const
    {
        return val;
    }




	//! Construct for the specified scanner
	/**
	 *  Use the default error handling, or optionally provide an error
	 *  handler, which will not be deleted upon destruction.
	 */
	Parser(Scanner* scan, Errors* err = 0);
	~Parser();
	void SemErr(const wchar_t* msg);    //!< Handle semantic error

	void calcEntry();
	void Expr(scalar& val);
	void Term(scalar& val);
	void Factor(scalar& val);
	void Func(scalar& val);

	void Parse();                       //!< Execute the parse operation

}; // end Parser

} // namespace
} // namespace
} // namespace


#endif // COCO_calcEntryPARSER_H__

