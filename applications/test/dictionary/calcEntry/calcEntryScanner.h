

#ifndef COCO_calcEntrySCANNER_H__
#define COCO_calcEntrySCANNER_H__

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <wchar.h>

// io.h and fcntl are used to ensure binary read from streams on windows
#if _MSC_VER >= 1300
#include <io.h>
#include <fcntl.h>
#endif

#if _MSC_VER >= 1400
#define coco_swprintf swprintf_s
#elif _MSC_VER >= 1300
#define coco_swprintf _snwprintf
#else
// assume every other compiler knows swprintf
#define coco_swprintf swprintf
#endif

#define COCO_WCHAR_MAX    65535
#define MIN_BUFFER_LENGTH 1024
#define MAX_BUFFER_LENGTH (64*MIN_BUFFER_LENGTH)
#define HEAP_BLOCK_SIZE   (64*1024)


namespace Foam {
namespace functionEntries {
namespace calcEntryInternal {



// * * * * * * * * * *  Wide Character String Routines * * * * * * * * * * * //

//
// string handling, wide character
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! Create by copying str
wchar_t* coco_string_create(const wchar_t* str);

//! Create a substring of str starting at index and length characters long
wchar_t* coco_string_create(const wchar_t* str, int index, int length);

//! Create an uppercase string from str
wchar_t* coco_string_create_upper(const wchar_t* str);

//! Create an uppercase substring from str starting at index and length characters long
wchar_t* coco_string_create_upper(const wchar_t* str, int index, int length);

//! Create a lowercase string from str
wchar_t* coco_string_create_lower(const wchar_t* str);

//! Create a lowercase substring from str starting at index and length characters long
wchar_t* coco_string_create_lower(const wchar_t* str, int index, int length);

//! Create a string by concatenating str1 and str2
wchar_t* coco_string_create_append(const wchar_t* str1, const wchar_t* str2);

//! Create a string by concatenating a character to the end of str
wchar_t* coco_string_create_append(const wchar_t* str, const wchar_t ch);

//! Free storage and nullify the argument
void  coco_string_delete(wchar_t* &str);

//! The length of the str, or 0 if the str is NULL
int   coco_string_length(const wchar_t* str);

//! Return true if the str ends with the endstr
bool  coco_string_endswith(const wchar_t* str, const wchar_t* endstr);

//! Return the index of the first occurrence of ch.
//  Return -1 if nothing is found.
int   coco_string_indexof(const wchar_t* str, const wchar_t ch);

//! Return the index of the last occurrence of ch.
//  Return -1 if nothing is found.
int   coco_string_lastindexof(const wchar_t* str, const wchar_t ch);

//! Append str to dest
void  coco_string_merge(wchar_t* &dest, const wchar_t* str);

//! Compare strings, return true if they are equal
bool  coco_string_equal(const wchar_t* str1, const wchar_t* str2);

//! Compare strings, return 0 if they are equal
int   coco_string_compareto(const wchar_t* str1, const wchar_t* str2);

//! Simple string hashing function
int   coco_string_hash(const wchar_t* str);

//
// String conversions
// ~~~~~~~~~~~~~~~~~~

//! Convert wide string to double
double coco_string_toDouble(const wchar_t* str);

//! Convert wide string to float
float coco_string_toFloat(const wchar_t* str);

//
// String handling, byte character
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! Create by copying byte str
wchar_t* coco_string_create(const char* str);

//! Create a substring of byte str starting at index and length characters long
wchar_t* coco_string_create(const char* str, int index, int length);

//! Create a byte string by copying str
char* coco_string_create_char(const wchar_t* str);

//! Create a byte substring of str starting at index and length characters long
char* coco_string_create_char(const wchar_t* str, int index, int length);

//! Free storage and nullify the argument
void  coco_string_delete(char* &str);


//
// String conversions, byte character
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//! Convert byte string to double
double coco_string_toDouble(const char* str);

//! Convert byte string to float
float coco_string_toFloat(const char* str);

// * * * * * * * * * End of Wide Character String Routines * * * * * * * * * //



//! Scanner Token
class Token
{
public:
	int kind;       //!< token kind
	int pos;        //!< token position in the source text (starting at 0)
	int col;        //!< token column (starting at 1)
	int line;       //!< token line (starting at 1)
	wchar_t* val;   //!< token value
	Token *next;    //!< Peek tokens are kept in linked list

	Token();        //!< Construct null
	~Token();       //!< Destructor - cleanup allocated val
};


//! Scanner Buffer
//
//! This Buffer supports the following cases:
//! -# seekable stream (file)
//!    -# whole stream in buffer
//!    -# part of stream in buffer
//! -# non seekable stream (network, console)
class Buffer {
private:
	unsigned char *buf; //!< input buffer
	int bufCapacity;    //!< capacity of buf
	int bufStart;       //!< position of first byte in buffer relative to input stream
	int bufLen;         //!< length of buffer
	int fileLen;        //!< length of input stream (may change if the stream is no file)
	int bufPos;         //!< current position in buffer
	FILE* stream;       //!< input stream (seekable)
	bool isUserStream;  //!< was the stream opened by the user?

	int ReadNextStreamChunk();
	bool CanSeek();     //!< true if stream can be seeked otherwise false

public:
	static const int EoF = COCO_WCHAR_MAX + 1;

	Buffer(FILE*, bool isUserStream);
	Buffer(const unsigned char* buf, int len);
	Buffer(const char* buf, int len);
	Buffer(Buffer*);
	virtual ~Buffer();

	virtual void Close();
	virtual int Read();
	virtual int Peek();
	virtual wchar_t* GetString(int beg, int end);
	virtual int GetPos();
	virtual void SetPos(int value);
};


//! A Scanner buffer that handles UTF-8 characters
class UTF8Buffer : public Buffer {
public:
	UTF8Buffer(Buffer* b) : Buffer(b) {}
	virtual int Read();
};


//------------------------------------------------------------------------------
// StartStates
//------------------------------------------------------------------------------
//! maps characters to start states of tokens
class StartStates {
private:
	class Elem {
	public:
		int key, val;
		Elem *next;
		Elem(int key, int val) {
			this->key = key;
			this->val = val;
			next = NULL;
		}
	};

	Elem **tab;

public:
	StartStates() {
		tab = new Elem*[128];
		memset(tab, 0, 128 * sizeof(Elem*));
	}
	virtual ~StartStates() {
		for (int i = 0; i < 128; ++i) {
			Elem *e = tab[i];
			while (e != NULL) {
				Elem *next = e->next;
				delete e;
				e = next;
			}
		}
		delete [] tab;
	}

	void set(int key, int val) {
		Elem *e = new Elem(key, val);
		int k = ((unsigned int) key) % 128;
		e->next = tab[k];
		tab[k] = e;
	}

	int state(int key) {
		Elem *e = tab[((unsigned int) key) % 128];
		while (e != NULL && e->key != key) e = e->next;
		return e == NULL ? 0 : e->val;
	}
};


//------------------------------------------------------------------------------
// KeywordMap
//------------------------------------------------------------------------------
//! maps strings to integers (identifiers to keyword kinds)
class KeywordMap {
private:
	class Elem {
	public:
		wchar_t *key;
		int val;
		Elem *next;
		Elem(const wchar_t *key, int val) {
			this->key = coco_string_create(key);
			this->val = val;
			next = NULL;
		}
		virtual ~Elem() {
			coco_string_delete(key);
		}
	};

	Elem **tab;

public:
	KeywordMap() {
		tab = new Elem*[128];
		memset(tab, 0, 128 * sizeof(Elem*));
	}
	virtual ~KeywordMap() {
		for (int i = 0; i < 128; ++i) {
			Elem *e = tab[i];
			while (e != NULL) {
				Elem *next = e->next;
				delete e;
				e = next;
			}
		}
		delete [] tab;
	}

	void set(const wchar_t *key, int val) {
		Elem *e = new Elem(key, val);
		int k = coco_string_hash(key) % 128;
		e->next = tab[k]; tab[k] = e;
	}

	int get(const wchar_t *key, int defaultVal) {
		Elem *e = tab[coco_string_hash(key) % 128];
		while (e != NULL && !coco_string_equal(e->key, key)) e = e->next;
		return e == NULL ? defaultVal : e->val;
	}
};


//! A Coco/R Scanner
class Scanner {
private:
	static const unsigned char EOL = '\n';   // end-of-line character
	static const int eofSym = 0;             // end-of-file token id

	void *firstHeap;
	void *heap;
	void *heapTop;
	void **heapEnd;

	int noSym;        //!< noSym gets highest number, set in Parser
	int maxT;
	int charSetSize;  //!< unused?
	StartStates start;
	KeywordMap keywords;

	Token *t;         //!< current token
	wchar_t *tval;    //!< text of current token
	int tvalLength;   //!< length of text of current token
	int tlen;         //!< length of current token

	Token *tokens;    //!< list of tokens already peeked (first token is a dummy)
	Token *pt;        //!< current peek token

	int ch;           //!< current input character

	int pos;          //!< byte position of current character
	int line;         //!< line number of current character
	int col;          //!< column number of current character
	int oldEols;      //!< EOLs that appeared in a comment;

	void CreateHeapBlock();
	Token* CreateToken();
	void AppendVal(Token*);

	void Init();
	void NextCh();
	void AddCh();
	bool Comment0();
	bool Comment1();

	Token* NextToken();

public:
	//! scanner buffer
	Buffer *buffer;

	//! Attach scanner to an existing character buffer
	Scanner(const unsigned char* buf, int len);
	//! Attach scanner to an existing character buffer
	Scanner(const char* buf, int len);
	//! Open a file for reading and attach scanner
	Scanner(const wchar_t* fileName);
	//! Using an existing open file handle for the scanner
	Scanner(FILE* s);
	~Scanner();
	Token* Scan();
	Token* Peek();
	void ResetPeek();

}; // end Scanner

} // namespace
} // namespace
} // namespace


#endif // COCO_calcEntrySCANNER_H__

