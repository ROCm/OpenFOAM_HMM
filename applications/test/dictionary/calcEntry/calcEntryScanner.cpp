

#include <memory.h>
#include <string.h>
#include "calcEntryScanner.h"

namespace Foam {
namespace functionEntries {
namespace calcEntryInternal {


// * * * * * * * * * *  Wide Character String Routines * * * * * * * * * * * //

// string handling, wide character

wchar_t* coco_string_create(const wchar_t* str) {
	int len = coco_string_length(str);
	wchar_t* dest = new wchar_t[len + 1];
	if (len) {
		wcsncpy(dest, str, len);
	}
	dest[len] = 0;
	return dest;
}

wchar_t* coco_string_create(const wchar_t* str, int index, int length) {
	int len = coco_string_length(str);
	if (len) {
		len = length;
	}
	wchar_t* dest = new wchar_t[len + 1];
	if (len) {
		wcsncpy(dest, &(str[index]), len);
	}
	dest[len] = 0;
	return dest;
}

wchar_t* coco_string_create_upper(const wchar_t* str) {
	if (!str) { return NULL; }
	return coco_string_create_upper(str, 0, wcslen(str));
}


wchar_t* coco_string_create_upper(const wchar_t* str, int index, int len) {
	if (!str) { return NULL; }
	wchar_t* dest = new wchar_t[len + 1];

	for (int i = 0; i < len; i++) {
		const wchar_t ch = str[index + i];
		if ((L'a' <= ch) && (ch <= L'z')) {
			dest[i] = ch + (L'A' - L'a');
		}
		else {
			dest[i] = ch;
		}
	}
	dest[len] = L'\0';
	return dest;
}


wchar_t* coco_string_create_lower(const wchar_t* str) {
	if (!str) { return NULL; }
	return coco_string_create_lower(str, 0, wcslen(str));
}


wchar_t* coco_string_create_lower(const wchar_t* str, int index, int len) {
	if (!str) { return NULL; }
	wchar_t* dest = new wchar_t[len + 1];

	for (int i = 0; i < len; i++) {
		const wchar_t ch = str[index + i];
		if ((L'A' <= ch) && (ch <= L'Z')) {
			dest[i] = ch - (L'A' - L'a');
		}
		else {
			dest[i] = ch;
		}
	}
	dest[len] = L'\0';
	return dest;
}


wchar_t* coco_string_create_append(const wchar_t* str1, const wchar_t* str2) {
	int str1Len = coco_string_length(str1);
	int str2Len = coco_string_length(str2);

	wchar_t* dest = new wchar_t[str1Len + str2Len + 1];

	if (str1Len) { wcscpy(dest, str1); }
	if (str2Len) { wcscpy(dest + str1Len, str2); }

	dest[str1Len + str2Len] = 0;
	return dest;
}

wchar_t* coco_string_create_append(const wchar_t* str1, const wchar_t ch) {
	int len = coco_string_length(str1);
	wchar_t* dest = new wchar_t[len + 2];
	wcsncpy(dest, str1, len);   // or use if (len) { wcscpy(dest, str1); }
	dest[len] = ch;
	dest[len + 1] = 0;
	return dest;
}

void coco_string_delete(wchar_t* &str) {
	delete [] str;
	str = NULL;
}

int coco_string_length(const wchar_t* str) {
	return str ? wcslen(str) : 0;
}

bool coco_string_endswith(const wchar_t* str, const wchar_t* endstr) {
	int strLen = wcslen(str);
	int endLen = wcslen(endstr);
	return (endLen <= strLen) && (wcscmp(str + strLen - endLen, endstr) == 0);
}

int coco_string_indexof(const wchar_t* str, const wchar_t ch) {
	const wchar_t* fnd = wcschr(str, ch);
	return fnd ? (fnd - str) : -1;
}

int coco_string_lastindexof(const wchar_t* str, const wchar_t ch) {
	const wchar_t* fnd = wcsrchr(str, ch);
	return fnd ? (fnd - str) : -1;
}

void coco_string_merge(wchar_t* &dest, const wchar_t* str) {
	if (!str) { return; }
	wchar_t* newstr = coco_string_create_append(dest, str);
	delete [] dest;
	dest = newstr;
}

bool coco_string_equal(const wchar_t* str1, const wchar_t* str2) {
	return wcscmp(str1, str2) == 0;
}

int coco_string_compareto(const wchar_t* str1, const wchar_t* str2) {
	return wcscmp(str1, str2);
}

int coco_string_hash(const wchar_t* str) {
	int h = 0;
	if (!str) { return 0; }
	while (*str != 0) {
		h = (h * 7) ^ *str;
		++str;
	}
	if (h < 0) { h = -h; }
	return h;
}


double coco_string_toDouble(const wchar_t* str)
{
	return str ? wcstod(str, NULL) : 0;
}

float coco_string_toFloat(const wchar_t* str)
{
	return str ? wcstof(str, NULL) : 0;
}



//
// string handling, byte character
//

wchar_t* coco_string_create(const char* str) {
	int len = str ? strlen(str) : 0;
	wchar_t* dest = new wchar_t[len + 1];
	for (int i = 0; i < len; ++i) {
		dest[i] = (wchar_t) str[i];
	}
	dest[len] = 0;
	return dest;
}

wchar_t* coco_string_create(const char* str, int index, int length) {
	int len = str ? length : 0;
	wchar_t* dest = new wchar_t[len + 1];
	for (int i = 0; i < len; ++i) {
		dest[i] = (wchar_t) str[index + i];
	}
	dest[len] = 0;
	return dest;
}


char* coco_string_create_char(const wchar_t* str) {
	int len = coco_string_length(str);
	char *dest = new char[len + 1];
	for (int i = 0; i < len; ++i)
	{
		dest[i] = (char) str[i];
	}
	dest[len] = 0;
	return dest;
}

char* coco_string_create_char(const wchar_t* str, int index, int length) {
	int len = coco_string_length(str);
	if (len) {
		len = length;
	}
	char *dest = new char[len + 1];
	for (int i = 0; i < len; ++i) {
		dest[i] = (char) str[index + i];
	}
	dest[len] = 0;
	return dest;
}


void coco_string_delete(char* &str) {
	delete [] str;
	str = NULL;
}


double coco_string_toDouble(const char* str)
{
	return str ? strtod(str, NULL) : 0;
}

float coco_string_toFloat(const char* str)
{
	return str ? strtof(str, NULL) : 0;
}


// * * * * * * * * * End of Wide Character String Routines * * * * * * * * * //


Token::Token()
:
    kind(0),
    pos(0),
    col(0),
    line(0),
    val(NULL),
    next(NULL)
{}


Token::~Token() {
	coco_string_delete(val);
}


Buffer::Buffer(FILE* s, bool isUserStream) {
// ensure binary read on windows
#if _MSC_VER >= 1300
	_setmode(_fileno(s), _O_BINARY);
#endif
	stream = s; this->isUserStream = isUserStream;
	if (CanSeek()) {
		fseek(s, 0, SEEK_END);
		fileLen = ftell(s);
		fseek(s, 0, SEEK_SET);
		bufLen = (fileLen < MAX_BUFFER_LENGTH) ? fileLen : MAX_BUFFER_LENGTH;
		bufStart = INT_MAX; // nothing in the buffer so far
	}
	else {
		fileLen = bufLen = bufStart = 0;
	}
	bufCapacity = (bufLen > 0) ? bufLen : MIN_BUFFER_LENGTH;
	buf = new unsigned char[bufCapacity];
	if (fileLen > 0) SetPos(0);          // setup buffer to position 0 (start)
	else bufPos = 0; // index 0 is already after the file, thus Pos = 0 is invalid
	if (bufLen == fileLen && CanSeek()) Close();
}


Buffer::Buffer(Buffer* b) {
	buf = b->buf;
	bufCapacity = b->bufCapacity;
	b->buf = NULL;
	bufStart = b->bufStart;
	bufLen = b->bufLen;
	fileLen = b->fileLen;
	bufPos = b->bufPos;
	stream = b->stream;
	b->stream = NULL;
	isUserStream = b->isUserStream;
}


Buffer::Buffer(const unsigned char* buf, int len) {
	this->buf = new unsigned char[len];
	memcpy(this->buf, buf, len*sizeof(unsigned char));
	bufStart = 0;
	bufCapacity = bufLen = len;
	fileLen = len;
	bufPos = 0;
	stream = NULL;
}


Buffer::Buffer(const char* buf, int len) {
	this->buf = new unsigned char[len];
	memcpy(this->buf, buf, len*sizeof(unsigned char));
	bufStart = 0;
	bufCapacity = bufLen = len;
	fileLen = len;
	bufPos = 0;
	stream = NULL;
}


Buffer::~Buffer() {
	Close();
	if (buf != NULL) {
		delete [] buf;
		buf = NULL;
	}
}


void Buffer::Close() {
	if (!isUserStream && stream != NULL) {
		fclose(stream);
		stream = NULL;
	}
}


int Buffer::Read() {
	if (bufPos < bufLen) {
		return buf[bufPos++];
	} else if (GetPos() < fileLen) {
		SetPos(GetPos()); // shift buffer start to Pos
		return buf[bufPos++];
	} else if ((stream != NULL) && !CanSeek() && (ReadNextStreamChunk() > 0)) {
		return buf[bufPos++];
	} else {
		return EoF;
	}
}


int Buffer::Peek() {
	int curPos = GetPos();
	int ch = Read();
	SetPos(curPos);
	return ch;
}


wchar_t* Buffer::GetString(int beg, int end) {
	int len = 0;
	wchar_t *buf = new wchar_t[end - beg];
	int oldPos = GetPos();
	SetPos(beg);
	while (GetPos() < end) buf[len++] = (wchar_t) Read();
	SetPos(oldPos);
	wchar_t *res = coco_string_create(buf, 0, len);
	coco_string_delete(buf);
	return res;
}


int Buffer::GetPos() {
	return bufPos + bufStart;
}


void Buffer::SetPos(int value) {
	if ((value >= fileLen) && (stream != NULL) && !CanSeek()) {
		// Wanted position is after buffer and the stream
		// is not seek-able e.g. network or console,
		// thus we have to read the stream manually till
		// the wanted position is in sight.
		while ((value >= fileLen) && (ReadNextStreamChunk() > 0))
		{}
	}

	if ((value < 0) || (value > fileLen)) {
		wprintf(L"--- buffer out of bounds access, position: %d\n", value);
		::exit(1);
	}

	if ((value >= bufStart) && (value < (bufStart + bufLen))) { // already in buffer
		bufPos = value - bufStart;
	} else if (stream != NULL) { // must be swapped in
		fseek(stream, value, SEEK_SET);
		bufLen = fread(buf, sizeof(unsigned char), bufCapacity, stream);
		bufStart = value; bufPos = 0;
	} else {
		bufPos = fileLen - bufStart; // make Pos return fileLen
	}
}


// Read the next chunk of bytes from the stream, increases the buffer
// if needed and updates the fields fileLen and bufLen.
// Returns the number of bytes read.
int Buffer::ReadNextStreamChunk() {
	int freeLen = bufCapacity - bufLen;
	if (freeLen == 0) {
		// in the case of a growing input stream
		// we can neither seek in the stream, nor can we
		// foresee the maximum length, thus we must adapt
		// the buffer size on demand.
		bufCapacity = bufLen * 2;
		unsigned char *newBuf = new unsigned char[bufCapacity];
		memcpy(newBuf, buf, bufLen*sizeof(unsigned char));
		delete [] buf;
		buf = newBuf;
		freeLen = bufLen;
	}
	int read = fread(buf + bufLen, sizeof(unsigned char), freeLen, stream);
	if (read > 0) {
		fileLen = bufLen = (bufLen + read);
		return read;
	}
	// end of stream reached
	return 0;
}


bool Buffer::CanSeek() {
	return (stream != NULL) && (ftell(stream) != -1);
}


int UTF8Buffer::Read() {
	int ch;
	do {
		ch = Buffer::Read();
		// until we find a utf8 start (0xxxxxxx or 11xxxxxx)
	} while ((ch >= 128) && ((ch & 0xC0) != 0xC0) && (ch != EoF));
	if (ch < 128 || ch == EoF) {
		// nothing to do, first 127 chars are the same in ascii and utf8
		// 0xxxxxxx or end of file character
	} else if ((ch & 0xF0) == 0xF0) {
		// 11110xxx 10xxxxxx 10xxxxxx 10xxxxxx
		int c1 = ch & 0x07; ch = Buffer::Read();
		int c2 = ch & 0x3F; ch = Buffer::Read();
		int c3 = ch & 0x3F; ch = Buffer::Read();
		int c4 = ch & 0x3F;
		ch = (((((c1 << 6) | c2) << 6) | c3) << 6) | c4;
	} else if ((ch & 0xE0) == 0xE0) {
		// 1110xxxx 10xxxxxx 10xxxxxx
		int c1 = ch & 0x0F; ch = Buffer::Read();
		int c2 = ch & 0x3F; ch = Buffer::Read();
		int c3 = ch & 0x3F;
		ch = (((c1 << 6) | c2) << 6) | c3;
	} else if ((ch & 0xC0) == 0xC0) {
		// 110xxxxx 10xxxxxx
		int c1 = ch & 0x1F; ch = Buffer::Read();
		int c2 = ch & 0x3F;
		ch = (c1 << 6) | c2;
	}
	return ch;
}


Scanner::Scanner(const unsigned char* buf, int len) {
	buffer = new Buffer(buf, len);
	Init();
}


Scanner::Scanner(const char* buf, int len) {
	buffer = new Buffer(buf, len);
	Init();
}


Scanner::Scanner(const wchar_t* fileName) {
	FILE* stream;
	char *chFileName = coco_string_create_char(fileName);
	if ((stream = fopen(chFileName, "rb")) == NULL) {
		wprintf(L"--- Cannot open file %ls\n", fileName);
		::exit(1);
	}
	coco_string_delete(chFileName);
	buffer = new Buffer(stream, false);
	Init();
}


Scanner::Scanner(FILE* s) {
	buffer = new Buffer(s, true);
	Init();
}


Scanner::~Scanner() {
	char* cur = (char*) firstHeap;

	while (cur != NULL) {
		cur = *(char**) (cur + HEAP_BLOCK_SIZE);
		free(firstHeap);
		firstHeap = cur;
	}
	delete [] tval;
	delete buffer;
}


void Scanner::Init() {
	maxT = 13;
	noSym = 13;
	int i;
	for (i = 65; i <= 90; ++i) start.set(i, 1);
	for (i = 97; i <= 122; ++i) start.set(i, 1);
	for (i = 36; i <= 36; ++i) start.set(i, 5);
	start.set(45, 20);
	for (i = 48; i <= 57; ++i) start.set(i, 9);
	start.set(34, 2);
	start.set(46, 7);
	start.set(123, 14);
	start.set(125, 15);
	start.set(43, 21);
	start.set(42, 16);
	start.set(47, 17);
	start.set(40, 18);
	start.set(41, 19);
		start.set(Buffer::EoF, -1);


	tvalLength = 128;
	tval = new wchar_t[tvalLength]; // text of current token

	// HEAP_BLOCK_SIZE byte heap + pointer to next heap block
	heap = malloc(HEAP_BLOCK_SIZE + sizeof(void*));
	firstHeap = heap;
	heapEnd = (void**) (((char*) heap) + HEAP_BLOCK_SIZE);
	*heapEnd = 0;
	heapTop = heap;
	if (sizeof(Token) > HEAP_BLOCK_SIZE) {
		wprintf(L"--- Too small HEAP_BLOCK_SIZE\n");
		::exit(1);
	}

	pos = -1; line = 1; col = 0;
	oldEols = 0;
	NextCh();
	if (ch == 0xEF) { // check optional byte order mark for UTF-8
		NextCh(); int ch1 = ch;
		NextCh(); int ch2 = ch;
		if (ch1 != 0xBB || ch2 != 0xBF) {
			wprintf(L"Illegal byte order mark at start of file");
			::exit(1);
		}
		Buffer *oldBuf = buffer;
		buffer = new UTF8Buffer(buffer); col = 0;
		delete oldBuf; oldBuf = NULL;
		NextCh();
	}


	pt = tokens = CreateToken(); // first token is a dummy
}


void Scanner::NextCh() {
	if (oldEols > 0) {
		ch = EOL;
		oldEols--;
	}
	else {
		pos = buffer->GetPos();
		ch = buffer->Read(); col++;
		// replace isolated '\r' by '\n' in order to make
		// eol handling uniform across Windows, Unix and Mac
		if (ch == L'\r' && buffer->Peek() != L'\n') ch = EOL;
		if (ch == EOL) { line++; col = 0; }
	}

}


void Scanner::AddCh() {
	if (tlen >= tvalLength) {
		tvalLength *= 2;
		wchar_t *newBuf = new wchar_t[tvalLength];
		memcpy(newBuf, tval, tlen*sizeof(wchar_t));
		delete [] tval;
		tval = newBuf;
	}
	if (ch != Buffer::EoF) {
		tval[tlen++] = ch;
		NextCh();
	}
}



bool Scanner::Comment0() {
	int level = 1, pos0 = pos, line0 = line, col0 = col;
	NextCh();
	if (ch == L'/') {
		NextCh();
		for(;;) {
			if (ch == 10) {
				level--;
				if (level == 0) { oldEols = line - line0; NextCh(); return true; }
				NextCh();
			} else if (ch == buffer->EoF) return false;
			else NextCh();
		}
	} else {
		buffer->SetPos(pos0); NextCh(); line = line0; col = col0;
	}
	return false;
}

bool Scanner::Comment1() {
	int level = 1, pos0 = pos, line0 = line, col0 = col;
	NextCh();
	if (ch == L'*') {
		NextCh();
		for(;;) {
			if (ch == L'*') {
				NextCh();
				if (ch == L'/') {
					level--;
					if (level == 0) { oldEols = line - line0; NextCh(); return true; }
					NextCh();
				}
			} else if (ch == L'/') {
				NextCh();
				if (ch == L'*') {
					level++; NextCh();
				}
			} else if (ch == buffer->EoF) return false;
			else NextCh();
		}
	} else {
		buffer->SetPos(pos0); NextCh(); line = line0; col = col0;
	}
	return false;
}


void Scanner::CreateHeapBlock() {
	void* newHeap;
	char* cur = (char*) firstHeap;

	while (((char*) tokens < cur) || ((char*) tokens > (cur + HEAP_BLOCK_SIZE))) {
		cur = *((char**) (cur + HEAP_BLOCK_SIZE));
		free(firstHeap);
		firstHeap = cur;
	}

	// HEAP_BLOCK_SIZE byte heap + pointer to next heap block
	newHeap = malloc(HEAP_BLOCK_SIZE + sizeof(void*));
	*heapEnd = newHeap;
	heapEnd = (void**) (((char*) newHeap) + HEAP_BLOCK_SIZE);
	*heapEnd = 0;
	heap = newHeap;
	heapTop = heap;
}


Token* Scanner::CreateToken() {
	Token *t;
	if (((char*) heapTop + (int) sizeof(Token)) >= (char*) heapEnd) {
		CreateHeapBlock();
	}
	t = (Token*) heapTop;
	heapTop = (void*) ((char*) heapTop + sizeof(Token));
	t->val = NULL;
	t->next = NULL;
	return t;
}


void Scanner::AppendVal(Token *t) {
	int reqMem = (tlen + 1) * sizeof(wchar_t);
	if (((char*) heapTop + reqMem) >= (char*) heapEnd) {
		if (reqMem > HEAP_BLOCK_SIZE) {
			wprintf(L"--- Too long token value\n");
			::exit(1);
		}
		CreateHeapBlock();
	}
	t->val = (wchar_t*) heapTop;
	heapTop = (void*) ((char*) heapTop + reqMem);

	wcsncpy(t->val, tval, tlen);
	t->val[tlen] = L'\0';
}


Token* Scanner::NextToken() {
	while (ch == ' ' ||
			ch <= 31 || (ch >= 127 && ch <= 65535)
	) NextCh();
	if ((ch == L'/' && Comment0()) || (ch == L'/' && Comment1())) return NextToken();
	t = CreateToken();
	t->pos = pos; t->col = col; t->line = line;
	int state = start.state(ch);
	tlen = 0; AddCh();

	switch (state) {
		case -1: { t->kind = eofSym; break; } // NextCh already done
		case 0: { t->kind = noSym; break; }   // NextCh already done
		case 1:
			case_1:
			if ((ch >= L'0' && ch <= L':') || (ch >= L'A' && ch <= L'Z') || ch == L'_' || (ch >= L'a' && ch <= L'z')) {AddCh(); goto case_1;}
			else {t->kind = 1; break;}
		case 2:
			case_2:
			if (ch <= 9 || (ch >= 11 && ch <= 12) || (ch >= 14 && ch <= L'!') || (ch >= L'#' && ch <= L'[') || (ch >= L']' && ch <= 65535)) {AddCh(); goto case_2;}
			else if (ch == L'"') {AddCh(); goto case_4;}
			else if (ch == 92) {AddCh(); goto case_3;}
			else {t->kind = noSym; break;}
		case 3:
			case_3:
			if ((ch >= L' ' && ch <= L'~')) {AddCh(); goto case_2;}
			else {t->kind = noSym; break;}
		case 4:
			case_4:
			{t->kind = 2; break;}
		case 5:
			if ((ch >= L'A' && ch <= L'Z') || (ch >= L'a' && ch <= L'z')) {AddCh(); goto case_6;}
			else {t->kind = noSym; break;}
		case 6:
			case_6:
			if ((ch >= L'0' && ch <= L':') || (ch >= L'A' && ch <= L'Z') || ch == L'_' || (ch >= L'a' && ch <= L'z')) {AddCh(); goto case_6;}
			else {t->kind = 3; break;}
		case 7:
			case_7:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_8;}
			else {t->kind = noSym; break;}
		case 8:
			case_8:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_8;}
			else {t->kind = 4; break;}
		case 9:
			case_9:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_9;}
			else if (ch == L'E' || ch == L'e') {AddCh(); goto case_10;}
			else if (ch == L'.') {AddCh(); goto case_13;}
			else {t->kind = 4; break;}
		case 10:
			case_10:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_12;}
			else if (ch == L'+' || ch == L'-') {AddCh(); goto case_11;}
			else {t->kind = noSym; break;}
		case 11:
			case_11:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_12;}
			else {t->kind = noSym; break;}
		case 12:
			case_12:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_12;}
			else {t->kind = 4; break;}
		case 13:
			case_13:
			if ((ch >= L'0' && ch <= L'9')) {AddCh(); goto case_13;}
			else if (ch == L'E' || ch == L'e') {AddCh(); goto case_10;}
			else {t->kind = 4; break;}
		case 14:
			{t->kind = 5; break;}
		case 15:
			{t->kind = 6; break;}
		case 16:
			{t->kind = 9; break;}
		case 17:
			{t->kind = 10; break;}
		case 18:
			{t->kind = 11; break;}
		case 19:
			{t->kind = 12; break;}
		case 20:
			if (ch == L'.') {AddCh(); goto case_7;}
			else {t->kind = 8; break;}
		case 21:
			if (ch == L'.') {AddCh(); goto case_7;}
			else {t->kind = 7; break;}

	}
	AppendVal(t);
	return t;
}


// get the next token (possibly a token already seen during peeking)
Token* Scanner::Scan() {
	if (tokens->next == NULL) {
		return pt = tokens = NextToken();
	} else {
		pt = tokens = tokens->next;
		return tokens;
	}
}


// peek for the next token, ignore pragmas
Token* Scanner::Peek() {
	do {
		if (pt->next == NULL) {
			pt->next = NextToken();
		}
		pt = pt->next;
	} while (pt->kind > maxT); // skip pragmas

	return pt;
}


// make sure that peeking starts at the current scan position
void Scanner::ResetPeek() {
	pt = tokens;
}


} // namespace
} // namespace
} // namespace


