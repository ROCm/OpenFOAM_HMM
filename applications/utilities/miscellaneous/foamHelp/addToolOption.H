argList::addArgument("tool");
const wordList opts(helpType::dictionaryConstructorTablePtr_->sortedToc());

string note = "Valid <tool> options include:";
forAll(opts, i)
{
    note = note + ' ' + opts[i];
}

argList::notes.append(note);
