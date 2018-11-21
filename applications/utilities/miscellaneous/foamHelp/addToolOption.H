argList::addArgument("tool");

argList::notes.append("Valid <tool> options include:");
for (const word& tool : helpType::dictionaryConstructorTablePtr_->sortedToc())
{
    argList::notes.append("    " + tool);
}

argList::notes.append
(
    "\nNOTE the <tool> must actually appear *before* any options"
);
