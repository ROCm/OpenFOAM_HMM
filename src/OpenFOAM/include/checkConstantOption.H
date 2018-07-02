// Unless -constant is present, skip startTime if it is "constant"

if
(
    !args.found("constant")
 && (startTime < Times.size()-1)
 && (Times[startTime].name() == "constant")
)
{
    ++startTime;
}
