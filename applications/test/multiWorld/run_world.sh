#!/bin/sh
world="${1:?specify world/case}"
application="${2:-laplacianFoam}"

# case "$application" in
# (-test)
#     application="Test-multiWorld1"
#     ;;
# esac

if [ -z "$application" ] || ! command -v "$application" > /dev/null
then
    echo "No application: $application"
    exit 2
fi

worldCase="$(echo "$world" | tr '[:upper:]' '[:lower:]')"
worldName="$(echo "$world" | tr '[:lower:]' '[:upper:]')"

"$application" -case "$worldCase" -world "$worldName" 2>&1 | tee log.run_"$worldCase"
read dummy

# ----------------------------------------------------------------------------
