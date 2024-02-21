#!/bin/sh
for name in */*; do
    [ ! -f "$name" ] && continue

    dir="$( basename "$( dirname "$name" )" )"
    newname="$dir/${dir}_${name##*/}"

    if [ ! -e "$newname" ]; then
        mv "$name" "$newname"
    fi
done