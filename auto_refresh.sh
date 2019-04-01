#!/bin/bash

if git pull | grep -q "Already up to date"; then
    echo "auto_refresh.sh: No need to update."
else
    git pull
    sudo make install > auto_refresh.out 2> auto_refresh.err
fi
