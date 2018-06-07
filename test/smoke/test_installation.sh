#!/bin/bash

echo ""
echo "---------------------------------------------------------------------------------------"

VERSION="$(python -c "from main.version import __version__; print __version__")"

red=`tput setaf 1`
green=`tput setaf 2`
reset=`tput sgr0`

echo "Smoke test of installation:"

./cart2enst --help &> /dev/null

if [ $? != 0 ]
then
    echo "${red}Error: CARTtools $VERSION installation failed.${reset}"
    echo "---------------------------------------------------------------------------------------"
    echo ""
    exit 1
fi


echo "${green}CARTtools $VERSION installation succeeded!${reset}"
echo "---------------------------------------------------------------------------------------"
echo ""
