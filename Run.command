#!/bin/bash

cd "${0%/*}"

if [ -f pid.txt ]; then
	p_id=`cat pid.txt`
	kill -9 $p_id
	rm pid.txt  
fi

cd src
Rscript -e 'library(methods); shiny::runApp(launch.browser=TRUE)'
