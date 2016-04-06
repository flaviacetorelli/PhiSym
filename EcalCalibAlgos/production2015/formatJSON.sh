#!/bin/bash

cat ${1} | sed ':a;N;$!ba;s/\n/ /g' | sed 's:  ::g' | sed 's:,":\n":g' | sed 's:\[ \(.*\) \]:\[\1\]:g' 
