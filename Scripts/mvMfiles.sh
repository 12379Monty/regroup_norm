#!/bin/bash

Mfiles="`ls M*.*`"
for f in $Mfiles
do 
  mv "$f"  "_$f"
done
