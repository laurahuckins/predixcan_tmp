#!/bin/bash

for FILE in res* 

do 

  echo $FILE

  r2=$(head -7 $FILE| tail -1) 
  p=$(head -8 $FILE | tail -1)
  echo $FILE $r2 $p >> Res_Full
  
done

cat Res_Full | awk 's/res_//g;s/_1mb//g' | awk '$2>=0.01 || $3<=0.05' > Res_pass

