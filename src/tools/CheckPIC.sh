#!/bin/bash

for i in */*.o */*/*.o *.o */*.a *.a; do
  if readelf --relocs ${i} | awk '$3~/^R_/ && $5!~/^\.debug/{print $3}' | sort -u | grep "R_X86_64_32"; then
  	echo "${i} contains unrelocatable symbol."
  fi
done