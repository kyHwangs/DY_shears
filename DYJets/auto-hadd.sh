#! /usr/bin/env bash

echo $1
cd $1/results/ || exit 1

for s in data DYJets; do 
    echo $s
    ls dyjets-$s-?.root dyjets-$s-??.root
    hadd -f dyjets-$s.root dyjets-$s-?.root dyjets-$s-??.root
done

#for s in TT WToLNu; do 
#    echo $s
#    hadd -f dyjets-$s.root dyjets-$s-?.root
#done

#hadd -f tmp.root dyjets-data.root dyjets-data-smu.root; 
#mv tmp.root dyjets-data.root

