#!/bin/sh

mkdir -p out condor-out

mkdir -p lists

cd lists

CreateFileList.pl -run 4 -type 11 -n 100000 -embed DST_BBC_G4HIT DST_CALO_CLUSTER DST_CALO_G4HIT DST_TRUTH_JET DST_VERTEX

split -l40 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".list" dst_truth_jet.list "dst_truth_jet_"
split -l40 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".list" dst_calo_g4hit.list "dst_calo_g4hit_"
split -l40 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".list" dst_calo_cluster.list "dst_calo_cluster_"
split -l40 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".list" dst_vertex.list "dst_vertex_"
split -l40 --numeric-suffixes=1 --suffix-length=2 --additional-suffix=".list" dst_bbc_g4hit.list "dst_bbc_g4hit_"

if [ -f index.txt ]
then
    rm index.txt
fi

for i in {1..75..1}
do
  printf "%0*d\n" 2 $i >> index.txt
done
