#!/bin/sh

mkdir -p out condor-out/jet_bg out/jet_bg condor-out/calo_rho out/calo_rho

mkdir -p lists lists/calo_rho lists/jet_bg

DIR=$(cd "$(dirname "$0")"; pwd)

cd lists/calo_rho

echo "entering directory `pwd`"

CreateFileList.pl -run 4 -type 11 -embed DST_BBC_G4HIT DST_CALO_CLUSTER DST_CALO_G4HIT DST_TRUTH_JET DST_VERTEX

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

cd $DIR

cd lists/jet_bg

echo "entering directory `pwd`"

CreateFileList.pl -run 4 -type 4 DST_BBC_G4HIT DST_CALO_CLUSTER DST_CALO_G4HIT DST_VERTEX

split -l100 --numeric-suffixes=1 --suffix-length=3 --additional-suffix=".list" dst_calo_g4hit.list "dst_calo_g4hit_"
split -l100 --numeric-suffixes=1 --suffix-length=3 --additional-suffix=".list" dst_calo_cluster.list "dst_calo_cluster_"
split -l100 --numeric-suffixes=1 --suffix-length=3 --additional-suffix=".list" dst_vertex.list "dst_vertex_"
split -l100 --numeric-suffixes=1 --suffix-length=3 --additional-suffix=".list" dst_bbc_g4hit.list "dst_bbc_g4hit_"

if [ -f index.txt ]
then
    rm index.txt
fi

for i in {1..400..1}
do
  printf "%0*d\n" 3 $i >> index.txt
done
