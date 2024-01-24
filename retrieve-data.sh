#!/bin/bash
JGI_SESSION_TOKEN=""

# Download genome fasta file
curl --cookie jgi_session=$JGI_SESSION_TOKEN \
     --output download.20240122.160059.zip \
     -d "{\"ids\":{\"Phytozome-200\":[\"52b9c7b9166e730e43a34f15\"]}}" \
     -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/

# Unzip and move files to proper directory
unzip download.20240122.160059.zip
unzip download.20240122.160059/**/*.zip
mv download.20240122.160059/Phytozome/PhytozomeV9/Lusitatissimum/assembly/Lusitatissimum_200.fa \
   data/

