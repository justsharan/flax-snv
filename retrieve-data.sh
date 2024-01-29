#!/bin/bash
JGI_SESSION_TOKEN=""
FILENAME="output.zip"

# Download genome fasta file
curl --cookie jgi_session=$JGI_SESSION_TOKEN \
     --output $FILENAME \
     -d "{\"ids\":{\"Phytozome-200\":[\"53112aae49607a1be005599c\",\"53112ab149607a1be00559a3\"]}}" \
     -H "Content-Type: application/json" https://files-download.jgi.doe.gov/filedownload/

# Unzip and move files to proper directory
unzip $FILENAME
