#!/bin/bash

#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --job-name=zip_files
#SBATCH --mail-user=jeremy.rotzetter@students.unibe.ch
#SBATCH --mail-type=end
#SBATCH --output=/data/users/jrotzetter/ribosome-profiling/logs/output_zip_files_%j.o
#SBATCH --error=/data/users/jrotzetter/ribosome-profiling/logs/error_zip_files_%j.e
#SBATCH --partition=pall

# Script to zip a list of files specified in FILELIST. One line per file (also works with full paths to the files)

# Set the path to the file with the files to be zipped
FILELIST=$1 ## NOTE: there should be an empty line after the last file to be zipped otherwise said file will not be zipped

# Set the output directory path
OUTDIR=$2

# Set the output zip file name
OUTPUT_ZIP="$OUTDIR/$3.zip"

# Create the output directory if it doesn't exist
mkdir -p "$OUTDIR"

# Read the list of files from a file called FILELIST
while IFS= read -r file; do
  # Add each file to the zip file
  zip -rj9 "$OUTPUT_ZIP" "$file"
done < ${FILELIST}

# zip options:
# -r   recurse into directories
# -j   junk (don't record) directory name. Needed so the zipped files are not stored with their full path
# -9   compress better

# List the contents of the zip file
unzip -l "$OUTPUT_ZIP"

# Method to just gzip a list of files individually
# while IFS= read file; do
#     gzip -c "$file" > "${OUTDIR}/$(basename "$file").gz"
# done < ${FILELIST}

# How it works
# ...  < filelist redirects the contents of filelist to ....

# while IFS= read file; do ... done goes through the lines in filelist, stores the contents of the currently processed line in the variable file and executes ....

# IFS= modifies the internal file separator. This is needed to handle multiple, leading and trailing spaces properly.

# gzip -c "$file" > "OUTDIR/$(basename "$file").gz" compresses the currently processed file and stores the output in a file with the same name plus an .gz extension in the directory OUTDIR.

# Here basename "$file" extracts the bare filename from the file's path.
# Source: https://superuser.com/questions/577436/gzip-several-files-in-different-directories-and-copy-to-new-directory