#!/bin/bash

# first argument $1 should be rest of file directory to source location of fastq files NO FINAL SLASH
# second argument $2 should be rest of file directory to final location of fq files NO FINAL SLASH
# third argument $3 is the file.csv containing barcodes and filenames

# Define the base directory containing your fastq files from the Nanopore
BASE_DIR="$1"
# Define the destination directory for the concatenated files
DEST_DIR="$2"

# concatenate files, rename to the actual file name, move to be within the minimap2 destination directory

# Set the input txt file
barcode_file="$3"

# Read the txt file line by line
while IFS=$'\t' read -r barcode filename
do
    # Skip empty lines
    [ -z "$barcode" ] || [ -z "$filename" ] || {
        # Call your original command here; this is just an example
        echo "Processing barcode: $barcode with filename: $filename"
        
        # Navigate to the directory corresponding to the barcode
    DIR="$BASE_DIR/$barcode"
    
    # Check if the directory exists
    if [ -d "$DIR" ]; then

        # Navigate to the directory
        cd "$DIR"

        # Unzip all .gz files
        gunzip *.gz

        echo "Files for $barcode $filename unzipped"

        # Determine the output filename for this barcode
        output_fastq="$filename.fq"

        # Concatenate all .fastq files to the specified output filename
        cat *.fastq >> "$output_fastq"

        # Move the concatenated file to the destination directory
        cp "$BASE_DIR/$barcode/$output_fastq" "$DEST_DIR/"

        # Generate MD5 checksum for the concatenated file in the destination directory
        md5_dest=$(md5sum "${DEST_DIR}/${output_fastq}" | awk '{ print $1 }')
        md5_source=$(md5sum "${DIR}/${output_fastq}" | awk '{ print $1 }')

        # Compare the MD5 checksums
        if [ "$md5_dest" == "$md5_source" ]; then
            echo "MD5 checksum for $output_fastq matches between source and destination."
        else
            echo "WARNING: MD5 checksum for $output_fastq does NOT match between source and destination!"
        fi

        echo "Finished processing $barcode"
    else
        echo "Directory $DIR does not exist, skipping..."
    fi
    }
done < "$barcode_file"

echo "All barcodes in $1 have been concatenated and moved"