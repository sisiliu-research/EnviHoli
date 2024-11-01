#!/usr/bin/env python3
# Python script to remove duplicate records from a SAM file based on read names

import sys

# Check if the correct number of arguments is provided
if len(sys.argv) != 3:
    print("Usage: python remove_duplicates.py input.sam output.sam")
    sys.exit(1)

# Define input and output file names
input_sam_filename = sys.argv[1]
output_sam_filename = sys.argv[2]

# Create a set to track unique SN: values
unique_sn_values = set()

# Flag to indicate if the current line is part of the @SQ section
is_sq_section = False

# Open the input SAM file for reading
with open(input_sam_filename, "r") as input_sam:
    # Open the output SAM file for writing
    with open(output_sam_filename, "w") as output_sam:
        # Iterate through the lines in the SAM file
        for line in input_sam:
            # Check if the line starts with @SQ
            if line.startswith("@SQ"):
                fields = line.strip().split("\t")
                sequence_name = None
                # Extract the SN (sequence name) field
                for field in fields:
                    if field.startswith("SN:"):
                        sequence_name = field
                if sequence_name:
                    # Check if the sequence name is not already in the set (i.e., it's unique)
                    if sequence_name not in unique_sn_values:
                        unique_sn_values.add(sequence_name)
                        output_sam.write(line)  # Write the unique @SQ line
            else:
                # If the line is not an @SQ line, it's part of the data section
                if not is_sq_section:
                    is_sq_section = True  # Enter the data section
                
                # Write the current line to the output file (data section)
                output_sam.write(line)

print("SAM file with duplicate @SQ lines removed written to", output_sam_filename)




