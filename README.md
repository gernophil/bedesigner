# bedesigner
Design guides using a reference genome

## Manual
You need to run this script in a Python environment containing `pandas` and `pyfaidx`.

To run the tool run this in an Pandas and pyfaidx containing Python environment:
`python bedesigner.py -g [path_to_reference_genome] -i [path_to_input_file] -o [path_to_output_file]`

The input file needs to be a csv file with `","` as separator containing the columns `variant` with the following format of variants: [chromosome]\_[genomic_position]\_[alternative_base]

For example: chr1_123456_T
