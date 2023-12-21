# bedesigner
Design guides using a reference genome

## Manual
You need to run this script in a Python environment containing `pandas` and `pyfaidx`.

To run the tool run this in an Pandas containing Python environment:
`python bedesigner.py [path_to_reference_genome] [path_to_input_file] [path_to_output_file]`

The input file needs to be a csv file with `","` as separator containing the columns `variant` with the following format of variants: [chromosome]\_[genomic_position]\_[alternative_base]

For example: chr1_123456_T

The old version of this tool required a compiled version of [this project](https://github.com/dancooke/bioio). Since using `pyfaidx` that's not necessary anymore. Now, it should also work, if the genome is not indexed and with chromosome names that are not just numbers (without "chr"), but I haven't tested that yet. For the last point (`chr1` or `1`) you need to disable line 37 of the script, since it currently simply removes any `chr` from the chromosome part, but I am working on an improvement.
