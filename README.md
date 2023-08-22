# bedesigner
Design guides using a indexed reference genome

## manual
To run this script you need to compile [this project](https://github.com/dancooke/bioio) and copy the fasta binary into the root folder of this project. You might need to change the path to fasta.exe, if you are using Windows.

To run the tool run this in an Pandas containing Python environment:
`python bedesigner.py [path_to_reference_genome] [path_to_input_file] [path_to_output_file]`

The input file need to be a csv file with `","` as separator containing the columns `variant` with the following format of variants: [chromosome]_[genomic_position]_[alternative_base]

For example: chr1_123456_T

Currently it only works if the genome is indexed and the chromosome names are just numbers (without "chr"), but improvement is planed.
