# bedesigner
Design guides using a reference genome

## Manual
You need to run this script in a Python environment containing `pandas` and `pyfaidx`.
```
conda env create -f environment.yaml # will create a conda env with pandas and pyfaidx called bedesigner
pip install -r requirements.txt # will isntall pandas and pyfaidx using PyPI (highly recommended to do this in a virtual environment)
```

To run the tool run this in an Pandas and pyfaidx containing Python environment:
`python bedesigner.py -g path/to/reference-genome.fasta -i path/to/input-file.csv -o path/to/output-file.csv`

The input file needs to be a csv file with `","` as separator containing the columns `variant` with the following format of variants: [chromosome]\_[genomic_position]\_[alternative_base]

For example: chr1_123456_T


```
usage: python bedesigner.py [-h] -i <PATH> -g <PATH> -o <PATH> [-p {NGG,NG}] [-s WINDOW_START] [-e WINDOW_END] [-l GUIDE_LENGTH] [-n IGNORE_STRING] [-b {both,ABE,CBE}] [-a]

bedesigner

options:
  -h, --help            show this help message and exit
  -p {NGG,NG}, --pam-site {NGG,NG}
                        Sequence of the PAM site
  -s WINDOW_START, --window-start WINDOW_START
                        Starting position of editing window
  -e WINDOW_END, --window-end WINDOW_END
                        End position of editing window
  -l GUIDE_LENGTH, --guide-length GUIDE_LENGTH
                        Length of guide
  -n IGNORE_STRING, --ignore-string IGNORE_STRING
                        Substring to be ignored in the variant string
  -b {both,ABE,CBE}, --base-editor {both,ABE,CBE}
                        Base editor to design guides for
  -a, --all-possible    Design for all ALTs

required arguments:
  -i <PATH>, --input <PATH>
                        Path to the input file
  -g <PATH>, --ref-genome <PATH>
                        Path to the reference genome
  -o <PATH>, --output <PATH>
                        Path to the output file
```
