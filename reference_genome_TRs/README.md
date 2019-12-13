## Explaination for the files

- `code/tr_util.R` contains util functions for finding flanking microhomologies at ends of TR units
- `code/genome_util_lite.R` contains function to build a handy genome clousure
- `code/plot.R` conteins the code for generating plots in `plot/`
- `data` contains the reference genome assembly in `fasta` format
- `results/*.out` are output files directly from `trf` the [tandem repeat finder](https://tandem.bu.edu/trf/trf.html)
- `results/pars` contains the parameters used for the `trf` run
- `results/extract.awk` was used to extract relevent fields from the output files
- and `results/*.txt` are reformatted files output from `extract.awk`