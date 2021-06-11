# seqStats

This Perl program calculates basic summaries of fasta sequence files.


INSTALLATION

*seqStats* has no dependencies other than Perl itself.


USAGE

Program usage is as follows:

```
seqStats FASTA_input_file
```
```
Optional flags:
-c|contig INT     limit output to contigs larger than INT (default: 500 bp)
-d|distribution   return a file listing the sizes of all input sequences
```


EXAMPLE

```
seqStats test.fna -c 4 -d
```
