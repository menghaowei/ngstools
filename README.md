![logo](./ngstools_logo.png)

# ngstools
My own tools code for NGS data analysis (Next Generation Sequencing)

## parse-mpileup.py
parse `samtools mpileup` command output, also known as `.pileup` file, the file like:

```
chr1	10030	c	1	^6.	A
chr1	10031	t	1	.	A
chr1	10032	a	1	.	F
chr1	10033	a	1	.	F
chr1	10034	c	1	.	<
chr1	10035	c	1	.	F
chr1	10036	c	0	*	*
chr1	10037	t	1	.-1A	A
chr1	10038	a	0	*	*
chr1	10039	a	0	*	*
```
The `.pileup` format explain please check the HTML 
[pileup explain](http://samtools.sourceforge.net/pileup.shtml)

And convert `.pileup` file into `.bmat`format, the format like:

```
chr_name	chr_index	ref_base	A	G	C	T	del_count	insert_count	ambiguous_count	deletioninsertion	ambiguous	mut_num
chr1	10030	C	0	0	1	0	0	0	0	.	.	.	0
chr1	10031	T	0	0	0	1	0	0	0	.	.	.	0
chr1	10032	A	1	0	0	0	0	0	0	.	.	.	0
chr1	10033	A	1	0	0	0	0	0	0	.	.	.	0
chr1	10034	C	0	0	1	0	0	0	0	.	.	.	0
chr1	10035	C	0	0	1	0	0	0	0	.	.	.	0
chr1	10036	C	0	0	0	0	1	0	0	*	.	.	0
chr1	10037	T	0	0	0	1	1	0	0	A	.	.	0
chr1	10038	A	0	0	0	0	1	0	0	*	.	.	0
```
For help info, please run `python parse-mpileup.py -h`:

```
python parse-mpileup.py  -h
usage: parse-mpileup.py [-h] -i INPUT [-o OUTPUT] [-p THREADS] [-n MUTNUM]
                        [--TempDir TEMPDIR]

convert mpileup file to info file

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --Input INPUT
                        samtools mpileup format file
  -o OUTPUT, --Output OUTPUT
                        Output parsed file
  -p THREADS, --Threads THREADS
                        Multiple threads number, default=1
  -n MUTNUM, --MutNum MUTNUM
                        Only contain mutation info go to the output, set 0
                        mean output all site, default=0
  --TempDir TEMPDIR     Where to keep temp files, default is the same dir with
                        --Input
```
