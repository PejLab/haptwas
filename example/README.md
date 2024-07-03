# Command line examples


Here we provide fictitious data to demonstrate the
use of `afcn` from the command line.  The sections
are divided by `afcn` subcommands: `fit`, `predict`,
and `twas`.


## fit


🚧 Under Construction 🚧


## predict

Let's assume that we are in the directory consisting
of files:

```
- .
    |- parameters.bed
    |- simulated.vcf.bgz
    |- simulated.vcf.bgz.tbi
```

where the *simulated.vcf.bgz.tbi* file is the `tabix` index 
of *simulated.vcf.bgz.tbi*.  


### Default inputs

`afcn predict` requires the specification of a **vcf** and parameter **bed**
file.  Given these filenames gene expression predictions are generated
by


```
$ python -m afcn predict --vcf simulated.vcf.bgz --params parameters.bed
```

Under the results are printed to a newly created directory *afcn_predict*

```
- afcn_predict
    |- output.log
    |- output_by_haplotype.bed
    |- output_total.bed
```

The *output_by_haplotype.bed* and *output_total.bed* files contain
the linear scale gene expression predictions for each haplotype or
in total for each sample, respectively.  The *output.log* always records the input
command, start time, and end time, if errors occur they will also
be written to the log file.


### Customize scale and rounding

Suppose I wanted the data to be printed in log base 2 scale and
each value rounded to two decimal places.  This can be accommodated by
`afcn predict` optional arguments

* `--scale` : linear (default), log, log2, or log10
* `--decimals` : No rounding(default), or postive integer
    that deteremines the decimal place that predictions
    are rounded.

The desired predictions are generated by the following command
```
python -m afcn predict --vcf simulated.vcf.bgz --params parameters.bed \
--scale log2 --decimals 2
```

Note, that the meta data of the prediction bed files records the
scale that predictions are computed and decimal place for which
they are rounded.

### Customizing outputs, scale, and rounding

Suppose that I want to record the natural log values of gene
expression predictions to three decimal places in directory
*log_transformed_gene_expr_predictions* and file prefix
*todays_predictions*.  This is done by specifying the options

* `-p`: string for file prefix, e.g. `-p today` results in a
    total gene expression file *today_predictions_total.bed*
* `-o`: string / path of directory to write prediction files


This can be produced by

```
python -m afcn predict --vcf simulated.vcf.bgz --params parameters.bed \
--scale log --decimals 3 \
-p todays_predictions -o log_transformed_gene_expr_predictions
```

which will generate

```
- log_transformed_gene_expr_predictions
    |- todays_predictions.log
    |- todays_predictions_by_hap.bed
    |- todays_predictions_out.bed
```



## twas

🚧 Under Construction 🚧