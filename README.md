# Djerba - PanCuRx

![Djerba](./doc/djerba_logo_small.png)

Modular system to create reports from metadata and workflow output, see more at [ReadTheDocs](https://djerba.readthedocs.io/en/latest/).

The software in this directory contains plugins to djerba specific to PanCuRx and LBR at OICR.

## Loading the Environment

To enter the djerba environment on the OICR cluster, start with:

```
qrsh -P pcsi -l h_vmem=10G

module load djerba

export DJERBA_PRIVATE_DIR=/.mounts/labs/PCSI/users/
export DJERBA_SOURCE_DIR=/.mounts/labs/PCSI/users/fbeaudry/djerba-pancurx 
```

## Building the INI
Next, you'll need to make .INI file that contains the parameters of your djerba run. At the very least, the minimal INI should look like:

```
[core]

[sample_helper]
donor = BTC0025
tumour=BTC_0025_Sn_M_526
normal=BTC_0025_Ly_R

[file_helper]
[pancurx.summary]
[pancurx.somatic]
[pancurx.germline]
[pancurx.classification]
[pancurx.all_genes]
[pancurx.appendix]
[pancurx.slide]
```

where square brackets (e.g. `[pancurx.summary]`) indicate the specific plugins to run. If you do not wish to print the entire djerba report, it is possible to run djerba with only a subset of the plugins. 

Within a pluging, parameters are set by naming the parameter, putting an `=` sign and then followed by the value for the parameter. In the above, all values are left default except in the `sample_helper` plugin which requires the donor id, tumour sample id, and normal id. Many more parameters can be specified when necessary.

## Launching djerba
Finally, you can launch djerba with the following command to first make the report (`report`), then make a .pdf of the report (`render`):

```
python3 ${DJERBA_SOURCE_DIR}/src/bin/djerba.py --debug report --ini config.ini --out-dir report --no-archive

python3 ${DJERBA_SOURCE_DIR}/src/bin/djerba.py render --json report/report.json --out-dir report --pdf --no-archive
```

## Copyright and License

Copyright &copy; 2020-2024 by Genome Sequence Informatics, Ontario Institute for Cancer Research.

Licensed under the [GPL 3.0 license](https://www.gnu.org/licenses/gpl-3.0.en.html).
