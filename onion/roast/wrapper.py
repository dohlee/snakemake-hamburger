__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get('extra', '')

# Execute shell command.
shell("sed '1 s/^/roasted_/g' {snakemake.input.onion} > {snakemake.output[0]}")
