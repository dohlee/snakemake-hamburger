__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get('extra', '')

# Execute shell command.
shell("echo '__p_ic___kle___' > {snakemake.output[0]}")
