__author__ = "Dohoon Lee"
__copyright__ = "Copyright 2018, Dohoon Lee"
__email__ = "dohlee.bioinfo@gmail.com"
__license__ = "MIT"


from os import path

from snakemake.shell import shell

# Extract arguments.
extra = snakemake.params.get('extra', '')

def center_aligned_text(file_path, width):
    return open(file_path).read().strip().center(width)

file_paths = [
    snakemake.input.warm_bun,
    snakemake.input.pickle,
    snakemake.input.roasted_onion,
    snakemake.input.tomato,
    snakemake.input.cheese,
    snakemake.input.patty,
    snakemake.input.lettuce,
    snakemake.input.warm_bun
]

hamburger = '\n'.join([center_aligned_text(file_path, 20) for file_path in file_paths])

# Execute shell command.
shell("echo '" +
      hamburger +
      "' > {snakemake.output[0]} && "
      "rm {snakemake.input.warm_bun} "
      "{snakemake.input.pickle} "
      "{snakemake.input.roasted_onion} "
      "{snakemake.input.tomato} "
      "{snakemake.input.cheese} "
      "{snakemake.input.patty} "
      "{snakemake.input.lettuce}")


