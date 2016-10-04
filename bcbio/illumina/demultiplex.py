"""Demultiplex and fastq conversion from Illumina output directories.

Uses Illumina's bcl2fastq: http://support.illumina.com/downloads/bcl2fastq_conversion_software_184.ilmn
"""
import os
import subprocess
import time
from bcbio.log import logger

from bcbio import utils


def run_bcl2fastq(run_folder, ss_csv, config):
    """Run bcl2fastq for de-multiplexing and fastq generation.
    run_folder -- directory of Illumina outputs
    ss_csv -- Samplesheet CSV file describing samples.
    """
    bc_dir = os.path.join(run_folder, "Data", "Intensities", "BaseCalls")
    #output_dir = os.path.join(run_folder, "fastq")
    output_dir= utils.get_in(config, ("process", "storedir"))
    if not os.path.exists(os.path.join(output_dir, "Makefile")):
        cmd = " ".join(["configureBclToFastq.pl", "--use-bases-mask", "y*,I8,y*", "--mismatches", "1", "--fastq-cluster-count",
             "0", "--no-eamss", "--input-dir", bc_dir, "--output-dir", output_dir,
             "--adapter-sequence",
             "/usr/local/bcl2fastq-1.8.4/share/bcl2fastq-1.8.4/adapters/TruSeq_r1.fa",
             "--adapter-sequence",
             "/usr/local/bcl2fastq-1.8.4/share/bcl2fastq-1.8.4/adapters/TruSeq_r2.fa",
             "--sample-sheet", ss_csv, "--force"])
        print cmd
        logger.info("Configuring BclToFastq")
        logger.info(cmd)

        subprocess.check_call(cmd)

    with utils.chdir(output_dir):
        cores = str(utils.get_in(config, ("algorithm", "num_cores"), 1))
        cmd = ["make", "-j", cores]
        if "submit_cmd" in config["process"] and "bcl2fastq_batch" in config["process"]:
            _submit_and_wait(cmd, cores, config, output_dir)
        else:
            subprocess.check_call(cmd)
    return output_dir


def _submit_and_wait(cmd, cores, config, output_dir):
    """Submit command with batch script specified in configuration, wait until finished
    """

    batch_script = "submit_bcl2fastq.sh"
    if not os.path.exists(batch_script + ".finished"):
        if os.path.exists(batch_script + ".failed"):
            os.remove(batch_script + ".failed")
        with open(batch_script, "w") as out_handle:
            out_handle.write(config["process"]["bcl2fastq_batch"].format(
                cores=cores, bcl2fastq_cmd=" ".join(cmd), batch_script=batch_script, output_dir=output_dir))
        submit_cmd = utils.get_in(config, ("process", "submit_cmd"))
        subprocess.check_call(submit_cmd.format(batch_script=batch_script), shell=True)
        # wait until finished or failure checkpoint file
        while 1:
            if os.path.exists(batch_script + ".finished"):
                break
            if os.path.exists(batch_script + ".failed"):
                raise ValueError("bcl2fastq batch script failed: %s" %
                                 os.path.join(output_dir, batch_script))
            time.sleep(5)
