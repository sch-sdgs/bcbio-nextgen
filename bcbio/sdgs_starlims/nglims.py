
import collections
import glob
import gzip
import operator
import os
import subprocess
import joblib
import yaml
import json

from bcbio.log import logger
from bcbio import utils
from bcbio.distributed.transaction import file_transaction
from bcbio.sdgs_starlims.api import StarLimsApi
from bcbio.illumina import flowcell
from bcbio.pipeline.run_info import clean_name
from bcbio.workflow import template

json_data=open("/home/bioinfo/mparker/barcodes.json").read()

barcodes = json.loads(json_data)

# TODO: get new fields in starlims - genome, barcode? number of lanes?
# TODO: samples can appear twice if there has been a reanalysis - but no way of knowing? we need to have a reanalysis field

def prep_samples_and_config(run_folder, ldetails, fastq_dir, config):
    """Prepare sample fastq files and provide global sample configuration for the flowcell.

    Handles merging of fastq files split by lane and also by the bcl2fastq
    preparation process.
    """
    fastq_final_dir = utils.safe_makedir(os.path.join(fastq_dir, "merged"))

    cores = utils.get_in(config, ("algorithm", "num_cores"), 1)
    ldetails = joblib.Parallel(cores)(joblib.delayed(_prep_sample_and_config)(x, fastq_dir, fastq_final_dir)
                                      for x in _group_same_samples(ldetails))
    config_file = _write_sample_config(run_folder, [x for x in ldetails if x])
    return config_file, fastq_final_dir

def _prep_sample_and_config(ldetail_group, fastq_dir, fastq_final_dir):
    """Prepare output fastq file and configuration for a single sample.

    Only passes non-empty files through for processing.
    """
    files = []
    print "->", ldetail_group[0]["name"], len(ldetail_group)
    for read in ["R1", "R2"]:
        fastq_inputs = sorted(list(set(reduce(operator.add,
                                              (_get_fastq_files(x, read, fastq_dir) for x in ldetail_group)))))
        if len(fastq_inputs) > 0:
            logger.info("Merging FastQ files for" + ldetail_group[0]["name"])
            files.append(_concat_bgzip_fastq(fastq_inputs, fastq_final_dir, read, ldetail_group[0]))
        else:
            # TODO: need to fail if number of fastq files does not match expected
            return False
    if len(files) > 0:
        if _non_empty(files[0]):
            out = ldetail_group[0]
            out["files"] = files
            return out

def _non_empty(f):
    with gzip.open(f) as in_handle:
        for line in in_handle:
            return True
    return False


def _prepare_sample(data, run_folder):
    """Extract passed keywords from input LIMS information.
    """
    want = set(["description", "files", "genome_build", "name", "analysis", "upload", "algorithm"])
    out = {}
    for k, v in data.items():
        if k in want:
            out[k] = _relative_paths(v, run_folder)
    if "algorithm" not in out:
        analysis, algorithm = _select_default_algorithm(out.get("analysis"))
        out["algorithm"] = algorithm
        out["analysis"] = analysis
    description = "%s-%s" % (out["name"], clean_name(out["description"]))
    out["name"] = [out["name"], description]
    out["description"] = description
    return out

def _write_sample_config(run_folder, ldetails):
    """Generate a bcbio-nextgen YAML configuration file for processing a sample.
    """
    out_file = os.path.join(run_folder, "%s.yaml" % os.path.basename(run_folder))
    with open(out_file, "w") as out_handle:
        fc_name, fc_date = flowcell.parse_dirname(run_folder)
        out = {"details": sorted([_prepare_sample(x, run_folder) for x in ldetails],
                                 key=operator.itemgetter("name", "description")),
               "fc_name": fc_name,
               "fc_date": fc_date}
        yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)
    return out_file


def _relative_paths(xs, base_path):
    """Adjust paths to be relative to the provided base path.
    """
    if isinstance(xs, basestring):
        if xs.startswith(base_path):
            return xs.replace(base_path + "/", "", 1)
        else:
            return xs
    elif isinstance(xs, (list, tuple)):
        return [_relative_paths(x, base_path) for x in xs]
    elif isinstance(xs, dict):
        out = {}
        for k, v in xs.items():
            out[k] = _relative_paths(v, base_path)
        return out
    else:
        return xs


def _get_fastq_files(ldetail, read, fastq_dir):
    """Retrieve fastq files corresponding to the sample and read number.
    """
    print "READ " + read

    return glob.glob(os.path.join(fastq_dir, "Project_%s" % ldetail["project_name"],
                                  "Sample_%s" % ldetail["name"],
                                  "%s_*_%s_*.fastq.gz" % (ldetail["name"], read)))

def _concat_bgzip_fastq(finputs, out_dir, read, ldetail):
    """Concatenate multiple input fastq files, preparing a bgzipped output file.
    """
    out_file = os.path.join(out_dir, "%s_%s.fastq.gz" % (ldetail["name"], read))

    print "INPUTS: " + str(finputs)
    print "OUTPUT: " + out_file

    if not utils.file_exists(out_file):
        with file_transaction(out_file) as tx_out_file:
            subprocess.check_call("zcat %s | bgzip -c > %s" % (" ".join(finputs), tx_out_file), shell=True)
    return out_file

def _group_same_samples(ldetails):
    """Move samples into groups -- same groups have identical names.
    """
    sample_groups = collections.defaultdict(list)
    for ldetail in ldetails:
        sample_groups[ldetail["name"]].append(ldetail)


    for i in sorted(sample_groups.values(), key=lambda xs: xs[0]["name"]):
        print i

    return sorted(sample_groups.values(), key=lambda xs: xs[0]["name"])


def _select_default_algorithm(analysis):
    """Provide default algorithm sections from templates or standard
    """
    if not analysis or analysis == "Standard":
        return "Standard", {"aligner": "bwa", "platform": "illumina", "quality_format": "Standard",
                            "recalibrate": False, "realign": False, "mark_duplicates": True,
                            "variantcaller": False}
    elif "variant" in analysis:
        try:
            config, _ = template.name_to_config(analysis)
        except ValueError:
            config, _ = template.name_to_config("freebayes-variant")
        return "variant", config["details"][0]["algorithm"]
    else:
        return analysis, {}

def get_runinfo(run_folder, storedir):
    """Retrieve flattened run information for a processed directory from starLIMS webservices
    """

    star = StarLimsApi()
    run_id = os.path.basename(run_folder)
    starlims_info = star.run_details(run_id)
    fc_name, fc_date = flowcell.parse_dirname(run_folder)
    out = []
    lanes = []
    if "D00461" in run_folder:
        lanes = [1,2]
    if "M00969" in run_folder:
        lanes = [1]
    for lane in lanes:
        samples = []
        for i in starlims_info:
            if "RunNumber" in i:
                worklist = i["RunNumber"].replace("H", "")
                item = {}
                sample_id = i["CONTAINERID"]
                position = i["POSITION_IN_RUN"]
                if sample_id not in samples:
                    samples.append(sample_id)
                    if position is not None:
                        item["name"] = sample_id
                        position = position.lstrip('0')
                        bc_index = barcodes[str(position)]["sequence"]
                        item["regions"] = i["BED"]
                        item["analysis"] = "variant"
                        item["bc_index"] = bc_index
                        item["lane"] = lane  # lane number
                        item["name"] = sample_id
                        item["genome_build"] = "GRCh37"  # this is the genome build
                        item["description"] = "description"
                        item["project_name"] = worklist  # worklist
                        item["upload"] = {"method": "stalims", "run_id": fc_name,
                                          "fc_name": fc_name, "fc_date": fc_date,
                                          "dir": storedir, "library_method": i["CAPTUREMETHOD"]}

                        for k in ["lab_association", "private_libs", "researcher", "researcher_id", "sample_id",
                                  "galaxy_library", "galaxy_role"]:
                            item["upload"][k] = item.pop(k, "")

                        out.append(item)
    return out