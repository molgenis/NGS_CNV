#!/usr/bin/env python
import os
import argparse


def get_parameter_values():
    """Define CLI parameters and return their set values.

    Returns
    -------
    dict
        Set parameter values
    """
    make_hc_args = argparse.ArgumentParser()
    make_hc_args.add_argument("-b", "--bed-file", type=str, required=True, dest="bed-file", help="Path to BED file to use as starting point for the High Confident BED file")
    make_hc_args.add_argument("-m", "--minimum-samples", type=int, default=30, dest="minimum-samples", help="Minimum number of samples each population should have")
    make_hc_args.add_argument("-n", "--name-prefix", type=str, default="makehc", dest="name-prefix", help="Name prefix to use for the created script files")
    make_hc_args.add_argument("-o", "--output-directory", type=str, required=True, dest="output-directory", help="Path to write created scripts to")
    make_hc_args.add_argument("-p", "--project-directory", type=str, required=True, dest="project-directory", help="Path to directory to use for creating the High Confident BED file")
    make_hc_args.add_argument("-s", "--samples-to-use", type=str, required=True, dest="samples-to-use", help="Path to file with samples to use for the two populations")
    make_hc_args.add_argument("-u", "--umcu-scripts-dir", type=str, required=True, dest="umcu-scripts-dir", help="Path to directory containing the UMCU scripts required to create the High Confident BED file")
    make_hc_args.add_argument("-v", "--umcu-scripts-venv", type=str, required=True, dest="umcu-scripts-venv", help="Path to directory containing thge UMCU Python Virtual ENV")
    return vars(make_hc_args.parse_args())


def check_parameter_values(parameter_values):
    """Check that the set parameters are correct.

    Parameters
    ----------
    parameter_values : dict
        Set parameter values

    Returns
    -------
    incorrect_parameters : list of str
        Names of incorrect parameters
    """
    incorrect_parameters = []

    # Check that the project directory exists.
    if not os.path.isdir(parameter_values["project-directory"]):
        incorrect_parameters.append("-p / --project-directory")

    # Check that the BED file exists.
    if not os.path.isfile(parameter_values["bed-file"]):
        incorrect_parameters.append("-b / --bed-file")

    # Check that the samples list file exists.
    if not os.path.isfile(parameter_values["samples-to-use"]):
        incorrect_parameters.append("-s / --samples-to-use")

    # Check that the directory with UMCU scripts exists.
    if not os.path.isdir(parameter_values["umcu-scripts-dir"]):
        incorrect_parameters.append("-u / --umcu-scripts-dir")

    # Check that the directory with the UMCU Python Virtual ENV exists.
    if not os.path.isdir(parameter_values["umcu-scripts-venv"]):
        incorrect_parameters.append("-v / --umcu-scripts-venv")

    # Check that the output directory exists
    if not os.path.isdir(parameter_values["output-directory"]):
        incorrect_parameters.append("-o / --output-directory")

    return incorrect_parameters


def read_sample_list(sample_list_file_location):
    """Read the file with samples to use for the two make and two female populations

    Parameters
    ----------
    sample_list_file_location : str
        Path to the file with samples to use

    Returns
    -------
    population_samples : dict
        Samples for the four populations
    """
    population_labels = ["F1", "F2", "M1", "M2"]
    population_samples = {"F1": [], "F2": [], "M1": [], "M2":[]}
    try:
        with open(sample_list_file_location, 'r') as samplefile:
            for fileline in samplefile:
                filelinedata = fileline.strip().split("\t")
                if filelinedata[1] in population_labels:
                    population_samples[filelinedata[1]].append(filelinedata[0])
        return population_samples
    except IOError:
        print(f"Could not read sample list {sample_list_file_location}")


def check_populations(populationsamples, samplethreshold):
    """Make sure that the populations sizes are the same and have a set minimum number of samples.

    Parameters
    ----------
    populationsamples : dict
        Samples to use to construct the populations
    samplethreshold : int
        Minimum number of samples to use for

    Returns
    -------
    populationsamples : dict
        Modified population samples or empty dict
    """
    population_sizes = [len(populationsamples["F1"]), len(populationsamples["F2"]), len(populationsamples["M1"]), len(populationsamples["M2"])]
    if population_sizes[0] == population_sizes[1] == population_sizes[2] == population_sizes[3]:
        return populationsamples
    else:
        smallest_population = min(population_sizes)
        if smallest_population > samplethreshold:
            print(f"...Populations are not of equal size. Filtering populations to retain {smallest_population} samples...")
            populationsamples["F1"] = populationsamples["F1"][0:smallest_population]
            populationsamples["F2"] = populationsamples["F2"][0:smallest_population]
            populationsamples["M1"] = populationsamples["M1"][0:smallest_population]
            populationsamples["M2"] = populationsamples["M2"][0:smallest_population]
            return populationsamples
        else:
            print(f"...Populations smaller than set threshold {samplethreshold}...")
            return {}


def read_script_template(script_file_loc):
    """Read and return a script to use as a template.

    Parameters
    ----------
    script_file_loc : str
        Path to the script file used as a template

    Returns
    -------
    script_template : str
        Read script to use as a template
    """
    try:
        with open(script_file_loc, 'r') as templatefile:
            script_template = templatefile.read()
        return script_template
    except IOError:
        print(f"Could not read script template {script_file_loc}")


def modify_template(scripttemplate, projectdir, bedfile, slicedbedfile, umcuscripts, umcuvenv, outputdir):
    """Modify the template script to insert user supplied values.

    Parameters
    ----------
    scripttemplate : str
        Script template for creating the HighConfident BED file
    projectdir : str
        Path to directory to use for creating the High Confident BED file
    bedfile : str
        Path to BED file to start with
    slicedbedfile
        Path to sliced BED file
    umcuscripts : str
        Path to directory with UMCU Python scripts
    umcuvenv : str
        Path to directory with UMCU Python Virtual Env

    Returns
    -------
    hc_script1 : str
        Modified script for part1
    """
    hc_script = scripttemplate
    hc_script = hc_script.replace("\"${hcdir}\"", projectdir)
    hc_script = hc_script.replace("\"${bedfile}\"", bedfile)
    hc_script = hc_script.replace("\"${slicedbedfile}\"", slicedbedfile)
    hc_script = hc_script.replace("\"${umcuscriptsdir}\"", umcuscripts)
    hc_script = hc_script.replace("\"${umcuvenvdir}\"", umcuvenv)
    hc_script = hc_script.replace("\"${generatesambamba}\"", outputdir)
    return hc_script


def add_population_samples(hcscript, hcsamples, projectdir):
    """Add the population samples to the first script.

    Parameters
    ----------
    hcscript : str
        Script for the first part
    hcsamples : dict
        Samples to use for the four populations
    projectdir : str
        Path to the directory used for constructing the High Confident BED file

    Returns
    -------
    hcscript : str
        Modified scripts containing commands to soft link the samples
    """
    #Make the links
    f1_links = "\n".join([f"ln -s {x} {projectdir}population1/female/{x.split('/')[-1]}" for x in hcsamples["F1"]])
    f2_links = "\n".join([f"ln -s {x} {projectdir}population2/female/{x.split('/')[-1]}" for x in hcsamples["F2"]])
    m1_links = "\n".join([f"ln -s {x} {projectdir}population1/male/{x.split('/')[-1]}" for x in hcsamples["M1"]])
    m2_links = "\n".join([f"ln -s {x} {projectdir}population2/male/{x.split('/')[-1]}" for x in hcsamples["M2"]])

    # Add the links to the script
    hcscript = hcscript.replace("#<LINK F1>", f1_links)
    hcscript = hcscript.replace("#<LINK F2>", f2_links)
    hcscript = hcscript.replace("#<LINK M1>", m1_links)
    hcscript = hcscript.replace("#<LINK M2>", m2_links)
    return hcscript


def write_script(output_file_location, script_to_write):
    """Write the script to an output file.

    Parameters
    ----------
    output_file_location : str
        Path to write output file to
    script_to_write : str
        Script contents to write to the output file

    Returns
    -------
    file_written : bool
        True if output file has succesfully been written, False if not
    """
    file_written = False
    try:
        with open(output_file_location, 'w') as output_file:
            output_file.write(script_to_write)
        file_written = True
    except IOError:
        print(f"Could not write script file {output_file_location}")
    finally:
        return file_written


def copy_generate_sambamba(outdirloc):
    try:
        gsscript = ""
        with open("generate_sambamba.py", 'r') as gsoriginal:
            gsscript = gsoriginal.read()

        with open(f"{outdirloc}generate_sambamba.py", 'w') as gscopy:
            gscopy.write(gsscript)
    except IOError:
        print(f"Could not copy generate_sambamba.py to {outdirloc}")


def get_sliced_bed_file(projectdir, bedfile):
    """Make the file name and path for the sliced BED file.

    Parameters
    ----------
    projectdir : str
        Path to directory to use for creating the High Confident BED file
    bedfile : str
        Path to BED file to use

    Returns
    -------
    str
        Path to the sliced BED file
    """
    bedfilename = bedfile.split("/")[-1]
    return f"{projectdir}sliced_{bedfilename}"


def add_dir_slash(directory_location):
    """Adds a trailing / to a directory path.

    Parameters
    ----------
    directory_location : str
        Path to directory

    Returns
    -------
    directory_location : str
        Path to directory with trailing /
    """
    if not directory_location.endswith("/"):
        directory_location = f"{directory_location}/"
    return directory_location


def main():
    """Do the actual work."""
    make_hc_parameters = get_parameter_values()
    incorrectparameters = check_parameter_values(make_hc_parameters)

    if len(incorrectparameters) == 0:
        # Make sure paths to directories end in /
        project_directory = add_dir_slash(make_hc_parameters["project-directory"])
        umcu_scripts_directory = add_dir_slash(make_hc_parameters["umcu-scripts-dir"])
        umcu_venv_directory = add_dir_slash(make_hc_parameters["umcu-scripts-venv"])
        output_directory = add_dir_slash(make_hc_parameters["output-directory"])
        output_prefix = make_hc_parameters["name-prefix"]

        # Copy the generate_sambamba.py script to 
        copy_generate_sambamba(output_directory)

        hc_population_samples = read_sample_list(make_hc_parameters["samples-to-use"])
        hc_population_samples = check_populations(hc_population_samples, make_hc_parameters["minimum-samples"])

        if len(hc_population_samples) == 4:
            # Read the two template scripts
            part1_template = read_script_template("make_hcbedfile_part1.sh")
            part2_template = read_script_template("make_hcbedfile_part2.sh")

            # Get the path for the sliced BED file
            sliced_bed_file = get_sliced_bed_file(project_directory, make_hc_parameters["bed-file"])

            # Make the functional scripts
            part1_script = modify_template(part1_template, project_directory, make_hc_parameters["bed-file"], sliced_bed_file, umcu_scripts_directory, umcu_venv_directory, output_directory)
            part1_script = add_population_samples(part1_script, hc_population_samples, project_directory)
            part2_script = modify_template(part2_template, project_directory, make_hc_parameters["bed-file"], sliced_bed_file, umcu_scripts_directory, umcu_venv_directory, output_directory)

            script1_output_path = f"{output_directory}{output_prefix}_part1.sh"
            script2_output_path = f"{output_directory}{output_prefix}_part2.sh"

            # Write the functional scripts
            wrote_script1 = write_script(script1_output_path, part1_script)
            print(f"Wrote the script for the first part?: {wrote_script1}")
            wrote_script2 = write_script(script2_output_path, part2_script)
            print(f"Wrote the script for the second part?: {wrote_script2}")
        else:
            print("There were too few samples to use for the male and female populations")
    else:
        print(f"The following parameters are incorrect: {incorrectparameters}")


if __name__ == "__main__":
    main()
