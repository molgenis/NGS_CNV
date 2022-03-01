#!/usr/bin/env python
import argparse


def get_params():
    etop_args = argparse.ArgumentParser()
    etop_args.add_argument("-b", "--bed-to-bam", type=str, required=True, dest="bed-to-bam", help="Path to BED to BAM file")
    etop_args.add_argument("-e", "--experiment-data", type=str, required=True, dest="experiment-data", help="Path to RD3 Experiment data")
    etop_args.add_argument("-o", "--outfile", type=str, required=True, dest="outfile", help="Path to write output file to")
    etop_args.add_argument("-s", "--samples-data", type=str, required=True, dest="samples-data", help="Path to RD3 samples table")
    etop_args.add_argument("-u", "--subjects-data", type=str, required=True, dest="subjects-data", help="Path to RD3 subjects table")
    return vars(etop_args.parse_args())


def determine_field_indices(header_line_data, fields_to_select):
    """Determine and return the location of the required fields in the file.

    Parameters
    ----------
    header_line_data : list of str
        Column names of the header line
    fields_to_select : list of str
        

    Returns
    -------
    field_indices : list of int
        Column indices of the required columns
    """
    field_indices = {}
    print(header_line_data)
    for fts in fields_to_select:
        if fts in header_line_data:
            field_indices[fts] = header_line_data.index(fts)
    return field_indices


def read_experiment_data(experiment_fileloc):
    """Read the RD3 experiment table.

    Parameters
    ----------
    experiment_fileloc : str
        Path to RD3 experiment file

    Returns
    -------
    rd3_ex_data : dict
        
    """
    header_fields = ["experimentID", "sample"]
    rd3_ex_data = {}

    try:
        with open(experiment_fileloc, 'r') as experimentfile:
            headerline = next(experimentfile)
            data_indices = determine_field_indices(headerline.strip().split("\t"), header_fields)

            for fileline in experimentfile:
                filelinedata = fileline.strip().split("\t")
                expid = filelinedata[data_indices["experimentID"]]
                sampleid = filelinedata[data_indices["sample"]]
                rd3_ex_data[sampleid] = expid
    except IOError:
        print("Could not read RD3 experiment file...")
    finally:
        return rd3_ex_data


def read_samples_data(samples_fileloc):
    header_fields = ["sampleID", "subject"]
    rd3_sam_data = {}

    try:
        with open(samples_fileloc, 'r') as samplesfile:
            headerline = next(samplesfile)
            data_indices = determine_field_indices(headerline.strip().split("\t"), header_fields)
            maxindex = get_max_col_index(data_indices)

            for fileline in samplesfile:
                filelinedata = fileline.strip().split("\t")
                if len(filelinedata) > maxindex:
                    sampleid = filelinedata[data_indices["sampleID"]]
                    subjectid = filelinedata[data_indices["subject"]]
                    if sampleid != "" and subjectid != "":
                        rd3_sam_data[sampleid] = subjectid
    except IOError:
        print("Could not read RD3 samples file...")
    finally:
        return rd3_sam_data


def read_subjects_data(subjects_fileloc):
    header_fields = ["SubjectID", "Claimed sex"]
    rd3_sub_data = {}

    try:
        with open(subjects_fileloc, 'r') as subjectsfile:
            headerline = next(subjectsfile)
            data_indices = determine_field_indices(headerline.strip().split("\t"), header_fields)
            maxindex = get_max_col_index(data_indices)

            for fileline in subjectsfile:
                filelinedata = fileline.strip().split("\t")
                if len(filelinedata) > maxindex:
                    subjectid = filelinedata[data_indices["SubjectID"]]
                    claimedsex = filelinedata[data_indices["Claimed sex"]]
                    if subjectid != "" and claimedsex != "":
                        rd3_sub_data[subjectid] = claimedsex
    except IOError:
        print("Could not read RD3 subjects file...")
    finally:
        return rd3_sub_data


def get_bam_paths(bedtobam_fileloc):
    bampath_data = {}
    try:
        with open(bedtobam_fileloc, 'r') as btbfile:
            next(btbfile)
            for fileline in btbfile:
                filelinedata = fileline.strip().split("\t")
                samplename = filelinedata[1].split("/")[-1].split(".")[0]
                if samplename not in bampath_data:
                    bampath_data[samplename] = filelinedata[1]
    except IOError:
        print("Could not read BED to BAM file...")
    finally:
        return bampath_data


def write_sex_to_bam(experimentdata, samplesdata, subjectsdata, bampathdata, outfileloc):
    file_written = False
    try:
        with open(outfileloc, 'w') as outfile:
            outfile.write("ESample\tBAM file\tSample sex\n")
            for vsename in experimentdata:
                bamname = experimentdata[vsename]
                pnumber = get_data_element(vsename, samplesdata)
                samplesex = get_data_element(pnumber, subjectsdata)
                bamfilepath = get_data_element(bamname, bampathdata)
                if bamfilepath == "":
                    bamfilepath = "NA"
                outfile.write(f"{bamname}\t{bamfilepath}\t{samplesex}\n")
        file_written = True
    except IOError:
        print("Could not write sex to BAM file")
    finally:
        return file_written


def get_data_element(keytouse, datatouse):
    datavalue = ""
    if keytouse in datatouse:
        datavalue = datatouse[keytouse]
    return datavalue


def get_max_col_index(colindices):
    """Return the maximum column index.

    Parameters
    ----------
    colindices : dict
        Column names and their indices

    Returns
    -------
    maxcolindex : int
        Last column index
    """
    maxcolindex = 0
    for colname in colindices:
        if colindices[colname] > maxcolindex:
            maxcolindex = colindices[colname]
    return maxcolindex


def main():
    etop_params = get_params()

    print("...[READING RD3 EXPERIMENT DATA]...")
    experiment_data = read_experiment_data(etop_params["experiment-data"])

    print("...[READING RD3 SAMPLES DATA]...")
    samples_data = read_samples_data(etop_params["samples-data"])

    print("...[READING RD3 SUBJECTS DATA]...")
    subjects_data = read_subjects_data(etop_params["subjects-data"])

    print("...[READING PATHGS TO BAM FILES]...")
    bampath_data = get_bam_paths(etop_params["bed-to-bam"])

    print("...[WRITING OUTPUT FILE LINKING SEX TO BAM FILE]...")
    wrote_file = write_sex_to_bam(experiment_data, samples_data, subjects_data, bampath_data, etop_params["outfile"])
    print(f"...[WROTE OUTPUT FILE SUCCESFULLY?: {wrote_file}]...")


if __name__ == "__main__":
    main()
