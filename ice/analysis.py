import argparse
import os
import traceback

import pandas as pd
import datetime
import json

from ice.classes.ice_result import ICEResult
from ice.classes.sanger_analysis import SangerAnalysis
from ice.utility.misc import version
from ice.utility.sequence import is_nuc_acid

from .__version__ import __version__


def single_sanger_analysis(control_path, sample_path, base_outputname, guide, donor=None, verbose=False,
                           allprops=False):
    if control_path is None or not os.path.exists(control_path):
        raise Exception('Control @ {} not found'.format(control_path))

    if not os.path.exists(sample_path):
        raise Exception('Experiment sample @ {} not found'.format(sample_path))

    base_dir = os.path.join(*os.path.split(os.path.abspath(base_outputname))[:-1])
    if verbose:
        print('Base dir: %s' % base_dir)

    if not os.path.exists(base_dir):
        os.makedirs(base_dir, exist_ok=True)

    sa = SangerAnalysis(verbose=verbose)
    sa.debug = False
    sa.allprops = allprops
    sa.initialize_with(control_path=control_path,
                       edited_path=sample_path,
                       gRNA_sequences=guide,
                       indel_max_size=20,
                       base_outputname=base_outputname,
                       donor=donor,
                       allprops=allprops)
    try:
        sa.analyze_sample()
        return sa.results.to_json(sa.guide_targets, sa.warnings)

    except Exception as e:
        results = ICEResult()
        print('Exception Caught!')
        traceback.print_exc()
        return results.to_json(sa.guide_targets, [str(e)])


def single_sanger_analysis_cli():
    parser = argparse.ArgumentParser(description='Analyze Sanger reads to Infer Crispr Edit outcomes')
    parser.add_argument('--control', dest='control', required=True, help='The wildtype / unedited ab1 file (REQUIRED)')
    parser.add_argument('--edited', dest='edited', required=True, help='The edited ab1 file (REQUIRED)')
    parser.add_argument('--target', dest='target', required=True, help='Target sequence(s) (17-23 bases, RNA or DNA)')
    parser.add_argument('--out', dest='out', default=None, help='Output base path (Defaults to ./results/single)')
    parser.add_argument('--donor', dest='donor', default=None, help='Donor DNA sequence for HDR (Optional)')
    parser.add_argument('--verbose', dest='verbose', action='store_true')
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('--allprops', dest='allprops', action='store_true', default=False,
                        help="Output all Edit Proposals, even if they have zero contribution")

    args = parser.parse_args()

    assert os.path.isfile(args.control)
    assert os.path.isfile(args.edited)

    out_dir = os.path.abspath(args.out) if args.out else os.path.join(os.path.abspath('.'), 'results', 'single')
    os.makedirs(os.path.dirname(out_dir), exist_ok=True)

    print('Synthego ICE (https://synthego.com)')
    print('Version: {}'.format(__version__))

    single_sanger_analysis(control_path=args.control,
                           sample_path=args.edited,
                           base_outputname=out_dir,
                           guide=args.target,
                           donor=args.donor,
                           verbose=args.verbose,
                           allprops=args.allprops)


def multiple_sanger_analysis(definition_file, output_dir, data_dir=None, verbose=False, single_line=None, allprops=False):
    input_df = pd.read_excel(definition_file)
    input_df = input_df.rename(columns={"Donor Sequence": "Donor", "Control": "Control File", "Experiment": "Experiment File"})

    results = []
    fails = []
    jobs = []
    n = 0

    for m, experiment in input_df.iterrows():
        label = experiment['Label']
        base_outputname = os.path.join(output_dir, f'{n}-{label}')

        control_sequence_file = experiment['Control File']
        edit_sequence_file = experiment['Experiment File']
        guide = experiment['Guide Sequence']
        donor = experiment['Donor'] if 'Donor' in experiment and is_nuc_acid(experiment['Donor']) else None

        try:
            if pd.isnull(control_sequence_file):
                raise IOError(f"Control filepath not specified at line {n+1} in definition file")
            if pd.isnull(edit_sequence_file):
                raise IOError(f"Edit filepath not specified at line {n+1} in definition file")

            control_sequence_path = os.path.join(data_dir, control_sequence_file)
            edit_sequence_path = os.path.join(data_dir, edit_sequence_file)

            if single_line is not None and n != single_line:
                continue

            print("-" * 50, "analyzing", n, experiment['Label'], guide)

            job_args = (control_sequence_path, edit_sequence_path, base_outputname, guide)
            job_kwargs = {'verbose': verbose, 'allprops': allprops, 'donor': donor}
            result = single_sanger_analysis(*job_args, **job_kwargs)
            jobs.append((experiment, result))

        except Exception as e:
            fails.append(experiment)
            print("Single Sanger analysis failed", e)
            traceback.print_exc()

        n += 1

    for job in jobs:
        r = job[1]
        experiment = job[0]
        donor = experiment['Donor'] if 'Donor' in experiment and is_nuc_acid(experiment['Donor']) else None

        if r is not None:
            tmp = [experiment['Label'], r['ice'], r['ice_d'], r['rsq'], r['hdr_pct'], r['ko_score'], r['guides'],
                   r['notes'], experiment['Experiment File'], experiment['Control File'], donor]
        else:
            tmp = [experiment['Label'], 'Failed', '', '', '', '', '', '', '', '']
        results.append(tmp)

    if results:
        input_df = pd.DataFrame(results)
        timestamp = datetime.datetime.now().strftime('%Y-%m-%d-%H%M%S')
        out_file = os.path.join(output_dir, f"ice.results.{timestamp}.xlsx")

        header = ["sample_name", "ice", 'ice_d', "r_squared", "hdr_pct", "ko_score", "guides", "notes",
                  "experiment_file", "control_file", "donor"]
        input_df.columns = header

        out_dict = [dict(zip(header, r)) for r in results]
        with open(out_file.replace('.xlsx', '.json'), 'w') as f:
            json.dump(out_dict, f, ensure_ascii=False)

        # ✅ FIX: Removed writer.close()
        with pd.ExcelWriter(out_file, engine='openpyxl') as writer:
            input_df.to_excel(writer, sheet_name="Results")
            metadata = pd.DataFrame.from_dict([{'version': __version__}])
            metadata.to_excel(writer, sheet_name='Metadata')

        return out_dict
    else:
        print("None of the samples were able to be analyzed")
        return False


def multiple_sanger_analysis_cli():
    parser = argparse.ArgumentParser(description='Analyze Sanger reads to infer crispr edit outcomes')
    parser.add_argument('--in', dest='input', required=True, help='Input definition file in Excel xlsx format')
    parser.add_argument('--out', dest='out', default=None, help='Output directory path (defaults to .)')
    parser.add_argument('--data', dest='data', required=True, help='Data path, where .ab1 files are located')
    parser.add_argument('--verbose', dest='verbose', action='store_true', help='Display verbose output')
    parser.add_argument('--line', dest='line', default=None, type=int,
                        help="Only run specified line in the Excel xlsx definition file")
    parser.add_argument('--allprops', dest='allprops', action='store_true', default=False,
                        help="Output all Edit Proposals, even if they have zero contribution")
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))

    args = parser.parse_args()
    out_dir = os.path.abspath(args.out) if args.out else os.path.abspath(os.path.curdir)
    data_dir = os.path.abspath(args.data)
    os.makedirs(out_dir, exist_ok=True)

    print('Synthego ICE (https://synthego.com)')
    print('Version: {}'.format(__version__))

    multiple_sanger_analysis(args.input, output_dir=out_dir, data_dir=data_dir,
                             verbose=args.verbose, single_line=args.line, allprops=args.allprops)
