import os
from itertools import combinations
from collections import defaultdict
import numpy as np
from scipy.optimize import nnls
from scipy.stats.stats import pearsonr
from sklearn import linear_model
import re
from ice.classes.edit_proposal_creator import EditProposalCreator
from ice.classes.guide_target import GuideTarget
from ice.classes.ice_result import ICEResult
from ice.classes.pair_alignment import PairAlignment
from ice.classes.proposal_base import ProposalBase
from ice.classes.sanger_object import SangerObject
from ice.outputs.create_discordance_indel_files import generate_discordance_indel_files
from ice.outputs.create_json import write_individual_contribs, write_contribs_json, write_all_proposals_json
from ice.outputs.create_trace_files import generate_trace_files
from ice.utility.sequence import RNA2DNA, reverse_complement


class SangerAnalysis:
    """
    Runs Sanger sequencing analysis to infer CRISPR edit outcomes.
    """

    iced_correction_factor = 1.41  
    HDR_OVERLAP_FILTER_CUTOFF = 3  

    def __init__(self, verbose=False):
        self.control_sample = None  
        self.edited_sample = None  
        self.gRNA_sequences = []
        self.guide_targets = []
        self.valid_configuration = False
        self.donor_odn = None
        self.recombination_changed_bases = []
        self.proposals = None
        self.donor_alignment = None
        self.alignment_window = None  
        self.inference_window = None  
        self.alignment_sequence = None
        self.indel_max_size = None
        self.base_outputname = None
        self.results = ICEResult()
        self.r_squared_correction = True
        self.verbose = verbose
        self.debug = False
        self.warnings = []
        self.allprops = False  

    def initialize_with(self, control_path, edited_path, gRNA_sequences,
                        alignment_window=None, inference_window=None,
                        indel_max_size=5, base_outputname=None,
                        donor=None, allprops=False):

        if os.path.exists(control_path):
            self.control_sample = SangerObject()
            self.control_sample.initialize_from_path(control_path)
        else:
            raise Exception('Control path does not exist')

        if os.path.exists(edited_path):
            self.edited_sample = SangerObject()
            self.edited_sample.initialize_from_path(edited_path)
        else:
            raise Exception('Edited path does not exist')

        if gRNA_sequences is not None and isinstance(gRNA_sequences, str):
            self.gRNA_sequences = [RNA2DNA(seq.strip().upper()) for seq in gRNA_sequences.split(",")]

        if base_outputname:
            self.base_outputname = base_outputname + "."
        self.inference_window = inference_window
        self.indel_max_size = indel_max_size
        self.donor_odn = donor
        self.allprops = allprops

        return True

    def analyze_sample(self):
        self.find_targets()

        if self.donor_odn:
            self.check_recombination()

        self.find_alignment_window()
        alignment = PairAlignment(self.control_sample.primary_base_calls, self.edited_sample.primary_base_calls)
        self.alignment = alignment

        alignment.align_all()
        alignment.align_with_window(self.alignment_window)

        if not alignment.has_alignment:
            raise Exception("No alignment found between control and edited sample")

        aln_file = self.base_outputname + "all.txt"
        alignment.write_aln(alignment.all_aligned_clustal, to_file=aln_file)

        self._generate_edit_proposals()
        self._calculate_inference_window()
        self._generate_coefficient_matrix()
        self._generate_outcomes_vector()
        self.analyze_and_rank()
        self.calculate_discordance()
        self.simple_discordance_algorithm()

        indel_file = self.base_outputname + "indel.json"
        generate_discordance_indel_files(self, self.results, to_file=indel_file)
        write_contribs_json(self, self.base_outputname + "contribs.json")

