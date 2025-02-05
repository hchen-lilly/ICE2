from io import StringIO
import json
import io

from Bio.Align import PairwiseAligner
from Bio import AlignIO


class PairAlignment:
    '''
    Uses a list of lists and two dictionaries to store data about a pairwise alignment
    '''

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.name1 = "control"
        self.name2 = "edited"
        self.alignment_pairs = None
        self.all_aligned_seqs = None
        self.aln_seqs = None

        self.sample_to_control = {}
        self.control_to_sample = {}
        self.all_aligned_clustal = None
        self.aligned_clustal = None
        self.last_aligned_pair_idx = None

        self.aln_clustal = None

    def __str__(self):
        return self.all_aligned_seqs[0] + "\n" + self.all_aligned_seqs[1]

    @staticmethod
    def write_aln(aln_clustal, to_file=None):
        if to_file is not None:
            with open(to_file, 'w') as out_file:
                print(aln_clustal, file=out_file)

    @staticmethod
    def write_json(aln_seqs, to_file):
        if aln_seqs:
            out_dict = {'control': aln_seqs[0], 'edited': aln_seqs[1]}
            with open(to_file, 'w') as f:
                json.dump(out_dict, f)

    def align_list_to_clustal(self, aln, name1, name2):
        my_aln = ">{}\n".format(name1) + aln[0] + "\n" + ">{}\n".format(name2) + aln[1]
        f = StringIO(my_aln)
        aln_objs = list(AlignIO.parse(f, "fasta"))
        
        output_handle = io.StringIO()
        AlignIO.write(aln_objs, output_handle, "clustal")
        alignment_txt = output_handle.getvalue()
        
        return alignment_txt

    def align_all(self):
        seq1 = self.seq1
        seq2 = self.seq2

        match_bonus = 2
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = match_bonus
        aligner.mismatch_score = -1
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -1

        alignments = aligner.align(seq1, seq2)
        aln = alignments[0]
        
        self.all_aligned_seqs = (str(aln[0]), str(aln[1]))
        self.all_aligned_clustal = self.align_list_to_clustal((str(aln[0]), str(aln[1])), "control", "edited")

    @property
    def has_alignment(self):
        return not not self.alignment_pairs

    def align_with_window(self, alignment_window):
        seq1 = self.seq1
        seq2 = self.seq2
        aw = alignment_window
        window_size = aw[1] - aw[0]

        match_bonus = 2
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = match_bonus
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -1

        alignments = aligner.align(seq1[aw[0]:aw[1]], seq2)
        
        try:
            aln = alignments[0]
        except Exception as e:
            self.aln_clustal = "No alignment found\nCtrl {} to {}:{}\nEdited:{}\n".format(
                aw[0], aw[1], seq1, seq2)
            return False, ('No alignment found ' + str(e))

        aln_score = aln.score
        aln_score_normalized = (aln_score / (window_size * match_bonus)) * 100
        self.aln_clustal = self.align_list_to_clustal((str(aln[0]), str(aln[1])), "control_aln_window", "edited")

        if aln_score_normalized < 50:
            return False, "Poor alignment upstream of cutsite: %0.2f percent of full score" % (aln_score_normalized)

        self.aln_seqs = (str(aln[0]), str(aln[1]))
        return True, "Alignment succeeded"


class DonorAlignment(PairAlignment):
    """
    Adds functionality that is useful in the specific case of pair alignment between control sequence and donor sequence
    """
    MATCH_BONUS = 2

    def __init__(self, control_seq, donor_seq):
        super(DonorAlignment, self).__init__(control_seq, donor_seq)
        self.align_ssodn()
        self.hdr_indel_size = self._calc_hdr_indel_size()

    @property
    def aligned_control_seq(self):
        return self.all_aligned_seqs[0]

    @property
    def aligned_donor_seq(self):
        return self.all_aligned_seqs[1]

    def _calc_hdr_indel_size(self):
        insert_total_len = self.aligned_control_seq.strip('-').count('-')
        deletion_total_len = self.aligned_donor_seq.strip('-').count('-')
        return insert_total_len - deletion_total_len

    @property
    def control_seq(self):
        return self.seq1

    @property
    def donor_seq(self):
        return self.seq2

    def align_ssodn(self):
        print('starting to align ssodn')
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = self.MATCH_BONUS
        aligner.mismatch_score = -1
        aligner.open_gap_score = -6
        aligner.extend_gap_score = -0.5

        alignments = aligner.align(self.control_seq, self.donor_seq)
        alignment = alignments[0]

        self.all_aligned_seqs = (str(alignment[0]), str(alignment[1]))
        self.all_aligned_clustal = self.align_list_to_clustal((str(alignment[0]), str(alignment[1])), 'control', 'donor')

