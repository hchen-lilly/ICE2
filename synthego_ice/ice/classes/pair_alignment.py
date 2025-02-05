from io import StringIO
import json
import io
from Bio.Align import PairwiseAligner
from Bio import AlignIO


class PairAlignment:
    """
    Uses a list of lists and two dictionaries to store data about a pairwise alignment.
    """

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.name1 = "control"
        self.name2 = "edited"
        self.all_aligned_seqs = None
        self.all_aligned_clustal = None
        self.aln_clustal = None

    @staticmethod
    def write_aln(aln_clustal, to_file=None):
        if to_file:
            with open(to_file, 'w') as out_file:
                print(aln_clustal, file=out_file)

    @staticmethod
    def write_json(aln_seqs, to_file):
        if aln_seqs:
            with open(to_file, 'w') as f:
                json.dump({'control': aln_seqs[0], 'edited': aln_seqs[1]}, f)

##    def align_list_to_clustal(self, aln, name1, name2):
##        my_aln = f">{name1}\n{aln[0]}\n>{name2}\n{aln[1]}"
##        f = StringIO(my_aln)
##        aln_objs = list(AlignIO.parse(f, "fasta"))
##
##        output_handle = io.StringIO()
##        AlignIO.write(aln_objs, output_handle, "clustal")
##        return output_handle.getvalue()
    def align_list_to_clustal(self, aln, name1, name2):
        my_aln = f">{name1}\n{aln[0]}\n>{name2}\n{aln[1]}"
        f = StringIO(my_aln)
        
        # Read the alignment as a MultipleSeqAlignment object
        aln_objs = list(AlignIO.read(f, "fasta"))

        # Ensure there's at least one alignment
        if not aln_objs:
            return "No valid alignment found."

        # Write the alignment in Clustal format
        output_handle = io.StringIO()
        AlignIO.write([aln_objs], output_handle, "clustal")
        return output_handle.getvalue()

    def align_all(self):
        aligner = PairwiseAligner()
        aligner.mode = "local"
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -3
        aligner.extend_gap_score = -1

        alignments = list(aligner.align(self.seq1, self.seq2))

        if not alignments:
            raise ValueError("No valid alignment found")

        aln = alignments[0]
        self.all_aligned_seqs = (str(aln[0]), str(aln[1]))
        self.all_aligned_clustal = self.align_list_to_clustal(self.all_aligned_seqs, "control", "edited")

