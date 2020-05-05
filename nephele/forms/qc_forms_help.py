"""
Defines the help text for the QC pipeline forms.
"""

run_cutadapt_desc = """
    If this box is checked we will attempt to remove the adapters from your samples using Cutadapt and the parameters provided.
    If this box is not checked <strong>adapter trimming will NOT be attempted.</strong>
    """

run_qual_trimming_desc = """
    If this box is checked we will attempt to remove low quality sequences from the end of your reads using using Trimmomatic and the parameters provided.
    If this box is not checked <strong>quality trimming will NOT be attempted.</strong>
    """

run_flash2_merge_desc = """
    If this box is checked we will attempt to merge your read pairs using flash2 and the parameters provided.
    If this box is not checked <strong>merging of read pairs will NOT be attempted.</strong>
    """

error_rate_desc = """
    Maximum allowed error rate.
    """

indel_desc = """
    Allow insertions or deletions of bases when matching adapters.
    """

overlap_desc = """
    Require at least this many bases of overlap between read and adapter for an adapter to be found.
    """

wildcard_read_desc = """
    Interpret IUPAC wildcards (e.g., N) in reads.
    """

wildcard_adapter_desc = """
    Interpret IUPAC wildcards (e.g., N) in adapters.
    """

primer_f_desc = """
    Must be a valid IUPAC sequence. Maximum 250bp.
    """

primer_r_desc = """
    Must be a valid IUPAC sequence. Maximum 250bp.
    """

adapter_desc = """
    Sequence of an adapter ligated to the 3' end. The adapter and any subsequent bases
    are trimmed. If a <code>$</code> is appended, the adapter is only found if it is at the end of
    the read. If your sequence of interest is "framed" by a 5' and a 3' adapter, use this
    parameter to define a "linked" primer - see https://cutadapt.readthedocs.io for complete
    details.  Must be a valid IUPAC sequence. Maximum 250bp.
    """

front_desc = """
    Sequence of an adapter ligated to the 5' end. The adapter and any preceding bases are
    trimmed. Partial matches at the 5' end are allowed. If a <code>^</code> character is prepended,
    the adapter is only found if it is at the beginning of the read.  Must be a valid IUPAC sequence. Maximum 250bp.
    """

anywhere_desc = """
    Sequence of an adapter that may be ligated to the 5' or 3' end. Both types of matches
    as described under <code>3' adapter</code> and <code>front 5' adapter</code> are allowed. If the first base of the read is
    part of the match, the behavior is as with <code>front 5'</code>, otherwise as with <code>3' adapter</code>. This
    option is mostly for rescuing failed library preparations - do not use if you know which
    end your adapter was ligated to.  Must be a valid IUPAC sequence. Maximum 250bp.
    """

window_size_desc = """
    Sliding window size for quality trimming, cutting once the average quality within the window
    falls below the <code>required quality</code>. Specifies the number of bases to average across.
    """

req_qual_desc = """
    Specifies the average quality required in the sliding window.
    """

lead_qual_desc = """
    Remove low quality bases from the beginning (5' end) of the read. As long as a base has a
    value below this threshold the base is removed and the next base will be investigated.
    """

trail_qual_desc = """
    Remove low quality bases from the end of the read. As long as a base has a value below this
    threshold the base is removed and the next base (which as Trimmomatic is starting from the 3'
    end would be base preceding the just removed base) will be investigated.
    """

minlen_desc = """
    Remove reads that fall below the specified minimal length.
    """

avg_qual_desc = """
    Drop the read if the average quality across the entire read is below the specified level.
    """

min_overlap_desc = """
    The minimum required overlap length between two reads to provide a confident overlap.
    """

max_overlap_desc = """
    Maximum overlap length expected in approximately 90% of read pairs. Overlaps longer than the
    maximum overlap parameter are still considered as good overlaps, but the mismatch density
    (explained below) is <em>only calculated over the first <code>max_overlap</code> bases in the overlapped region</em>
    rather than the entire overlap.
    """

min_outie_overlap_desc = """
    The minimum required overlap length between two reads to provide a confident overlap in an "outie" scenario.
    """

max_mismatch_density_desc = """
    Maximum allowed ratio between the number of mismatched base pairs and the overlap length.
    """
