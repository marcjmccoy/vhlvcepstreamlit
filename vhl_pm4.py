import re


######################### VHL CONSTANTS #######################################

# Exon boundaries for NM_000551.4 (VHL); adjust if transcript changes
VHL_EXON_BOUNDS = [
    (1, 340),    # Exon 1: c.1 - c.340
    (341, 463),  # Exon 2: c.341 - c.463
    (464, 642)   # Exon 3: c.464 - c.642
]

# Functional domains for canonical pVHL213 as per GN078:
# Beta (β) domain (AA 63–155), Alpha (ɑ) domain (AA 156–192),
# Second Beta (β) domain (AA 193–204)
CRITICAL_DOMAINS = [
    (63, 155),   # 1st beta domain
    (156, 192),  # alpha domain
    (193, 204)   # 2nd beta domain
]

C_TERM = (205, 213)          # C-terminal tail of pVHL213
VHL_PROTEIN_LENGTH = 213     # AA length of canonical pVHL213


######################### cDNA-LEVEL HELPERS ##################################


def parse_cdna_interval_and_op(hgvs):
    """
    Parse cDNA interval and operation from a c.HGVS string, e.g.:
        NM_000551.4(VHL):c.189_191del
        NM_000551.4(VHL):c.190_195dup
        NM_000551.4(VHL):c.190_195delinsTGC
        NM_000551.4(VHL):c.191G>C

    Returns:
        dict with keys:
            'start': int or None
            'end': int or None
            'op': one of ['del', 'dup', 'ins', 'delins', 'sub', None]
            'op_seq': optional inserted sequence string (for dup/ins/delins)
    """
    # Generic cDNA pattern capturing interval and op
    m = re.search(
        r"c\.(\d+)(?:_(\d+))?"
        r"(?:(delins|del|dup|ins)|([ACGT]>[ACGT]))?"
        r"([ACGT]*)?",
        hgvs
    )
    if not m:
        return {"start": None, "end": None, "op": None, "op_seq": None}

    start = int(m.group(1))
    end = int(m.group(2)) if m.group(2) else start

    op = None
    op_seq = None

    if m.group(3):  # delins / del / dup / ins
        op = m.group(3)
        op_seq = m.group(5) if m.group(5) else None
    elif m.group(4):  # substitution, like 191G>C
        op = "sub"

    return {"start": start, "end": end, "op": op, "op_seq": op_seq}


def cdna_to_codon(cdna_pos):
    """
    Convert cDNA position (1-based, coding only) to codon index (1-based).
    """
    return ((cdna_pos - 1) // 3) + 1 if cdna_pos else None


def in_critical_domain(codon):
    """
    True if codon lies within any of the VHL critical domains (AA 63–204).
    """
    return any(start <= codon <= end for start, end in CRITICAL_DOMAINS)


def entirely_before_54(start_aa, end_aa):
    """
    Returns True if the entire affected AA interval is <54.
    PM4 does not apply to in-frame indels prior to codon 54 that
    do not alter the Met54 VHL p19 codon and beyond.
    """
    if start_aa is None or end_aa is None:
        return False
    return end_aa < 54


def is_exon_boundary(cdna_start, cdna_end):
    """
    Returns True if the interval completely covers any known exon
    (i.e., the full exon is within the deleted/affected region).
    """
    if cdna_start is None or cdna_end is None:
        return False
    for exon_start, exon_end in VHL_EXON_BOUNDS:
        if cdna_start <= exon_start and cdna_end >= exon_end:
            return True
    return False


######################### cDNA → PROTEIN INFERENCE ############################


def infer_inframe_type_from_cdna(cdna_info):
    """
    Given parsed cDNA info (start, end, op, op_seq), infer whether the variant
    is an in-frame deletion/insertion/indel affecting coding sequence.

    Returns:
        ('inframe_del' | 'inframe_ins' | 'inframe_indel' | 'other')
    """
    start = cdna_info["start"]
    end = cdna_info["end"]
    op = cdna_info["op"]
    op_seq = cdna_info["op_seq"]

    if start is None or end is None or op is None:
        return "other"

    length = end - start + 1

    # Deletions
    if op == "del":
        return "inframe_del" if length % 3 == 0 else "other"

    # Duplications: assume length equals duplicated inserted length
    if op == "dup":
        return "inframe_ins" if length % 3 == 0 else "other"

    # Pure insertions (between nucleotides)
    if op == "ins":
        ins_len = len(op_seq) if op_seq else 0
        return "inframe_ins" if ins_len % 3 == 0 and ins_len > 0 else "other"

    # Delins: need net effect on CDS length to be multiple of 3
    if op == "delins":
        del_len = length
        ins_len = len(op_seq) if op_seq else 0
        net = ins_len - del_len
        return "inframe_indel" if net % 3 == 0 else "other"

    # Substitutions are not PM4-relevant here
    return "other"


def cdna_interval_to_codon_interval(cdna_start, cdna_end):
    """
    Map a cDNA interval to the corresponding codon interval.
    For multi-nucleotide events, this gives the minimal codon span.
    """
    if cdna_start is None or cdna_end is None:
        return None, None
    aa_start = cdna_to_codon(cdna_start)
    aa_end = cdna_to_codon(cdna_end)
    return aa_start, aa_end


def affects_beta_alpha_domains(start_aa, end_aa):
    """
    Returns True if the in-frame change overlaps the B/alpha/second-B domains
    (AA 63–204) as defined in GN078.
    """
    if start_aa is None or end_aa is None:
        return False
    dom_start, dom_end = 63, 204
    return not (end_aa < dom_start or start_aa > dom_end)


############################# PM4 CLASSIFIER ##################################


def classify_vhl_pm4(hgvs_str):
    """
    Implement VHL-specific PM4 per GN078 using cDNA-only HGVS.

      - Moderate strength when:
          * In-frame insertion/deletion/indel in the B and alpha domains
            (AA 63–204; beta, alpha, second beta domains), OR
          * Stop-loss variant that adds significant additional amino acids
            beyond the canonical stop codon of VHL.

      - PM4 does NOT apply to in-frame indels prior to codon 54 that do not
        alter the Met54 VHL p19 codon and beyond.

    Expects an input like:
        'NM_000551.4(VHL):c.189_191del'
        'NM_000551.4(VHL):c.190_195dup'
    """

    # 1) Parse cDNA and operation
    cdna_info = parse_cdna_interval_and_op(hgvs_str)
    cdna_start = cdna_info["start"]
    cdna_end = cdna_info["end"]

    if cdna_start is None or cdna_end is None:
        return {
            "strength": None,
            "context": "Unable to parse cDNA interval; PM4 not applied for VHL."
        }

    # 2) Infer in-frame type from cDNA
    vt = infer_inframe_type_from_cdna(cdna_info)

    # Handle stop-loss separately if you later extend this to support p. notation.
    if vt == "other":
        return {
            "strength": None,
            "context": "Variant is not an in-frame indel at the cDNA level; VHL PM4 not applied."
        }

    # 3) Convert cDNA interval to codon interval
    aa_start, aa_end = cdna_interval_to_codon_interval(cdna_start, cdna_end)

    if aa_start is None or aa_end is None:
        return {
            "strength": None,
            "context": "Unable to map cDNA interval to codons; PM4 not applied for VHL."
        }

    # 4) Exclude events fully prior to codon 54
    if entirely_before_54(aa_start, aa_end):
        return {
            "strength": None,
            "context": "In-frame change entirely prior to codon 54 and not affecting Met54 or beyond; VHL PM4 not applied."
        }

    # 5) In-frame indels in B / alpha / second-B domains (AA 63–204)
    if vt in ["inframe_del", "inframe_ins", "inframe_indel"]:
        if affects_beta_alpha_domains(aa_start, aa_end):
            return {
                "strength": "PM4",
                "context": "In-frame insertion/deletion within VHL beta/alpha domains (AA 63–204); assign PM4 (Moderate)."
            }
        else:
            return {
                "strength": None,
                "context": "In-frame insertion/deletion outside VHL beta/alpha domains (AA 63–204); VHL PM4 not applied."
            }

    # Fallback
    return {
        "strength": None,
        "context": "Variant type is not an in-frame indel in AA 63–204; VHL PM4 not applied."
    }
