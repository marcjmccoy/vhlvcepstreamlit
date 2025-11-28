import re


######################### VHL CONSTANTS #######################################
#TODO eventually need to productionize and make constants OOP call

# Exon boundaries for NM_000551.4 (VHL); adjust if transcript changes
VHL_EXON_BOUNDS = [
    (1, 340),    # Exon 1: c.1 - c.340 [web:22]
    (341, 463),  # Exon 2: c.341 - c.463 [web:24]
    (464, 642)   # Exon 3: c.464 - c.642 [web:21][web:25]
]

# Functional domains for canonical pVHL213 as per GN078:
# Beta (β) domain (AA 63–155, Nuclear Export),
# Alpha (ɑ) domain (AA 156–192, Elongin C binding),
# Second Beta (β) domain (AA 193–204). [web:1][web:2][web:29]
CRITICAL_DOMAINS = [
    (63, 155),   # 1st beta domain
    (156, 192),  # alpha domain
    (193, 204)   # 2nd beta domain
]

C_TERM = (205, 213)  # C-terminal tail of pVHL213

VHL_PROTEIN_LENGTH = 213  # AA length of canonical pVHL213


######################### cDNA-LEVEL HELPERS ##################################

def parse_cdna(hgvs):
    """
    Extract cDNA start/end positions from an HGVS string like:
        'NM_000551.4(VHL):c.191G>C' or '...:c.191_193del'
    Returns:
        (cdna_start, cdna_end) as ints or (None, None) on failure.
    """
    m = re.search(r"c\.(\d+)(?:_(\d+))?", hgvs)
    if m:
        start = int(m.group(1))
        end = int(m.group(2)) if m.group(2) else start
        return start, end
    m_single = re.search(r"c\.(\d+)", hgvs)
    if m_single:
        pos = int(m_single.group(1))
        return pos, pos
    return None, None


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


def in_cterm(codon):
    """
    True if codon lies in the C-terminal tail (AA 205–213).
    """
    return codon is not None and 205 <= codon <= 213


def predict_nmd(codon):
    """
    Very simple VHL-specific NMD heuristic: treat nonsense around codons
    55–136 as likely to trigger NMD (approximate midcoding region). [web:2][web:42]
    """
    return codon is not None and 55 <= codon <= 136


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


######################### PROTEIN-LEVEL PARSING ###############################

def parse_protein_change(hgvs_protein):
    """
    Parse the protein-level HGVS fragment from a full HGVS string like:
        'NM_000551.4(VHL):c.191G>C (p.Arg64Pro)'

    Returns:
        dict with keys:
            'type': one of ['missense', 'inframe_del', 'inframe_ins',
                            'inframe_indel', 'stop_loss', 'other']
            'start_aa': int or None (first affected AA, 1-based)
            'end_aa': int or None (last affected AA, 1-based)
    """
    m = re.search(r"\((p\.[^)]+)\)", hgvs_protein)
    if not m:
        return {"type": "other", "start_aa": None, "end_aa": None}
    p = m.group(1)

    # Simple missense / synonymous / nonsense etc: p.Arg64Pro, p.Arg64=
    missense_like = re.match(r"p\.([A-Z*])(\d+)([A-Z*]|=)", p)
    if missense_like:
        aa_pos = int(missense_like.group(2))
        # Identify stop-loss: original is stop, new is not stop
        if missense_like.group(1) == "*" and missense_like.group(3) != "*":
            return {"type": "stop_loss", "start_aa": aa_pos, "end_aa": aa_pos}
        return {"type": "missense", "start_aa": aa_pos, "end_aa": aa_pos}

    # In-frame deletion: p.Lys171del or p.Lys171_Glu173del
    del_match = re.match(r"p\.([A-Z*])(\d+)(?:_([A-Z*])(\d+))?del", p)
    if del_match:
        start = int(del_match.group(2))
        end = int(del_match.group(4)) if del_match.group(4) else start
        return {"type": "inframe_del", "start_aa": start, "end_aa": end}

    # In-frame insertion: p.Lys171_Glu172insGly or p.Lys171_Glu172insGlyVal
    ins_match = re.match(r"p\.([A-Z*])(\d+)_([A-Z*])(\d+)ins", p)
    if ins_match:
        start = int(ins_match.group(2))
        end = int(ins_match.group(4))
        # A pure insertion is between start and end; treat as affecting that interval
        return {"type": "inframe_ins", "start_aa": start, "end_aa": end}

    # In-frame indel (delins) explicitly marked as in-frame:
    # e.g., p.Lys171_Glu173delinsAsnVal (CDS multiple of 3)
    indel_match = re.match(r"p\.([A-Z*])(\d+)_([A-Z*])(\d+)delins", p)
    if indel_match:
        start = int(indel_match.group(2))
        end = int(indel_match.group(4))
        return {"type": "inframe_indel", "start_aa": start, "end_aa": end}

    # Stop-loss explicit with extension length, e.g. p.*214Leuext*?
    stop_loss_match = re.match(r"p\.\*?(\d+)[A-Z]+ext\*?", p)
    if stop_loss_match:
        start = int(stop_loss_match.group(1))
        return {"type": "stop_loss", "start_aa": start, "end_aa": start}

    # Fallback
    return {"type": "other", "start_aa": None, "end_aa": None}


def affects_beta_alpha_domains(start_aa, end_aa):
    """
    Returns True if the in-frame change overlaps the B/alpha/second-B domains
    (AA 63–204) as defined in GN078. [web:1][web:2]
    """
    if start_aa is None or end_aa is None:
        return False
    dom_start, dom_end = 63, 204
    return not (end_aa < dom_start or start_aa > dom_end)


def entirely_before_54(start_aa, end_aa):
    """
    Returns True if the entire affected AA interval is <54.
    PM4 does not apply to in-frame indels prior to codon 54 that
    do not alter the Met54 VHL p19 codon and beyond. [web:2]
    """
    if start_aa is None or end_aa is None:
        return False
    return end_aa < 54


############################# PM4 CLASSIFIER ##################################

def classify_vhl_pm4(hgvs_str):
    """
    Implement VHL-specific PM4 per GN078:

      - Moderate strength when:
          * In-frame insertion/deletion/indel in the B and alpha domains
            (AA 63–204; beta, alpha, second beta domains), OR
          * Stop-loss variant that adds significant additional amino acids
            beyond the canonical stop codon of VHL. [web:2][web:34]

      - PM4 does NOT apply to in-frame indels prior to codon 54 that do not
        alter the Met54 VHL p19 codon and beyond. [web:2]

    Returns:
        dict: {"strength": <str or None>, "context": <str>}
    """

    # Extract and parse protein-level HGVS
    prot_info = parse_protein_change(hgvs_str)
    vt = prot_info["type"]
    start_aa = prot_info["start_aa"]
    end_aa = prot_info["end_aa"]

    # If cannot parse, do not apply PM4
    if vt == "other" or start_aa is None:
        return {
            "strength": None,
            "context": "Unable to confidently parse protein-level change; PM4 not applied for VHL."
        }

    # Exclude events fully prior to codon 54 and not touching/altering Met54+
    if entirely_before_54(start_aa, end_aa):
        return {
            "strength": None,
            "context": "In-frame change entirely prior to codon 54 and not affecting Met54 or beyond; VHL PM4 not applied."
        }

    # In-frame indels in B / alpha / second-B domains (AA 63–204)
    if vt in ["inframe_del", "inframe_ins", "inframe_indel"]:
        if affects_beta_alpha_domains(start_aa, end_aa):
            return {
                "strength": "PM4",
                "context": "In-frame insertion/deletion within VHL beta/alpha domains (AA 63–204); assign PM4 (Moderate)."
            }
        else:
            return {
                "strength": None,
                "context": "In-frame insertion/deletion outside VHL beta/alpha domains (AA 63–204); VHL PM4 not applied."
            }

    # Stop-loss variants adding significant additional amino acids
    if vt == "stop_loss":
        # GN078 emphasizes “significant” C-terminal extension, citing multiple
        # pathogenic cases and functional data, but does not specify an exact
        # numeric cutoff. [web:2]
        # Here, treat any well-defined extension beyond the canonical stop
        # that creates a new tail (start_aa > VHL_PROTEIN_LENGTH) as qualifying
        # for PM4, leaving finer downgrading to other criteria if needed.
        if start_aa > VHL_PROTEIN_LENGTH:
            return {
                "strength": "PM4",
                "context": "Stop-loss variant adding a C-terminal extension beyond VHL codon 213; assign PM4 (Moderate)."
            }
        else:
            return {
                "strength": None,
                "context": "Stop-loss extension cannot be confidently interpreted as adding a significant C-terminal tail; VHL PM4 not applied."
            }

    # Missense, synonymous, frameshift, early nonsense, etc. – PM4 not relevant
    return {
        "strength": None,
        "context": "Variant type is not an in-frame indel in AA 63–204 or a qualifying stop-loss extension; VHL PM4 not applied."
    }
