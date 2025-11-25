import re

# Exon boundaries for NM_000551.4 (VHL); adjust as needed
VHL_EXON_BOUNDS = [
    (1, 118),    # Exon 1: c.1 - c.118
    (119, 340),  # Exon 2: c.119 - c.340
    (341, 639)   # Exon 3: c.341 - c.639
]

CRITICAL_DOMAINS = [
    (63, 154),    # 1st Beta domain – Nuclear Export 114-155
    (155, 192),   # Alpha domain – Elongin C binding 157-172
    (193, 204)    # 2nd Beta domain
]
C_TERM = (205, 213)

def parse_cdna(hgvs):
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
    return ((cdna_pos - 1) // 3) + 1 if cdna_pos else None

def in_critical_domain(codon):
    return any(start <= codon <= end for start, end in CRITICAL_DOMAINS)

def in_cterm(codon):
    return codon is not None and 205 <= codon <= 213

def predict_nmd(codon):
    return codon is not None and 55 <= codon <= 136

def is_exon_boundary(cdna_start, cdna_end):
    """Returns True if the interval OR region completely covers any known exon boundary (full exon is deleted)."""
    for exon_start, exon_end in VHL_EXON_BOUNDS:
        if cdna_start <= exon_start and cdna_end >= exon_end:
            return True
    return False


############################# PVS1 ###################################################
def classify_vhl_pvs1(
    hgvs_str,
    exon_skipping=False,
    cryptic_disrupts_rf=None,
    cryptic_preserves_rf=None,
    duplication_type=None, # "tandem", "not_in_tandem", None
    initiation_codon=None
):
    # Type parsing and position
    cdna_start, cdna_end = parse_cdna(hgvs_str)
    codon_start = cdna_to_codon(cdna_start) if cdna_start else None
    codon_end   = cdna_to_codon(cdna_end) if cdna_end else codon_start
    min_codon = min(codon_start, codon_end) if codon_start and codon_end else None

    critical_domain = in_critical_domain(min_codon) if min_codon else False
    cterm = in_cterm(min_codon)
    nmd = predict_nmd(min_codon)

    vt = None
    canonical = re.search(r"c\.(\d+)(\+1|\+2|\-1|\-2)", hgvs_str)
    cryptic   = re.search(r"c\.(\d+)(\+\d+|\-\d+)", hgvs_str)

    # Canonical GT-AG splice or exon-skipping assigned Very Strong
    if canonical or exon_skipping:
        return {"strength": "PVS1", "context": "Canonical GT-AG +1/-1 splice site or exon skipping in VHL—assign PVS1 (Very Strong)."}
    
    # Cryptic splice logic
    if cryptic:
        pos_num = int(cryptic.group(2)[1:])
        if abs(pos_num) > 2:
            vt = "cryptic_splice"

    # Explicit exon deletion (whole-exon, not detected as small indel)
    elif re.search(r"(del)", hgvs_str):
        if cdna_start is not None and cdna_end is not None and is_exon_boundary(cdna_start, cdna_end):
            vt = "exon_deletion"
        else:
            size = (cdna_end - cdna_start + 1) if (cdna_start and cdna_end) else 1
            vt = "frameshift" if size % 3 != 0 else "inframe_del"

    elif re.search(r"(dup)", hgvs_str):
        size = (cdna_end - cdna_start + 1) if (cdna_start and cdna_end) else 1
        # Duplication logic: must be >=1 exon and completely within gene
        if duplication_type == "tandem":
            vt = "duplication_frameshift" if size % 3 != 0 else "duplication_inframe"
        elif duplication_type == "not_in_tandem":
            vt = "duplication_not_in_tandem"
        else:
            vt = "duplication_unknown"
    elif re.search(r"(fs)", hgvs_str) or re.search(r"(ins)", hgvs_str):
        vt = "frameshift"
    elif re.search(r"([A,T,G,C]>)", hgvs_str):
        if "*" in hgvs_str or "Ter" in hgvs_str or "stop" in hgvs_str:
            vt = "nonsense"
        else:
            vt = "missense"

    # Initiation codon logic
    if initiation_codon is None:
        if re.search(r"(Met1|M1)", hgvs_str):
            initiation_codon = "Met1"
        elif re.search(r"(Met54|M54)", hgvs_str):
            initiation_codon = "Met54"

    # EARLY truncation logic: exclude any variant that starts before codon 54
    if vt in ["nonsense", "frameshift", "duplication_frameshift"]:
        if codon_start is not None and codon_end is not None:
            if codon_start < 54 and codon_end <= 54:
                return {"strength": None, "context": "Truncation/frameshift before codon 54 does not activate NMD and is outside critical domains—VHL PVS1 not scored."}
            if codon_start < 54 and codon_end > 54:
                return {"strength": None, "context": "Truncation event starts prior to codon 54—PVS1 not scored per VHL guidance."}

    # Exon deletion logic (all exons scored PVS1)
    if vt == "exon_deletion":
        return {"strength": "PVS1", "context": "Any exon deletion in VHL; all exons are critical, assign PVS1 (Very Strong)."}

    # Frameshift/nonsense/duplication logic after codon 54
    if vt in ["nonsense", "frameshift", "duplication_frameshift"]:
        if nmd:
            if critical_domain:
                return {"strength": "PVS1", "context": "Truncation after codon 54, NMD predicted, within a critical domain (AA 63-204)—assign PVS1 (Very Strong)."}
            elif cterm:
                return {"strength": "PVS1_Moderate", "context": "Truncation or NMD in minimal function region (AA 205-213)—assign PVS1_Moderate."}
            else:
                return {"strength": "PVS1", "context": "Truncation after codon 54, NMD predicted, region outside defined critical domains—assign PVS1 (Very Strong)."}
        else:
            if cterm:
                return {"strength": "PVS1_Moderate", "context": "Truncation beyond codon 54, not predicted for NMD, minimal function region (AA 205-213)—assign PVS1_Moderate."}
            else:
                return {"strength": "PVS1_Moderate", "context": "Truncation beyond codon 54, not predicted for NMD, outside critical domains—assign PVS1_Moderate."}

    # In-frame deletions or insertions – not scored
    if vt in ["inframe_del", "duplication_inframe"]:
        return {"strength": None, "context": "In-frame deletion or duplication; PVS1 not applied."}

    # Duplication not in tandem (proven or presumed)
    if vt == "duplication_not_in_tandem":
        return {"strength": "PVS1_Strong", "context": "Reading frame presumed disrupted and NMD predicted to occur (AA 55–136), not in tandem—assign PVS1_Strong."}
    if vt == "duplication_unknown":
        return {"strength": None, "context": "No or unknown impact on reading frame and NMD—PVS1 not applied."}

    # Cryptic splice enhanced logic
    if vt == "cryptic_splice" and min_codon:
        if cryptic_disrupts_rf is True:
            if nmd:
                if critical_domain:
                    return {"strength": "PVS1", "context": "Cryptic splice, frame disrupted, NMD predicted, in critical domain (AA 63-204)—assign PVS1."}
                elif 55 <= min_codon <= 62:
                    return {"strength": "PVS1_Strong", "context": "Cryptic splice, frame disrupted, NMD, AA 55-62—assign PVS1_Strong."}
                elif cterm:
                    return {"strength": "PVS1_Moderate", "context": "Cryptic splice, frame disrupted, NMD, AA 205-213—assign PVS1_Moderate."}
            else:
                if critical_domain or 55 <= min_codon <= 62:
                    return {"strength": "PVS1_Strong", "context": "Cryptic splice, frame disrupted, not NMD, critical domain or AA 55-62—assign PVS1_Strong."}
                elif cterm:
                    return {"strength": "PVS1_Supporting", "context": "Cryptic splice, frame disrupted, not NMD, AA 205-213—assign PVS1_Supporting."}
        if cryptic_preserves_rf is True:
            if critical_domain:
                return {"strength": "PVS1_Strong", "context": "Cryptic splice, preserves frame, critical domain (AA 63-204)—assign PVS1_Strong."}
            elif cterm or (55 <= min_codon <= 62):
                return {"strength": "PVS1_Moderate", "context": "Cryptic splice, preserves frame, AA 205-213 or AA 55-62—assign PVS1_Moderate."}
        return {"strength": "PVS1_Supporting", "context": "Cryptic splice with insufficient RF/NMD info—assign PVS1_Supporting."}

    # Start codon
    if initiation_codon is not None:
        if initiation_codon == "Met1":
            return {"strength": None, "context": "Met1 loss does not impact p19 isoform, PVS1 not applied."}
        elif initiation_codon == "Met54":
            return {"strength": "PVS1", "context": "Met54 loss truncates all isoforms prior to domains—assign PVS1."}

    # Missense or synonymous – not LoF
    if vt == "missense":
        return {"strength": None, "context": "Missense variant—VHL PVS1 applies only to loss-of-function variants."}
    if re.search(r"(synonymous|silent)", hgvs_str):
        return {"strength": None, "context": "Synonymous (silent) change; does not affect protein function—VHL PVS1 not scored."}

    # Fallback
    details = []
    if vt is None:
        details.append("Could not determine variant type from HGVS string.")
    if min_codon is None:
        details.append("Could not extract coding position from HGVS string.")
    if not details:
        details.append("Variant does not match any pathogenic loss-of-function mechanism scored by VHL PVS1 decision tree.")
    context = "; ".join(details)
    return {"strength": None, "context": context}