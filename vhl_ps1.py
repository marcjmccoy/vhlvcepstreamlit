import re

############################# PS1 ###################################################
def classify_vhl_ps1(test_variant):
    """
    Assess PS1 for VHL using ACMG/AMP VCEP rules.
    Returns: dict with 'strength' (str or None), 'context' (str)
    """

    # Complete list of VHL VCEP pathogenic variants (with IDs)
    pathogenic_variants = [
        {"Preferred Variant Title":"NM_000551.4(VHL):c.191G>C (p.Arg64Pro)", "ClinVar_ID":2226, "CAID":"CA020089"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.263G>A (p.Trp88Ter)", "ClinVar_ID":182978, "CAID":"CA020197"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.208G>T (p.Glu70Ter)", "ClinVar_ID":428806, "CAID":"CA16602179"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.463+1G>A", "ClinVar_ID":526679, "CAID":"CA16621909"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.194C>G (p.Ser65Trp)", "ClinVar_ID":43597, "CAID":"CA020099"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.586A>T (p.Lys196Ter)", "ClinVar_ID":196284, "CAID":"CA020507"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.583C>T (p.Gln195Ter)", "ClinVar_ID":428794, "CAID":"CA70052558"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.500G>A (p.Arg167Gln)", "ClinVar_ID":2216, "CAID":"CA020454"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.477del (p.Glu160fs)", "ClinVar_ID":182959, "CAID":"CA020404"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.422dup (p.Asn141fs)", "ClinVar_ID":411979, "CAID":"CA16611276"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.408del (p.Phe136fs)", "ClinVar_ID":43601, "CAID":"CA020343"},
        {"Preferred Variant Title":"NM_000551.4(VHL):c.341-2A>G", "ClinVar_ID":223194, "CAID":"CA357004"}
    ]

    # Parse protein change if present
    protein_match = re.search(r'\((p\.[A-Za-z]{3}\d+[A-Za-z]{3})\)', test_variant)
    protein_change = protein_match.group(1) if protein_match else None

    # Assign variant type
    vt = None
    # Splice site
    if re.search(r"(\+|-)", test_variant):
        vt = "splice"
    # Frameshift
    elif "fs" in test_variant:
        vt = "frameshift"
    # Nonsense
    elif "Ter" in test_variant or "*" in test_variant or "stop" in test_variant:
        vt = "nonsense"
    # In-frame indels
    elif "del" in test_variant and not "fs" in test_variant:
        vt = "inframe_del"
    elif "dup" in test_variant and not "fs" in test_variant:
        vt = "inframe_dup"
    elif "ins" in test_variant and not "fs" in test_variant:
        vt = "inframe_ins"
    # Missense
    elif protein_change and not ("fs" in protein_change or "Ter" in protein_change):
        vt = "missense"
    else:
        vt = "other"

    # PS1 applies only to missense variants matching a previously established pathogenic missense protein change
    if vt == "missense":
        found = False
        for var in pathogenic_variants:
            ref_match = re.search(r'\((p\.[A-Za-z]{3}\d+[A-Za-z]{3})\)', var["Preferred Variant Title"])
            ref_protein_change = ref_match.group(1) if ref_match else None
            # Only consider missense protein changes from the reference list
            if ref_protein_change and ('fs' not in ref_protein_change and 'Ter' not in ref_protein_change):
                if protein_change == ref_protein_change:
                    found = True
                    return {
                        'strength': 'PS1',
                        'context': (
                            f"Missense variant with amino acid substitution {protein_change} matches previously established VHL VCEP/ClinGen pathogenic missense "
                            f'variant "{var["Preferred Variant Title"]}" [ClinVar_ID: {var["ClinVar_ID"]}, CAID: {var["CAID"]}]. '
                            "PS1 (Strong) applies only if interpretation is by VHL VCEP."
                        )
                    }
        # No match: valid missense, but not established
        return {
            'strength': None,
            'context': (
                f"Missense substitution ({protein_change}) but no matching VHL VCEP/ClinGen pathogenic missense in the reference list; PS1 does not apply."
            )
        }

    # All other types return context per ACMG/AMP PS1 guideline
    elif vt == "splice":
        return {'strength': None,
                'context': "This is a splice site variant; PS1 applies only to missense substitutions, not variants mainly affecting splicing."}
    elif vt == "frameshift":
        return {'strength': None,
                'context': "This is a frameshift variant; PS1 applies only to missense substitutions, not frameshift/loss-of-function variants."}
    elif vt == "nonsense":
        return {'strength': None,
                'context': "This is a nonsense/stop-gain variant; PS1 is not applicable; only scored for missense substitutions."}
    elif vt in ["inframe_del", "inframe_dup", "inframe_ins"]:
        return {'strength': None,
                'context': f"This is an in-frame {vt.split('_')[1]} variant; PS1 is not applicable."}
    else:
        details = []
        if not protein_change:
            details.append("No explicit amino acid change in input HGVS.")
        details.append("PS1 can only be considered for variants matching established pathogenic VHL VCEP/ClinGen variants.")
        return {'strength': None, 'context': "; ".join(details)}