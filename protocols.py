"""OncoAgent-GBM: open-access lab-protocol templates (grounded, no private SOPs).

These are GENERIC, textbook-level assay workflows written as TEMPLATES. Every
quantitative parameter (seeding density, drug concentrations, timepoints,
incubation, replicate count) is a PLACEHOLDER tagged 'verify against your SOP' —
nothing numeric is asserted as fact, and no private/mentor protocol is used.

Each entry carries:
  * `oa_refs`  : REAL, commonly-cited OPEN-ACCESS method papers (verify the
                 PMID/DOI yourself — the app also fetches live PMC results).
  * `pmc_query`: a query the app runs LIVE against PubMed/PMC so the citations
                 shown are real by construction, not fabricated.

Use: adapt the template to YOUR validated SOP, fill the placeholders, and cite
the open-access source you actually followed (after verifying it).
"""

from __future__ import annotations

VERIFY = "[verify against your SOP — not asserted by this tool]"

# REAL landmark open-access method papers (verify PMID/DOI before citing).
_SCRATCH_REF = ("Liang CC, Park AY, Guan JL. In vitro scratch assay: a "
                "convenient and inexpensive method for analysis of cell "
                "migration in vitro. Nat Protoc. 2007;2(2):329-333. "
                "doi:10.1038/nprot.2007.30. PMID: 17406593.")
_CLONO_REF = ("Franken NA, Rodermond HM, Stap J, Haveman J, van Bree C. "
              "Clonogenic assay of cells in vitro. Nat Protoc. "
              "2006;1(5):2315-2319. doi:10.1038/nprot.2006.339. PMID: 17406473.")


PROTOCOL_LIBRARY = {
    "3D spheroid (perimeter / area time-course)": {
        "summary": ("Grow uniform multicellular spheroids, treat, and track "
                    "perimeter/area over time as a 3D growth/invasion readout."),
        "steps": [
            f"Seed single-cell suspension into ULA round-bottom plates at {VERIFY} cells/well.",
            f"Centrifuge/incubate to form one spheroid per well over {VERIFY} h.",
            f"Confirm a single compact spheroid before treatment (brightfield).",
            f"Add treatment (NSC/TX) at {VERIFY} concentrations; include vehicle control.",
            f"Image at fixed timepoints {VERIFY} with the SAME magnification/settings.",
            "Segment each spheroid; measure area and perimeter in ImageJ/Fiji.",
            f"Normalise to t0 per spheroid; use {VERIFY} spheroids per condition as replicates.",
            "Report mean ± SD/SEM; keep technical vs biological replicates separate.",
        ],
        "oa_refs": [],
        "pmc_query": "glioblastoma spheroid invasion assay protocol",
    },
    "Dose-response cell viability (OD / % survival)": {
        "summary": ("Colorimetric/metabolic viability (MTT, resazurin/alamarBlue) "
                    "across a drug dose range to estimate relative survival."),
        "steps": [
            f"Seed cells in 96-well plates at {VERIFY} cells/well; let adhere {VERIFY} h.",
            f"Prepare a dose series {VERIFY} (log-spaced) + vehicle + blank wells.",
            f"Treat for {VERIFY} h.",
            f"Add viability reagent ({VERIFY}); incubate {VERIFY} h.",
            "Read absorbance/fluorescence; SUBTRACT blank wells.",
            "Normalise to vehicle control = 100%; label the axis 'absorbance' or '% of control', not 'viability', unless validated.",
            f"Use {VERIFY} technical replicates and ≥ the number of independent biological repeats your stats need.",
            "Fit % vs log[dose] (4-parameter logistic) for IC50; if the curve plateaus above the top dose, report 'IC50 > max dose' — do not extrapolate.",
        ],
        "oa_refs": [],
        "pmc_query": "MTT resazurin cell viability dose response protocol",
    },
    "Scratch / wound-healing (migration)": {
        "summary": ("2D migration: a scratch is made in a confluent monolayer and "
                    "gap closure is tracked over time."),
        "steps": [
            f"Grow cells to a confluent monolayer in {VERIFY}-well plates.",
            f"(Optional) pre-treat with a mitosis inhibitor to separate migration from proliferation {VERIFY}.",
            "Make a straight scratch with a sterile pipette tip; wash off debris with PBS.",
            f"Add treatment at {VERIFY} concentrations + vehicle control.",
            f"Image the SAME field at fixed timepoints {VERIFY} (mark reference points).",
            "Measure wound area (ImageJ/Fiji wound-healing tool) at each time.",
            "Report % wound closure vs t0; note whether proliferation was controlled for.",
            f"Replicates: {VERIFY} fields/well, {VERIFY} independent experiments.",
        ],
        "oa_refs": [_SCRATCH_REF],
        "pmc_query": "in vitro scratch wound healing migration assay protocol",
    },
    "Clonogenic (colony-formation)": {
        "summary": ("Long-term proliferative/survival capacity: single cells form "
                    "colonies after treatment; colonies are counted."),
        "steps": [
            f"Seed a low, counted number of single cells {VERIFY} per well/dish.",
            f"Treat with drug {VERIFY} (and/or irradiate) + vehicle control.",
            f"Incubate undisturbed {VERIFY} days until control colonies have ≥ ~50 cells.",
            "Fix and stain (e.g. crystal violet); call stain intensity 'absorbance', not 'viability'.",
            "Count colonies of ≥ 50 cells (manual or ImageJ ColonyArea).",
            "Compute plating efficiency (PE) and surviving fraction (SF) relative to control.",
            f"Replicates: {VERIFY} dishes/condition, {VERIFY} independent experiments.",
            "Report mean ± SD/SEM; one-way ANOVA + Dunnett vs control if multiple doses.",
        ],
        "oa_refs": [_CLONO_REF],
        "pmc_query": "clonogenic colony formation assay protocol glioma",
    },
    "In silico: docking + molecular dynamics methodology": {
        "summary": ("Reproducible structure-based docking and MD notes — define "
                    "the site, validate by redocking, report everything."),
        "steps": [
            "Get the receptor from RCSB PDB; record the PDB ID and resolution.",
            "Prepare: remove waters/ions as justified, add H, assign protonation at the stated pH, keep metals if catalytic.",
            "DEFINE the binding site (co-crystal ligand or catalytic residues) — avoid a whole-protein blind box.",
            "VALIDATE by redocking the native ligand; report RMSD ≤ 2.0 Å before trusting predictions.",
            f"Dock with a FIXED random seed; report seed, Vina version, exhaustiveness, box center/size {VERIFY}.",
            "Treat scores as RELATIVE rankings (kcal/mol), not experimental affinities or selectivity.",
            f"MD (if run): report force field, water model, box, ions, equilibration, production length {VERIFY}; a flat FEL ≠ proven stable binding.",
            "Keep computational observations separate from biological interpretation.",
        ],
        "oa_refs": [],
        "pmc_query": "molecular docking validation redocking RMSD protocol",
    },
}


def list_protocols() -> list[str]:
    return list(PROTOCOL_LIBRARY.keys())


def get_protocol(name: str) -> dict | None:
    return PROTOCOL_LIBRARY.get(name)
