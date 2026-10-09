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

from urllib.parse import quote_plus

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
    "MTT / metabolic viability (endpoint)": {
        "summary": ("Tetrazolium (MTT/MTS/XTT) reduction as a metabolic-activity "
                    "readout; a surrogate for relative viable-cell number."),
        "steps": [
            f"Seed cells in 96-well plates at {VERIFY} cells/well; let adhere {VERIFY} h.",
            f"Treat (NSC/TX) at {VERIFY} concentrations + vehicle + cell-free blanks.",
            f"Incubate {VERIFY} h.",
            f"Add MTT/MTS/XTT reagent ({VERIFY}); incubate {VERIFY} h protected from light.",
            "For MTT: remove medium, dissolve formazan in DMSO/solubiliser; read ~570 nm (ref ~630 nm).",
            "Subtract blanks; normalise to vehicle = 100%. Call it 'metabolic activity / absorbance', not 'viability', unless validated against a direct count.",
            f"Replicates: {VERIFY} technical wells and >= the independent biological repeats your stats need.",
        ],
        "oa_refs": [],
        "pmc_query": "MTT MTS tetrazolium metabolic viability assay protocol",
    },
    "Western blot (protein expression)": {
        "summary": ("Semi-quantitative detection of a protein (e.g. CDC25A/B/C, "
                    "cell-cycle / apoptosis markers) by SDS-PAGE + immunoblot."),
        "steps": [
            f"Lyse cells in cold lysis buffer + protease/phosphatase inhibitors; clarify by centrifugation {VERIFY}.",
            f"Quantify protein (BCA/Bradford); load equal amounts {VERIFY} ug/lane.",
            f"Separate by SDS-PAGE ({VERIFY}% gel) and transfer to PVDF/nitrocellulose.",
            f"Block ({VERIFY}); incubate primary antibody {VERIFY} dilution, overnight 4 C.",
            f"Wash; incubate HRP secondary {VERIFY}; develop by ECL.",
            "Probe a loading control (e.g. GAPDH / beta-actin / total protein stain) on the SAME membrane.",
            "Quantify by densitometry; normalise target to the loading control; report mean +/- SD/SEM.",
            f"Replicates: {VERIFY} independent biological experiments (not just re-exposures of one blot).",
            "Keep full uncropped blots for the thesis annex.",
        ],
        "oa_refs": [],
        "pmc_query": "western blot immunoblot protocol quantification loading control",
    },
    "qRT-PCR (gene expression)": {
        "summary": ("Relative mRNA quantification by reverse transcription + "
                    "real-time PCR, normalised to reference genes."),
        "steps": [
            "Extract total RNA; check integrity and purity (A260/280 ~2.0); DNase-treat if needed.",
            f"Reverse-transcribe {VERIFY} ng RNA to cDNA (include a no-RT control).",
            f"Run qPCR with validated primers {VERIFY}; include no-template controls; run in technical triplicate.",
            f"Pick >= 2 stable reference genes {VERIFY} (validate stability, do not assume).",
            "Confirm single products (melt curve / gel) and acceptable efficiency.",
            "Quantify by 2^-ddCt (or efficiency-corrected); state the method and reference genes used.",
            f"Replicates: technical triplicate within {VERIFY} independent biological repeats.",
        ],
        "oa_refs": [],
        "pmc_query": "qRT-PCR relative gene expression delta delta Ct reference gene protocol",
    },
    "Flow cytometry: cell-cycle (PI)": {
        "summary": ("DNA-content cell-cycle distribution (G0/G1, S, G2/M) by "
                    "propidium-iodide staining of fixed cells — relevant to CDC25 / G2-M."),
        "steps": [
            f"Treat cells {VERIFY}; harvest including floating cells; wash in PBS.",
            "Fix in cold 70% ethanol, added dropwise while vortexing; store {VERIFY}.",
            "Wash; treat with RNase A; stain with propidium iodide.",
            "Acquire on a flow cytometer using a LINEAR scale; gate out doublets (area vs width).",
            "Model the DNA histogram (e.g. in the cytometer software) into G0/G1, S, G2/M fractions.",
            "Report % per phase; a sub-G1 peak is suggestive of fragmentation, not proof of apoptosis on its own.",
            f"Replicates: {VERIFY} independent biological experiments.",
        ],
        "oa_refs": [],
        "pmc_query": "cell cycle analysis propidium iodide flow cytometry protocol",
    },
    "Flow cytometry: apoptosis (Annexin V / PI)": {
        "summary": ("Distinguish viable / early-apoptotic / late-apoptotic-necrotic "
                    "cells by Annexin V-FITC + PI co-staining."),
        "steps": [
            f"Treat cells {VERIFY}; harvest ALL cells (keep floaters — do not use trypsin steps that strip the signal if avoidable).",
            "Wash; resuspend in Annexin-V binding buffer.",
            f"Stain with Annexin V-FITC + PI {VERIFY}; incubate in the dark at RT.",
            "Include single-stain and unstained controls for compensation/gating; a positive apoptosis inducer is a good control.",
            "Acquire promptly; quadrant-gate: AnnV-/PI- viable, AnnV+/PI- early, AnnV+/PI+ late, AnnV-/PI+ necrotic/debris.",
            "Report % per quadrant; interpret cautiously — this is a snapshot, not a mechanism.",
            f"Replicates: {VERIFY} independent biological experiments.",
        ],
        "oa_refs": [],
        "pmc_query": "annexin V propidium iodide apoptosis flow cytometry protocol",
    },
    "Caspase-3/7 activity": {
        "summary": ("Luminescent/fluorescent caspase-3/7 cleavage as an apoptosis "
                    "effector readout."),
        "steps": [
            f"Seed cells in opaque/white plates {VERIFY}; treat {VERIFY} with a positive control.",
            f"Add caspase-3/7 substrate reagent ({VERIFY}); incubate {VERIFY}.",
            "Read luminescence/fluorescence; include no-cell and vehicle controls.",
            "Normalise to cell number / viability where appropriate (a drop in cells can confound raw signal).",
            "Confirm with an orthogonal apoptosis readout before claiming apoptosis.",
            f"Replicates: {VERIFY} technical wells, {VERIFY} independent experiments.",
        ],
        "oa_refs": [],
        "pmc_query": "caspase 3 7 activity apoptosis assay protocol",
    },
    "Immunofluorescence / ICC": {
        "summary": ("Localise/visualise a protein in fixed cells by antibody staining "
                    "and fluorescence microscopy."),
        "steps": [
            f"Grow cells on coverslips/chamber slides; treat {VERIFY}.",
            f"Fix ({VERIFY}, e.g. 4% PFA); permeabilise if the target is intracellular; block.",
            f"Incubate primary antibody {VERIFY}; wash; incubate fluorophore secondary {VERIFY}.",
            "Counterstain nuclei (e.g. DAPI); mount with anti-fade.",
            "Include a no-primary control; image all conditions with IDENTICAL acquisition settings.",
            "Quantify objectively (ImageJ/Fiji) if making quantitative claims; do not eyeball.",
            f"Replicates: {VERIFY} fields and {VERIFY} independent experiments.",
        ],
        "oa_refs": [],
        "pmc_query": "immunofluorescence immunocytochemistry staining protocol",
    },
    "Transwell invasion / migration (Boyden chamber)": {
        "summary": ("Directed migration (bare membrane) or invasion (Matrigel-coated) "
                    "through a porous insert toward a chemoattractant."),
        "steps": [
            "For invasion: coat the insert membrane with a matrix (e.g. Matrigel) and let it set; migration uses a bare membrane.",
            f"Seed cells in serum-free medium in the upper chamber {VERIFY}; add chemoattractant (e.g. serum) below.",
            f"Add treatment {VERIFY} + vehicle control.",
            f"Incubate {VERIFY}; remove non-migrated cells from the upper side with a swab.",
            "Fix and stain cells on the lower side; image multiple random fields.",
            "Count migrated/invaded cells (manual or ImageJ); control for any proliferation difference over the assay window.",
            f"Replicates: {VERIFY} inserts/condition, {VERIFY} independent experiments.",
        ],
        "oa_refs": [],
        "pmc_query": "transwell Boyden chamber invasion migration assay protocol",
    },
    "Comet assay (DNA damage)": {
        "summary": ("Single-cell gel electrophoresis measuring DNA strand breaks; "
                    "damaged DNA migrates as a comet tail."),
        "steps": [
            f"Treat cells {VERIFY}; embed in low-melting agarose on slides.",
            "Lyse; choose alkaline (single + double-strand breaks) or neutral (double-strand) conditions and STATE which.",
            f"Electrophorese {VERIFY}; neutralise; stain DNA.",
            "Score >= 50-100 comets/sample with validated software (e.g. % tail DNA, tail moment).",
            "Include a positive control (e.g. H2O2) and a vehicle control.",
            f"Replicates: {VERIFY} independent experiments; score blinded if possible.",
        ],
        "oa_refs": [],
        "pmc_query": "comet assay single cell gel electrophoresis DNA damage protocol",
    },
    "EdU / BrdU proliferation": {
        "summary": ("Label actively replicating (S-phase) cells by thymidine-analogue "
                    "incorporation to measure proliferation directly."),
        "steps": [
            f"Treat cells {VERIFY}; pulse with EdU/BrdU {VERIFY} before harvest.",
            "EdU: fix, permeabilise, run the click reaction with the fluorescent azide. BrdU: fix, denature DNA, stain with anti-BrdU.",
            "Counterstain total DNA/nuclei; include a no-label control.",
            "Read by microscopy or flow; report % labelled (S-phase) cells vs control.",
            f"Replicates: {VERIFY} independent experiments.",
        ],
        "oa_refs": [],
        "pmc_query": "EdU BrdU incorporation proliferation S phase assay protocol",
    },
    "In silico ADMET pre-screen (SwissADME / ProTox)": {
        "summary": ("Computational physchem / drug-likeness / early tox flags on a "
                    "candidate (e.g. NSC95397) — a prioritisation filter, not proof."),
        "steps": [
            "Get a clean, correct structure (canonical SMILES); record the source and any tautomer/charge choice.",
            "Run SwissADME: physchem, Lipinski/Veber, BOILED-Egg (GI absorption / BBB), solubility, CYP flags.",
            "Run ProTox (or equivalent): predicted toxicity class / endpoints — record model + confidence.",
            "Report predictions as FLAGS to prioritise, never as measured ADMET or safety.",
            "State the exact tool version/date; a prediction is not an experimental result.",
            "Keep computational observations separate from any biological interpretation.",
        ],
        "oa_refs": [],
        "pmc_query": "SwissADME ProTox in silico ADMET prediction drug-likeness",
    },
}


def list_protocols() -> list[str]:
    return list(PROTOCOL_LIBRARY.keys())


def get_protocol(name: str) -> dict | None:
    return PROTOCOL_LIBRARY.get(name)


# ============================================================
# "Name any protocol -> get real sources" finder
# ============================================================
def protocol_pubmed_query(name: str) -> str:
    """Build a PubMed query that favours real METHOD/PROTOCOL papers for whatever
    assay/experiment the user types. No fabrication — this only shapes a search.
    """
    name = (name or "").strip()
    if not name:
        return ""
    return (f'("{name}"[Title/Abstract]) AND '
            f'(protocol[Title/Abstract] OR method*[Title/Abstract] OR '
            f'assay[Title/Abstract] OR procedure[Title/Abstract])')


def protocol_source_links(name: str) -> list[tuple[str, str]]:
    """Return REAL search-engine links for the typed protocol name across the
    main protocol/method repositories. These are deterministic public search
    URLs — nothing invented; each opens the live results for the user's term.
    """
    q = quote_plus((name or "").strip())
    if not q:
        return []
    return [
        ("PubMed", f"https://pubmed.ncbi.nlm.nih.gov/?term={q}+protocol"),
        ("Europe PMC (open access)", f"https://europepmc.org/search?query={q}%20protocol"),
        ("PMC (free full text)", f"https://www.ncbi.nlm.nih.gov/pmc/?term={q}+protocol"),
        ("protocols.io", f"https://www.protocols.io/search?q={q}"),
        ("Bio-protocol", f"https://bio-protocol.org/en/search?q={q}"),
        ("Nature Protocols", f"https://www.nature.com/search?q={q}&journal=nprot"),
        ("Springer Protocols", f"https://experiments.springernature.com/search?q={q}"),
    ]
