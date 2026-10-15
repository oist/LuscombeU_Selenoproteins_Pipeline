# Agent instructions: selenocysteine (Sec) annotation of Oikopleura genomes

These instructions tell an agent how to annotate selenoprotein genes and their Sec codons in a new _Oikopleura_ genome assembly, using the pipeline, models and checks developed on _O. dioica_ OKI2018_I68 (GCF_907165135.1) and six diploid assemblies. Follow them in order. Report what was checked and what was not.

## 1. Goal and deliverables

For each genome, produce:

1. A Sec codon layer (GFF3, one `stop_codon_redefined_as_selenocysteine` feature per Sec codon) with the evidence for each codon. This is the main product: annotation pipelines should take Sec codons as high-scoring evidence, not Selenoprofiles gene models (Ticó et al. 2026, PLoS Comput Biol 22:e1013885, Fig 5).
2. A table of the selenoprotein genes found, by family, compared with the expected _Oikopleura_ set.
3. If an annotation exists: the status of each Sec codon in it (classes of Ticó et al. 2026) and a list of corrections.
4. A short report stating the checks done, the genes still unresolved, and why.

## 2. Resources

All in the `oist/LuscombeU_Selenoproteins_Pipeline` repository (branch `selenoprotein-annotation`) unless noted.

| Resource                        | Path                                                                                                                                      | Use                                                                                                                                                                                                    |
| ------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Pipeline                        | repository root                                                                                                                           | Selenoprofiles, miniprot cross-check, SECIS, tRNA-Sec, annotation checks                                                                                                                               |
| _O. dioica_ SECIS model         | `assets/secis/oikopleura_dioica_secis.cm`                                                                                                 | SECIS search; no gathering threshold, always use a bit-score cut-off                                                                                                                                   |
| Expectation                     | `assets/expectations/oikopleura_dioica.csv`, lineage `Oikopleura`                                                                         | expected genes per family (17 per haploid genome)                                                                                                                                                      |
| Reference Sec codons            | `OKI2018_Sec_codons.gff3` (outputs of this project)                                                                                       | 16 Sec codons of OKI2018_I68, with evidence                                                                                                                                                            |
| Appendicularian SECIS model     | `assets/secis/appendicularia_secis.cm`                                                                                                    | 3′ regions in other appendicularians and of known genes (held out: 69 of 106 genes at 20 bits, against 50 with the _O. dioica_ model); genome-wide in _O. dioica_ its background reaches about 34 bits |
| Merged Rfam + _O. dioica_ model | `RF00031_plus_Odioica.cm` (outputs)                                                                                                       | second SECIS opinion; scores the SELENOW PAR element higher than the _O. dioica_ model (21.9 against 11.3 bits)                                                                                        |
| Scripts                         | `bin/readthrough_screen.py`, `bin/rank_readthrough.py`, `bin/sec_annotation_check.py`, `bin/expectation_check.py`, `bin/compare_genes.py` | usable outside the pipeline                                                                                                                                                                            |

## 3. Expected selenoproteome of _Oikopleura dioica_

17 genes per haploid genome, in 9 families. Positions are OKI2018_I68; in other assemblies genes stay on the same chromosome arm but may move by up to about 0.8 Mb.

| Gene                | Location (OKI2018) | Sec codon                             | Found by       | Notes                                                                                                     |
| ------------------- | ------------------ | ------------------------------------- | -------------- | --------------------------------------------------------------------------------------------------------- |
| DIO (chr1)          | NC_137291.1        | 3,441,967–969 (−)                     | manual         | missed in OKI2018 by Selenoprofiles; found with the OSKA2016 protein; RefSeq protein starts after the Sec |
| DIO (PAR)           | NC_137293.1        | 9,648,087–089 (−)                     | both           |                                                                                                           |
| GPX (chr1-a)        | NC_137291.1        | 4,438,200–202 (−)                     | Selenoprofiles | RefSeq translates the Sec codon as X                                                                      |
| GPX (chr1-b)        | NC_137291.1        | 11,590,011–013 (+)                    | Selenoprofiles | RefSeq protein starts after the Sec                                                                       |
| GPX (chr2)          | NC_137292.1        | 9,603,099–101 (+)                     | both           | RefSeq translates the Sec codon as X                                                                      |
| GPX (PAR)           | NC_137293.1        | 3,900,591–593 (+)                     | both           |                                                                                                           |
| SELENOF             | NC_137291.1        | 12,150,879–881 (+)                    | Selenoprofiles | Selenoprofiles misses a 40-nt intron; RefSeq model correct                                                |
| SELENOO             | NC_137293.1        | 2,223,744–746 (−)                     | miniprot       | C-terminal Sec (…KVT-U-SS*); Selenoprofiles labels it unaligned                                           |
| SELENOS             | PAR (LOC145705538) | — in OKI2018                          | SECIS search   | C-terminal Sec (…RRGPSRGG-U-G*); OKI2018 has a 13-nt deletion that shifts the frame                       |
| SELENOT             | NC_137294.1 (XSR)  | 3,315,985–987 (+)                     | manual         | Sec exon (…AICTS-U-G) missed by Selenoprofiles and RefSeq; last intron is GA–AG                           |
| SELENOW (chr1-a)    | NC_137291.1        | 11,963,712–714 (−)                    | both           | no RefSeq model                                                                                           |
| SELENOW (chr1 pair) | NC_137291.1        | 13,732,976–978 and 13,733,472–474 (+) | Selenoprofiles | tandem pair about 200 nt apart; no RefSeq model (inside the 3′UTR of LOC145695605)                        |
| SELENOW (chr2)      | NC_137292.1        | 6,496,264–266 (+)                     | Selenoprofiles | RefSeq protein starts after the Sec                                                                       |
| SELENOW (PAR)       | NC_137293.1        | 16,697,745–747 (+)                    | Selenoprofiles | no RefSeq model (inside the 3′UTR of an eEF2K-like gene)                                                  |
| SEPHS2              | NC_137293.1 (PAR)  | 5,008,527–529 (−)                     | manual         | Sec exon (…GLKG-U-GCK) missed by Selenoprofiles and RefSeq                                                |
| TXNRD               | NC_137294.1 (XSR)  | 6,205,537–539 (−)                     | both           | C-terminal Sec (…GC-U-G*); RefSeq protein ends at the Sec                                                 |

SELENOK (the 18th gene, 10th family): XSR, Sec codon NC_137294.1:7,274,834–836 (+), three exons from 7,274,477, 94 aa ending …NMGVGGG-U-GR, SECIS 64 nt after the stop (28 bits); 26 nt downstream of a ciliogenesis and planar polarity effector 2-like gene (operon-like), not annotated by RefSeq and missed by Selenoprofiles and miniprot. TGA present in OSKA2016 (with a transcript, GSE64421) and OdB3; SELENOK with the same C-terminus in O. albicans, O. vanhoeffeni and O. longicauda.

Candidate 19th gene: thioredoxin-like protein 1 (OKI2018_I68 LOC145697077, chr2 about 13.04 Mb), C-terminus …KKVG-U-G\*, with the TGA present in OKI2018_I68, OSKA2016 and OdB3, and a SECIS-like element 78 nt after it (27–28 bits with the appendicularian model only). _O. albicans_ and _O. vanhoeffeni_ have Thr or Ala at that position. Needs transcript or protein evidence.

Also expected: one tRNA-Sec gene on PAR (G73), eEFsec, SecS, SBP2. PSTK is present but divergent: in OKI2018_I68 it is LOC145708732 on XSR (about 1.94 Mb, + strand), which RefSeq truncates to 110 aa; the OKI2018 transcriptome (Bliznina et al. 2021) gives five coding exons (1,940,258–381, 425–554, 598–636, 700–962, 1,941,026–195) and a 241-aa protein matching the whole PSTK profile (hmmsearch E 8e-25); an OSKA2016 transcript (GSE64421) encodes a 246-aa PSTK (E 3e-23). Selenoprofiles, miniprot and tblastn miss it; search six-frame ORFs with a PSTK profile HMM and tell it from KTI12 (LOC145708157, XSR 1.47 Mb) by comparing alignment scores to both. Subfamily names from `selenoprofiles orthology` (vertebrate anchors) are tentative: GPX4 ×4 (supported by a tree), SELENOW1 ×4 and SELENOW2 ×1, DIO1 ×2 and TXNRD3 ×1 (trees disagree). Not selenoproteins in _O. dioica_: the SELENOW-like UGA fragment on PAR (8 of 12 haplotypes), SELENOK hits (proline-rich spurious matches), AhpC (Thr/Cys), the chr1 SEPHS paralog (Ser), MSRB, MsrA, SELENOU (Cys forms).

## 4. Procedure

### Step 1. Prepare the inputs

- Genome FASTA. For diploid assemblies, note which sequences are the second haplotype (`_alt`, `_PAT`) and which are sex-specific (`XSR`, `YSR`) or mitochondrial (`chrM`, `chrMT`). Ignore any Sec call on mitochondrial sequences (TGA is tryptophan there).
- Annotation, if any, as GFF3 with mRNA (or transcript), exon and CDS features. A GFF3 without exon features cannot be screened for transcripts. To convert a GenBank flat file, keep mRNA exons and `transl_except`.

### Step 2. Run the pipeline

```bash
nextflow run oist/LuscombeU_Selenoproteins_Pipeline -r selenoprotein-annotation -profile <institution> \
  --input samplesheet.csv --outdir results \
  --selenoprofiles_data <selenoprofiles_data_dir> \
  --secis_cm ${projectDir}/assets/secis/oikopleura_dioica_secis.cm --secis_min_bits 15 \
  --readthrough_screen \
  --expectation ${projectDir}/assets/expectations/oikopleura_dioica.csv --expectation_lineage Oikopleura
# diploid assembly of a male: add --expectation_copies 2 --expectation_hemizygous TXNRD,SELENOT
```

Samplesheet: `sample,fasta,domain,profiles,annotation` with `domain` = `eukaryota` and the GFF3 in `annotation` when available. Keep the default profile sets (`eukarya,machinery`).

### Step 3. Read the results in this order

1. `expectation/<sample>.expectation_summary.tsv`. Missing families are expected for SELENOT, SEPHS2 and SELENOS (blind spots of Selenoprofiles 4.5.8). Any other missing or extra family needs Step 4.
2. `compare/<sample>.sec_sites.tsv` and `compare/<sample>.sec_codons.gff3`: the Sec codons, with support (Selenoprofiles, miniprot or both) and SECIS bits and E-value.
3. With an annotation: `annotation/<sample>.sec_annotation.tsv` (class per Sec codon) and `readthrough/<sample>.candidates.tsv` (transcripts ranked by SECIS).
4. `multiqc/multiqc_report.html` for the overview.

### Step 4. Close the known blind spots

For each, record the Sec codon coordinates, the motif, the splice sites and the SECIS score.

- **SELENOT** (X-specific region). Selenoprofiles predicts it as `unaligned` and starts at …YRNAV…. Translate the 3 kb upstream in all frames, with TGA as U, and find the in-frame `C..UG` motif (…AICTS-U-G in OKI2018). It ends an exon of about 68 nt, followed by a GT donor; the next exon begins …YRN…. The gene has 5 exons: introns 1–3 are GT–AG of 44–49 nt, intron 4 is GA–AG. The protein is about 198 aa with an ATG in exon 1.
- **SEPHS2** (PAR). Selenoprofiles predicts it as `unaligned`, starting after the Sec. Look upstream for `G.KG U G` with a TGA (…NFTGLKG-U-GCKVP…), in a first exon ending with a GT donor 46 nt before exon 2.
- **SELENOS** (PAR). Not predicted by Selenoprofiles. Find it by a SECIS hit (about 27 bits with the _O. dioica_ model) about 65–129 nt after the stop codon, and by homology of its N-terminus (…YGWYIV…, transmembrane) to SELENOS. In intact copies (OSKA2016, OdB3, and one OKI2018 allele in its transcriptome) the C-terminus is …RRGPSRGG-U-G*. Check for frameshifts: the OKI2018_I68 assembly carries a frameshifted allele.
- **SELENOO**. Selenoprofiles labels it `unaligned`; miniprot finds the C-terminal Sec (…KVT-U-SS*).
- **chr1 DIO**. Can be missed by Selenoprofiles: search with the protein from another _O. dioica_ assembly (miniprot, or tblastn of the gene ±3 kb, then align the Sec exon) before calling it absent.

For any other family flagged missing, apply the same approach: take the allele or homolog from another assembly, search with blastn or tblastn, align the Sec exon, and read the codon.

### Step 5. Verify each Sec codon

A Sec codon is accepted when:

- it lies in frame in a gene model with a conserved Sec motif of its family; and
- at least one of: miniprot support; the same codon in other haplotypes or assemblies; a SECIS hit in the 3′UTR (15 bits or more with the _O. dioica_ model, or with the merged model, in the first 500 nt after the stop codon); transcript evidence.

Rules for gene structure:

- Allow GA–AG introns (about 12% of _O. dioica_ introns) and very short introns (about 40 nt). When a boundary falls on a non-canonical pair, look for an equivalent placement (shifts that keep the spliced sequence) giving GT–AG or GA–AG.
- When two alleles differ by more than 10% in protein identity, project one allele's exon model onto the other haplotype before concluding: most low values come from different predicted exons.
- A gene missing from one haplotype is not missing from the genome: search that haplotype with the allele from the other one, then search 10–15 kb flanks to distinguish deletion from assembly gap; report gaps that need reads.

### Step 6. Compare with the existing annotation

Use the classes of Ticó et al. 2026 from `sec_annotation_check.py`: Well annotated (Sec declared, or X/other), Stop codon, Skipped, Downstream, Upstream, Out of frame, Spliced, Absent. In OKI2018_I68 RefSeq, only 3 of 17 genes are well annotated with Sec; the dominant errors are Downstream (protein starts after the Sec) and Absent, unlike vertebrate RefSeq where most errors are C-terminal Stop codon cases. Corrections should give Sec codon coordinates and evidence; for SELENOT and SEPHS2 the transcript models are already right and only the CDS start and the Sec exception change.

### Step 7. Screen for genes nobody predicted

With an annotation, `readthrough/<sample>.candidates.tsv` lists transcripts with a TGA at the edges or inside the CDS, ranked by SECIS. In OKI2018_I68 every candidate of 15 bits or more was a known selenoprotein. Without an annotation, search the genome with the _O. dioica_ SECIS model (`cmsearch -T 20`): the known SECIS elements score 33–56 bits; background hits reach about 28 bits, so treat 30 bits or more as candidates and check the gene upstream for a TGA in an open reading frame.

## 5. Interpretation rules

- Absence from Selenoprofiles output is not absence of the gene (SELENOT, SEPHS2, SELENOS, chr1 DIO were all missed somewhere).
- Diploid assemblies carry autosomes and PAR twice and the sex-specific regions once (TXNRD and SELENOT are on the X-specific region).
- Disrupted copies can occur in one assembly (the OKI2018_I68 SELENOS frameshift); compare with other _O. dioica_ assemblies before calling a gene lost.
- The Rfam SECIS model (RF00031, threshold 35 bits) misses _O. dioica_ SECIS elements; do not use its gathering threshold here.
- `selenoprofiles lineage` cannot be used: it needs the `ncbi_db` package and its expectation table covers vertebrates only. Use the pipeline's expectation check instead. `selenoprofiles orthology` needs `pyaln` and gives vertebrate subfamily names that are tentative for tunicates.
- Lowering exonerate's minimum intron length does not recover the missed Sec exons; do not spend time on it.

## 6. Known limitations of the tools

- Container images of the new Selenoprofiles and Secmarker modules are placeholders; the official Selenoprofiles Docker image is `maxtico/selenoprofiles_container` (only a `latest` tag; pin it by digest).
- The new local modules (annotation check, transcript screen, expectation check) have no nf-test tests yet; the test profile does not exercise them.
- No RNA-seq was used to define 3′UTRs or confirm the manual Sec exons; prefer transcript evidence when available.
