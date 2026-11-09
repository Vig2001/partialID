# Plan: literature review (related work) for the paper

Date: 2026-11-08. Status: **draft, awaiting approval**. Nothing has been
written into the paper yet.

## Goal

A related-work section of about four or five paragraphs, in three parts:

1. intersection-based partial identification and sensitivity analysis (two
   paragraphs);
2. the winner's curse and inference on intersection bounds (one or two
   paragraphs);
3. data fusion for estimation, focusing on convex combinations (one short
   paragraph).

Each part ends with a sentence placing our work against it. The CLeaR version
(due 16 Nov 2026, separate Overleaf project, copied by hand) will probably need
a shorter cut, so every paragraph is planned with a core set of references and
an optional set.

## Where it goes

- `partialID_Overleaf/Literature Review/` already has three files:
  `partial_identification.tex` and `data_fusion.tex` (both empty) and
  `df_sensitivity.tex` (three lines of notes in comments, on Yu et al.,
  Wang et al. and Lanners et al.). None is input in `compiled.tex`, and
  `Introduction/introduction.tex` is empty.
- Proposal:
  - Part 1 goes in `partial_identification.tex`, absorbing the notes in
    `df_sensitivity.tex`.
  - Part 2 goes in a new file, `winners_curse.tex`.
  - Part 3 goes in `data_fusion.tex`.
  - Add `\section{Related Work}` with the three `\input`s to `compiled.tex`,
    after the abstract and before "Convex Combination". This needs your OK
    (open question 1).

## Key to reference status

- **bib**: already in `references.bib`.
- **Z**: in Zotero, not in the `.bib` yet.
- **V**: not in Zotero. I checked that it exists and what it says, against the
  PDF on arXiv or a verified reference list (2026-11-08).
- **M**: from memory, not checked this session. Verify before citing.

Papers you may have missed are marked **[new]**.

---

## Part 1: intersection-based partial identification and sensitivity analysis

### Paragraph 1a: bands from assumptions (background, keep short)

Message: relaxing an untestable assumption gives a band indexed by a
sensitivity parameter. Different datasets, or different assumptions, give
different valid bands for the same causal quantity.

| Reference | Status | Use |
|---|---|---|
| Manski (1990), AER P&P 80 | V **[new]** | origin of nonparametric bounds |
| Horowitz & Manski (2000) | bib | bounds with missing data in RCTs |
| Tamer (2010), Annu. Rev. Econ.; Molinari (2020), Handbook of Econometrics 7A | M **[new]** | reviews of partial identification (cite one) |
| Rosenbaum (1987); Tan (2006); Zhao, Small & Bhattacharya (2019); Dorn & Guo (2023); Dorn, Guo & Kallus (2025) | bib | MSM line; DVDS is what the causal simulation uses |
| Ding & VanderWeele (2016); Bonvini & Kennedy (2022); Yadlowsky et al. (2022) | bib | other sensitivity models (optional) |
| Nguyen et al. (2017, 2018); Nie, Imbens & Wager (2021); Dahabreh et al. (2023); Asiaee et al. (2026) | bib | sensitivity for generalisation and transport (violations of study exchangeability) |

### Paragraph 1b: combining sources by intersecting identified sets

Message: intersecting identified sets from several valid restrictions or
sources is classical in econometrics. Recent causal data fusion does the same
with one sensitivity parameter per source, and that is our setting.

| Reference | Status | Use |
|---|---|---|
| Manski & Pepper (2000), Econometrica 68 | V **[new]** | monotone IV: bounds are a sup/inf over instrument values, the canonical intersection bound |
| Haile & Tamer (2003), JPE | M **[new]** | intersection bounds from an incomplete model (auctions) |
| Cross & Manski (2002), Econometrica; Ridder & Moffitt (2007), Handbook 6B | M **[new]** | partial identification from combining datasets |
| Fan, Sherman & Shum (2014), Econometrica 82 | V **[new]** | treatment effects under data combination |
| Sakaguchi (2025), arXiv 2508.12206 | V **[new]** | combining experimental and observational data sharpens bounds on distributional effects beyond either dataset alone |
| **Lanners et al. (2025)** | bib | closest work: RCT and OS, parameters rho and gamma, Boltzmann smoothing of the max/min |
| **Yu et al. (2025)** | bib (key needs fixing) | two-parameter sensitivity analysis, Fisher combination of p-values |
| Levis et al. (2023), covariate-assisted bounds | Z **[new to bib]** | intersecting pointwise in X gives a sharper target than intersecting at the ATE level (our ATE-level vs covariate-level point) |
| De Bartolomeis et al. (2024); Hussain et al. (2023); Fawkes et al. (2025); Demirel et al. (2025) | bib | optional sentence: the RCT used to test or bound the OS's bias, rather than intersecting |

Positioning sentence: two independent datasets, each with its own assumptions
and a valid band. The target is the intersection of the population bands. Our
methods take each dataset's bound estimates as given, so they do not depend on
the sensitivity model or the estimator. Against Lanners et al., reuse the draft
paragraph in `quality_reports/reviews/2026-11-07_lit-notes-lanners-yu.md`:
- their weights smooth the max so that an influence function exists, while
  ours are chosen to narrow the CI;
- they work at the covariate level, we work at the ATE level;
- for fixed alpha their inference is valid for the smoothed (wider) target. A
  large alpha approaches the hard max, which is the non-regular regime this
  paper is about.

Fixes to the notes in `df_sensitivity.tex`:
- "intersection-based sensitivity analysis was first introduced by Lanners" is
  not defensible given Manski & Pepper (2000) and the data-combination
  literature. Say "the closest work in causal data fusion, where treatments
  and outcomes are observed in both sources".
- "a fixed and finite smoothing parameter is not subject to the post-selective
  inference problem" is fine. Sharpen it: the weight is a fixed smooth function
  of the nuisances, its estimation error enters the influence function, and the
  inference is for the smoothed target.

---

## Part 2: the winner's curse and inference on intersection bounds

### Paragraph 2a: why the estimated intersection is a problem

Message: the max or min of estimated bounds is biased inwards (a winner's
curse) and not differentiable at ties. No regular or unbiased estimator
exists, the naive CI undercovers, and the standard bootstrap fails.

| Reference | Status | Use |
|---|---|---|
| Hirano & Porter (2012), Econometrica 80 | M **[new]** | no locally unbiased or regular estimators for non-differentiable functionals |
| Andrews (2000), Econometrica 68; Fang & Santos (2019), REStud 86 | M **[new]** (Fang & Santos is already in your To_Do) | bootstrap failure, and inference for directionally differentiable maps |
| Andrews, Kitagawa & McCloskey (2024), QJE 139 | bib | inference on winners: conditional and hybrid CIs |
| Zhong & Prentice (2008), Biostatistics | M **[new]** | winner's curse in GWAS (optional, for the name) |
| Berk et al. (2013), AoS; Lee, Sun, Sun & Taylor (2016), AoS; Fithian, Sun & Taylor (2014) | V **[new]** | selective inference: simultaneous (PoSI) and conditional |
| Zrnic & Fithian (2024), *Locally simultaneous inference*, arXiv 2212.09009 | V **[new]** (check the journal) | between simultaneous and conditional; an example is inference on the winner |
| Benjamini, Hechtlinger & Stark (2019), arXiv 1906.00505 | V **[new]** | CIs for selected parameters; for two candidates no correction is needed under conditions (see the caveat below) |
| Luedtke & van der Laan (2016), AoS | M **[new]** | non-regularity of E[max] when ties have positive mass ("exceptional laws"), relevant to the covariate-level point |

Caveat to state clearly in the text: *inference on the selected dataset's own
bound* (the winner's mean) and *inference on the min or max of the bounds*
(the intersection bound) are different targets. AKM, Benjamini et al. and
Han & McCloskey are about the first; CLR and our methods are about the second.
Benjamini et al.'s "no correction for two candidates" result does not carry
over to ours: M1 is the naive min and undercovers (0.898 at r = 1, g = 0).

### Paragraph 2b: inference on intersection bounds, and the remedies

| Reference | Status | Use |
|---|---|---|
| Imbens & Manski (2004) | bib | CIs for the parameter versus for the identified set |
| Stoye (2009), Econometrica 77 | V **[new]** | refinements of Imbens & Manski |
| Bugni, Gao, Obradović & Velez (2024), arXiv 2407.20386 | V **[new]** (optional) | power of these CIs; more efficient bound estimators need not give more power for one of them |
| **Chernozhukov, Lee & Rosen (2013)** | bib | precision-corrected bounds, half-median-unbiased estimates, adaptive inequality selection; M8 and M9 are this without the selection step |
| Romano, Shaikh & Wolf (2014), Econometrica 82; Andrews, Roth & Pakes (2023), REStud 90 | V **[new]** (optional) | moment-inequality inference (two-step and hybrid tests) |
| **Han & McCloskey (2024), arXiv 2403.00422** | V **[new]** | post-selection inference for interval-identified parameters. Their running example shows the estimated max of two bounds is not normal; conditional, projection and hybrid CIs |
| McCloskey (2024), Biometrika 111 | V **[new]** (optional) | hybrid CIs after model selection |
| Cox (1975), Biometrika 62 | M **[new]** | sample splitting for valid inference |
| Chernozhukov et al. (2018), Econometrics J. 21 | V **[new]** | cross-fitting (the origin of M7) |

Positioning sentence: there are four ways to handle a data-dependent choice,
and each pays differently:
- sample splitting pays in sample size (Cox; DML);
- simultaneous critical values pay in the critical value (CLR; PoSI);
- conditional or hybrid inference pays in complexity and instability
  (AKM; McCloskey; Han & McCloskey);
- smoothing pays in a wider target (Lanners).

Our contribution covers the continuum of convex combinations between the two
bands, with weights chosen to minimise CI width. It is valid by splitting (M3)
or by a sup over all weights (M9, which extends CLR), with cross-fitting (M7)
as an extension. Width comparisons follow the Corrections in CLAUDE.md: no
comparison with methods that undercover.

---

## Part 3: data fusion for estimation and convex combinations (one paragraph)

Message: in point-identified fusion, a convex combination of an unbiased and a
possibly biased estimator trades bias against variance for one target. The
weight is data-dependent, so valid inference is hard: without a bound on the
bias, any CI valid for all biases is about as long as the unbiased estimator's
own CI.

| Reference | Status | Use |
|---|---|---|
| Colnet et al. (2024), Stat. Sci. 39; Lin, Tarp & Evans (2024); Degtiar & Rose (2023) | bib | reviews (Lin et al. is what `convcomb.tex` cites as inspiration) |
| Bates & Granger (1969); Green & Strawderman (1991), JASA 86 | V **[new]** | classical roots: forecast combination; James–Stein combination of unbiased and biased estimators |
| Rosenman, Basse, Owen & Baiocchi (2023), Biometrics 79 | Z **[new to bib]** | shrinkage estimators |
| Rosenman, Dominici & Miratrix (2023), arXiv 2309.06727 | Z (optional) | empirical Bayes double shrinkage |
| Cheng & Cai (2021); Chen, Zhang & Ye (2021), arXiv 2109.10522 | Z / V | adaptive linear combination; minimax rates and adaptivity |
| Yang, Gao, Zeng & Wang (2023), JRSSB, elastic integrative analysis | Z **[new to bib]** | test-then-pool (a pretest, so a selection step) |
| **Oberst et al. (2023), arXiv 2205.10467** | V **[new]** | plug-in MSE-optimal convex weight; they say explicitly that they do not address inference |
| **Yang, Lin, Athey, Jordan & Imbens (2025), arXiv 2511.00727** | V **[new]** | convex weight chosen by cross-validation on the experimental folds; error bounds and bootstrap SEs, no coverage theory in the parts I read |
| **Armstrong, Kline & Sun (2025), Econometrica 93(6)** | V **[new]** | adaptive (soft-thresholding) weighting. They note that, without a bias bound, a CI valid for all biases has length close to the unbiased CI's (citing Armstrong & Kolesár, 2021) |
| **Lin, Bickel & Ding (2026), arXiv 2602.16310** | V **[new]** | the b-value: CIs indexed by a bias bound, i.e. a sensitivity-analysis view of combining estimators |
| van der Laan, Qiu, Tarp & van der Laan (2025), arXiv 2405.07186 | V **[new]** (optional) | adaptive TMLE, a data-adaptive pooled estimator with valid inference |
| Guerre & Wang (2026), arXiv 2609.19415 | V **[new]** (optional) | estimated weights create interaction bias. They note Gu & Chen (2023) use sample splitting to make the weights independent of the estimators, the same device as our split |
| Wang, Schröder, Frauen, Schweisthal, Hess & Feuerriegel (2025), ICLR, arXiv 2412.11511 | V **[new]** | the "Wang et al." paper in your notes: prediction-powered CIs for the ATE from a small unconfounded and a large confounded dataset, with an RCT + OS extension |
| Angelopoulos et al. (2023), Science 382; Demirel et al. (2024) | V / Z | PPI, and its use for generalisation (optional) |
| Kallus, Puli & Shalit (2018) | Z (optional) | bias correction rather than weighting |

Positioning sentence (my reading; check before relying on it): in the
point-identified setting the biased estimator's bias has unknown sign and
size. That is why Armstrong et al. say valid CIs cannot gain without a bias
bound, and why Lin et al. index CIs by one. In our setting both bands are
valid, so the looser bound only ever errs outwards. Every convex combination
of the bounds is then itself a valid, conservative bound (the "Conservative
Oracle" argument in `toy_example.tex`). Width gains with valid coverage are
therefore possible, and the only cost is the data-dependent weight.

More careful version, to use alongside or instead of the line above (added
2026-11-08). It makes explicit that a valid band is itself a bound on the
bias, set through the sensitivity parameter, and it does not compare a CI for
a point with a CI for a band:

> In both settings, validity requires bounding how wrong each source can be.
> Point-identified fusion does this with an explicit bias bound (Armstrong et
> al., 2025; Lin et al., 2026); partial-identification fusion builds it into
> each source's band through its sensitivity parameter. Given valid bands,
> every convex combination of the bounds is itself valid, so the remaining
> difficulty is only the data-dependent weight.

---

## Not cited (and why)

- Semenova, "Adaptive estimation of intersection bounds" (from memory): not
  found on alphaXiv, so I cannot confirm it exists. Drop it unless you know it.
- Search hits I did not read: "Inference conditional on selection: a review"
  (arXiv 2604.09779), "Shrinkage through multiple identifiability"
  (2604.18430), "Prediction-powered data fusion for treatment effect
  estimation" (2610.12332). Each could be relevant; read one before citing it.
- The A/B-testing winner's-curse papers (2026) are not needed.

## Steps (after approval)

1. You add the papers you want to Zotero (zotero-cli is read-only here),
   using the identifiers below, then re-export the `.bib`. Also remove the
   duplicate `horowitz_nonparametric_2000-1`, and fix `yu_using_nodate` to
   `yu_using_2025`. I then check the new keys against the draft.
2. Draft the three files in connected prose: British English, `\citet` and
   `\citep`, Methods notation (omega_0, omega_1, omega^ast), no "first" claims.
3. Add the section to `compiled.tex` (once you agree where), compile, and
   check for undefined citations.
4. Run the citation-checker agent on the new section, to check each claim
   against its source; then proofread.
5. Prepare a condensed version for CLeaR, for you to copy across.

Also: `Results/Non-Causal Toy/toy_example.tex:137` uses a bare `\cite`, which
should be `\citet`. Fix it in the same pass.

## Open questions

1. Placement: a separate "Related Work" section after the abstract, or folded
   into the Introduction once that is written?
2. Length: what is the CLeaR page budget? That decides how many of the
   optional references survive.
3. Keep the econometrics data-combination references (recommended: they are
   what make a "first" claim untenable) and the falsification sentence
   (optional)?

## Papers to add to Zotero

**By arXiv ID**, all checked on 2026-11-08. Zotero's "Add Item by
Identifier" accepts these, and several can be pasted at once. Where a
journal version exists, it is noted; you may prefer to add that instead.

| arXiv ID | Paper |
|---|---|
| 2403.00422 | Han & McCloskey, Inference for interval-identified parameters selected from an estimated set |
| 2305.14265 | Armstrong, Kline & Sun, Adapting to misspecification (journal version: Econometrica 93(6), 2025) |
| 2602.16310 | Lin, Bickel & Ding, Introducing the b-value |
| 2205.10467 | Oberst et al., Understanding the risks and rewards of combining unbiased and possibly biased estimators |
| 2511.00727 | Yang, Lin, Athey, Jordan & Imbens, Cross-validated causal inference |
| 2412.11511 | Wang et al., Constructing confidence intervals for average treatment effects from multiple datasets (ICLR 2025) |
| 2212.09009 | Zrnic & Fithian, Locally simultaneous inference (check the journal version) |
| 2508.12206 | Sakaguchi, The identification power of combining experimental and observational data for distributional treatment effect parameters |
| 2109.10522 | Chen, Zhang & Ye, Minimax rates and adaptivity in combining experimental and observational data |
| 1906.00505 | Benjamini, Hechtlinger & Stark, Confidence intervals for selected parameters (optional) |
| 1410.2597 | Fithian, Sun & Taylor, Optimal inference after model selection (optional) |
| 2407.20386 | Bugni, Gao, Obradović & Velez, On the power properties of inference for parameters with interval identified sets (optional) |
| 2405.07186 | van der Laan, Qiu, Tarp & van der Laan, Adaptive-TMLE (optional) |
| 2609.19415 | Guerre & Wang, Aggregating many estimators using estimated weights (optional) |

**By title**: journal papers whose details I took from verified reference
lists. I don't have their DOIs to hand, so search by title.

- Manski (1990), "Nonparametric bounds on treatment effects", AER P&P 80:319–323.
- Manski & Pepper (2000), "Monotone instrumental variables", Econometrica 68:997–1010.
- Fan, Sherman & Shum (2014), "Identifying treatment effects under data combination", Econometrica 82:811–822.
- Stoye (2009), "More on confidence intervals for partially identified parameters", Econometrica 77:1299–1315.
- Romano, Shaikh & Wolf (2014), "A practical two-step method for testing moment inequalities", Econometrica 82:1979–2002 (optional).
- Andrews, Roth & Pakes (2023), "Inference for linear conditional moment inequalities", REStud 90:2763–2791 (optional).
- McCloskey (2024), "Hybrid confidence intervals for informative uniform asymptotic inference after model selection", Biometrika 111:109–127 (optional).
- Berk, Brown, Buja, Zhang & Zhao (2013), "Valid post-selection inference", AoS 41:802–837.
- Lee, Sun, Sun & Taylor (2016), "Exact post-selection inference, with application to the lasso", AoS 44:907–927.
- Green & Strawderman (1991), "A James–Stein type estimator for combining unbiased and possibly biased estimators", JASA 86:1001–1006.
- Bates & Granger (1969), "The combination of forecasts", Operational Research Quarterly 20:451–468.
- Chernozhukov et al. (2018), "Double/debiased machine learning for treatment and structural parameters", Econometrics Journal 21:C1–C68.
- Angelopoulos, Bates, Fannjiang, Jordan & Zrnic (2023), "Prediction-powered inference", Science 382:669–674 (optional).

**Already in Zotero, not in the `.bib`**: add them to the collection that
exports to `references.bib` (it looks like "Vignesh").

- Levis et al. (2023), covariate-assisted bounds (`3PJVNE24`)
- Rosenman, Basse, Owen & Baiocchi, shrinkage estimators (`P5NDBAUH`, the
  preprint; the journal version is Biometrics 79(4), 2023)
- Yang, Gao, Zeng & Wang (2023), elastic integrative analysis (`2NLUMMGF`)
- Cheng & Cai (2021), adaptive combination (`DEGK6VM2`)
- Rosenman, Dominici & Miratrix (2023), EB double shrinkage (`2TUHNS5Q`, optional)
- Demirel et al. (2024), prediction-powered generalisation (`IGSUS8S5`, optional)
- Kallus, Puli & Shalit (2018), experimental grounding (`GRD3LRJI`, optional)

**From memory (M)**: Hirano & Porter (2012), Fang & Santos (2019),
Andrews (2000), Haile & Tamer (2003), Cross & Manski (2002),
Ridder & Moffitt (2007), Tamer (2010), Molinari (2020),
Zhong & Prentice (2008), Luedtke & van der Laan (2016), Cox (1975). Check
each exists as cited before adding it; I can look any of them up on request.
