# DNAPL CSAT Calculator — Version 2

Requires Python 3.10 or newer. Extract the folder, activate your environment, and run:

```bash
python -m pip install -r requirements.txt
python -m streamlit run CSAT_Calculator_V2.py
```

The Python filename remains unchanged for existing launch commands.

## Changes in Version 2

Calculation 4 now displays the screenshot inequality directly: `sum(CTobs,i / CTS,i) >= 1`. The unknown-composition soil result column is labeled **Saturation Index** in the UI and all exported tables. Calculation values and screening rules are unchanged.

- Solubility terminology is used consistently in controls, result tables, documentation and Excel exports.
- The known-composition soil table labels are **Saturation Index** and **Soil Concentration Threshold [CiT] (mg/kg)**.
- The unknown-composition soil table uses **Single-component threshold [CTS,i] (mg/kg)** and now includes a **Threshold met or exceeded** flag for each component. The flag compares its observed concentration against CTS,i; the separate sample CSAT sum still adds component indices.
- The input schema now requires **matrix**. Soil and groundwater calculations are separate, and matrix is retained in full sample metadata and exports.
- A new **1% Effective Solubility Threshold Not Knowing DNAPL Composition** card implements Calculation 6. It is always visible directly below the unknown-composition soil card. It is greyed out with its filter and export disabled until No is selected, groundwater records are uploaded, and valid Calculation 6 results are available. Inactive tables are empty. All three calculation boxes remain visible.
- The groundwater table includes component rows and a sum row for each sample, plus a sample-level possible-DNAPL flag. The Excel export includes that table, a sample-summary sheet, analyte properties and Calculation 6 settings.

## Example inputs

Select two components, PCE and TCE, for any of the included examples:

- `example_input.csv`: soil only, with the required matrix column. Use Yes for composition-adjusted thresholds or No for soil CSAT.
- `example_groundwater_input.csv`: groundwater only. Select No; one sample exceeds the 1% rule and one does not. Soil density and moisture rows are unnecessary for groundwater.
- `example_mixed_input.csv`: soil and groundwater together. With No selected, both soil CSAT and groundwater Calculation 6 results are available.

## Input schema and validation

Required columns:

`location, field_id, top_depth, bottom_depth, depth_units, lab_id, cas_no, analyte, result, units, matrix`

Matrix accepts **Soil**, **Ground Water**, or **Groundwater**, case-insensitively. Leading/trailing whitespace is ignored, and groundwater spellings are normalized to Groundwater. Other matrix values and a missing matrix column stop processing. Existing input files need a matrix column added before use.

Use one row per chemical or soil property per sample. Populate the full metadata key on every row: location, field_id, lab_id, top_depth, bottom_depth, depth_units, matrix. Depths must be finite numbers, with bottom_depth >= top_depth. Groundwater samples also retain the depth metadata, so samples from different intervals are not combined. Results from different matrices, sample IDs, lab IDs or depths are never pooled, even if the location name is the same.

Chemical matching uses **cas_no only**, regardless of the laboratory analyte name. CAS accepts whitespace, Unicode hyphens and digit-only forms. Supported reference chemicals remain PCE, TCE, cis-1,2-DCE, vinyl chloride, 1,1,1-TCA and 1,4-dioxane, with one to five unique components selected for a run. Unrecognized records are excluded with a warning. Supported but unselected chemicals are also excluded with a warning.

Units:

| Record | Required units |
|---|---|
| Soil chemical concentration | mg/kg |
| Groundwater concentration | mg/L |
| Soil moisture | % |
| Dry soil bulk density | g/cm3; g/cm³, g/cc and g/mL equivalents also accepted |

Units are case-insensitive and tolerate normal spacing. All groundwater concentration records must be in mg/L, including unrecognized chemicals. Other units, including ug/L, must be converted before uploading; no automatic unit conversion is performed. Soil moisture and density rows must use matrix Soil.

Each soil sample requires exactly one moisture result, one dry bulk density result and one result for each selected chemical. Groundwater samples need only chemical concentration rows: no moisture, density or Henry's constant is used in Calculation 6. For groundwater, only selected measured components are included; missing selected components generate a warning and are not treated as zero. At least one selected chemical must be measured per groundwater sample. If an unmeasured component may be in the DNAPL, alpha can be underestimated.

Duplicate sample/analyte results stop processing; no averaging occurs. Relevant concentrations must be finite, nonnegative numeric values. ND, `<` results and text qualifiers are rejected. Because there is no detection qualifier field, verify numeric results are detected concentrations rather than reporting limits. Soil concentrations must be on a dry-soil basis.

Soil parameter matching uses case-insensitive token/prefix patterns: `dry bulk density`, `dry bulk den`, `dry bulk dens.`, `bulk density, dry`, `dry density`, `% moisture`, `perc moisture`, `percent moisture`, `moisture content` and the prior `moisuter` typo. Dry density requires a dry designation. Wet bulk density, generic bulk density, particle density and moisture-holding capacity are not substituted. Soil-property CAS may be blank or a nonnumeric lab parameter code; an unknown numeric CAS is not guessed from the name.

## Known-composition soil workflow

Calculations 5 → 3 → 2 convert mass percentages into mole fractions, calculate effective solubility, and derive component soil thresholds:

`mi = (msi/MWi) / sum(msj/MWj)`

`Ci = mi*Si`

`CiT = Ci*H'i*theta_a/rho_b`

The **Saturation Index** is observed soil concentration / CiT. Below 1 is below threshold, 1 is at threshold, and above 1 is above threshold. Each component is assessed independently; these indices are not summed into Calculation 4 CSAT. A zero component fraction gives a zero threshold: its index is blank and its comparison is marked not assessable.

A composition below 100% includes the unknown balance using an editable estimated molecular weight (default: mean of selected molecular weights). Known components are not renormalized to 100%. Results depend on the entered composition and unknown-balance estimate. Effective Solubility is a subsection inside the known-composition soil box.

## Unknown-composition soil workflow

Calculation 4 uses single-component solubility:

`CTS,i = Si*H'i*theta_a/rho_b`

`sum(CTobs,i / CTS,i) >= 1`

The component Saturation Index is the observed concentration divided by its single-component threshold. The sample CSAT total is the sum of the component indices.

The component **Threshold met or exceeded** flag is Yes when its measured concentration >= CTS,i. The sample CSAT sum can reach 1 even when no individual component meets its threshold. Composition percentages are disabled and ignored in this workflow.

Soil assumptions apply to both soil workflows: sorption contribution Kd*rho_b = 0 and additive water contribution theta_w = 0. Measured moisture still reduces air-filled porosity:

`theta_a = 1 - rho_b/rho_s - w*rho_b/rho_w`

rho_b is dry bulk density, rho_s is particle density, rho_w is water density and w is dry-mass gravimetric moisture fraction. Wet-mass moisture is converted using `w_dry = w_wet/(1-w_wet)`. Dry bulk density is not corrected a second time. g/cm3 equals kg/L numerically; no additional factor of 1,000 is needed. Zero or negative air-filled porosity blocks soil calculation because the vapor-only threshold vanishes.

A Yes is a screening exceedance under these assumptions, not confirmation of DNAPL presence. Excluding sorbed and dissolved mass can create an exceedance without DNAPL. The soil-sum export retains the requested labels: CSAT >= 1 gives Yes; below 1 gives `No – attributed to dissolve phase`. A No does not establish absence.

## Groundwater Calculation 6

`alpha = sum(Cobs,i / Si)`

Cobs,i is the sampled groundwater concentration of component i (mg/L); Si is its single-component solubility (mg/L); alpha is the dimensionless cumulative mole-fraction estimate for the sample set; n is the number of selected measured components included in that groundwater sample.

The rule flags **DNAPL may be present = Yes only if alpha > 0.01**. Exactly 0.01 gives No. Comparison uses unrounded alpha, separately for each full sample key. The flag is repeated on component and sum rows for that sample; it does not test each individual ratio against 0.01. Component concentration and solubility cells are blank on the sum row, matching the supplied Calculation 6 example layout.

For example, ratios 0.004, 0.009, 0.007, 0.008 and 0.006 sum to 0.034: the sample exceeds the 1% rule even though no single ratio is greater than 0.01.

The calculation assesses possible DNAPL near the monitoring location. It assumes borehole dilution, dispersion and degradation affect all components equally. Samples are evaluated independently, without assuming spatially uniform DNAPL composition. Unmeasured components and aqueous-phase degradation can cause underestimation or false negatives. A non-exceedance means the rule was not exceeded; it does not demonstrate DNAPL absence.

## Properties, exports and checks

Property defaults remain as supplied in the original CVOC properties screenshot. Molecular weight and Solubility display two decimals; Henry's constant retains nine decimals. Koc is not shown. Formatting does not round intermediate calculations. 1,4-dioxane requires a technically justified numeric solubility because the supplied reference lists it as miscible.

Display filters do not restrict exports. Soil tables support field ID and analyte filters; groundwater supports field ID filtering, preserving all component and sum rows for each displayed sample. Exports retain full sample metadata, including matrix. Groundwater component tables include sum rows, and the separate summary sheet has one alpha and flag per sample. Exported text is not interpreted as Excel formulas.

Run the 37 calculation and validation tests:

```bash
python -m unittest discover -v
```

Run the additional Streamlit UI checks:

```bash
python check_ui_states.py
```

The tests cover soil calculations, Calculation 6 arithmetic from the screenshot, the strict 0.01 boundary, mixed matrices, CAS/name patterns, units, duplicate results, missing measured components and Excel exports. UI checks cover activation and disabling of the always-visible groundwater card, soil-only/groundwater-only/mixed files, workflow switching and invalid input.

Source: Calculations 2–6 in Kueper, B. H., & Davies, K. L. (2009). *Assessment and Delineation of DNAPL Source Zones at Hazardous Waste Sites* (EPA/600/R-09/119). U.S. Environmental Protection Agency, Office of Research and Development, National Risk Management Research Laboratory.
https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1006Y98.TXT
