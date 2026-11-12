# GROPlant licensing and attribution

Copyright (c) 2026 GROPlant contributors.

GROPlant is the public model/project name. The frozen scientific model identity
remains `strand_seiplant_v1`, with repeat sequences retained and 14 / 5 / 3 SRX
libraries. The repository and Python package name remain
`plant-pausing-bowtie2`; the display name does not imply a newly trained model.

## Which license applies?

| Material | License and scope |
| --- | --- |
| Original code in `model/src/`, `model/scripts/`, `model/tests/`, `figure/scripts/`, `figure/plotting/`, the root installation checker, configuration and associated documentation | [MIT](../LICENSE) |
| The three pretrained GROPlant checkpoint assets identified in `model/metadata/pretrained_models.json` | [CC BY 4.0](../LICENSES/CC-BY-4.0.txt) |
| Original project-produced plotting/evaluation tables, bootstrap results, figures and previews; original example predictions and dataset curation | [CC BY 4.0](../LICENSES/CC-BY-4.0.txt), to the extent the contributors hold applicable rights |
| Third-party reference DNA, original experimental data, annotations, dependencies and incorporated upstream material | Their original terms and notices, described below; this project does not relicense third-party rights |

Both selected licenses permit commercial use. For code, retain the MIT
copyright and license notice. For CC BY material, credit GROPlant contributors,
link to the source and the [CC BY 4.0 license](https://creativecommons.org/licenses/by/4.0/),
and indicate changes, as required by that license. No endorsement is implied.
The license texts govern; these sentences are an explanation, not extra terms.

Suggested attribution for a reused checkpoint, figure or data table:

> GROPlant contributors (2026), GROPlant, model strand_seiplant_v1
> (repeat-retained, 14/5/3 SRX), plant-pausing-bowtie2 repository,
> CC BY 4.0. Identify the exact release and describe any changes.

Repository: https://github.com/wangyiduo217-beep/plant-pausing-bowtie2

`CITATION.cff` supplies software citation metadata. Cite the source studies and
reference assemblies when scientifically appropriate. No manuscript DOI for
GROPlant is claimed.

## Third-party material

- The model was developed with reference to SeiPlant. Its upstream
  `Copyright (c) 2025 lvtongxuan-Bioinfo` MIT notice is retained verbatim in
  [LICENSES/SeiPlant-MIT.txt](../LICENSES/SeiPlant-MIT.txt). Retain that notice
  for incorporated or adapted upstream material.
- HiLearn 0.2.2 is an external Apache-2.0 dependency. Its own license remains
  applicable; GROPlant does not redistribute it under MIT or CC BY.
- The FASTA records in `model/examples/real_dna/*.fa` are extracts of the listed
  Ensembl Plants reference assemblies. They are not newly authored DNA and are
  not relicensed as original GROPlant content. Source URLs, assembly accessions,
  genomic coordinates and source identities are in the example manifest.
  See the [Ensembl data policy](https://www.ensembl.org/info/about/legal/disclaimer.html)
  and any applicable source-specific terms. Ensembl notes that third-party
  conditions can remain applicable.
- Original sequencing experiments and annotations retain their source
  attribution. Sample accessions and assembly provenance are recorded in
  `model/metadata/samples.tsv` and `model_identity.json`.
- Font binaries are not distributed. Use an appropriately licensed installed
  font; Arial is the requested publication style and Liberation Sans is a
  supported validation fallback.

These grants cover rights held by GROPlant contributors. They do not claim
ownership of public-domain facts, third-party sequence data or upstream works.
