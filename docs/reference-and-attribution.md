# References, provenance, and licensing status

## Selected project snapshot

The documented research snapshot is `strand_seiplant_v1`: repeat sequences are
retained, and the selected experimental inputs are 14 Arabidopsis, 5 wheat, and
3 maize SRX libraries. Species-specific models take a 1,024-bp DNA sequence and
produce two continuous GRO-seq scores, one for each reference strand. These are
project-specific GRO-seq models, rather than the histone-mark classifiers or
pretrained weights distributed by the reference study.

## SeiPlant reference and implementation

Lv T, Han Q, Li Y, Liang C, Ruan Z, Chao H, Chen M, Chen D. **Cross-species
prediction of histone modifications in plants via deep learning.** *Genome
Biology* **27**, 20 (2026).
[Article and DOI](https://doi.org/10.1186/s13059-025-03929-4).

The article links to the public
[compbioNJU/SeiPlant repository](https://github.com/compbioNJU/SeiPlant), a fork of
[Lv-BioInfo/SeiPlant](https://github.com/Lv-BioInfo/SeiPlant). The implementation
reviewed for this project is commit
[`1c88cfcfc67c37c4b9f5deb56175a09b8ffdd853`](https://github.com/compbioNJU/SeiPlant/tree/1c88cfcfc67c37c4b9f5deb56175a09b8ffdd853).
It informed the convolutional and residual backbone, dilated convolutions,
fixed cubic B-spline representation, and genome-scanning workflow. The
project-specific target construction, two-output regression head, and GRO-seq
evaluation are documented separately in this repository.

We also follow the reference repository's presentation pattern: installation
and quick-start instructions, separate model and execution code, a from-scratch
training guide, and a clearly identified plotting workflow. This organizational
reference does not identify our training labels, model weights, or GRO-seq
measurements as outputs of the reference paper.

The reference repository declares the
[MIT License](https://github.com/compbioNJU/SeiPlant/blob/1c88cfcfc67c37c4b9f5deb56175a09b8ffdd853/LICENSE),
with the notice `Copyright (c) 2025 lvtongxuan-Bioinfo`. That upstream attribution
and the upstream license apply to any SeiPlant material reused under those
terms. They do not select a license for newly authored project code.

The reference study's separate
[Zenodo archive](https://doi.org/10.5281/zenodo.15421964) contains its example
Arabidopsis reference and histone-model weights. Those histone weights are not
our GRO-seq checkpoints. Citing that archive does not claim that it contains
this project's trained weights or plotting data.

## HiLearn plotting dependency

The ROC/PR plotting workflow uses the separately installed
[`hilearn==0.2.2`](https://pypi.org/project/hilearn/0.2.2/) package, with its
published `ROC_plot` and `PR_curve` functions. This version is pinned to retain
the plotting and area-calculation behavior used for the selected result.

HiLearn is maintained by Yuanhua Huang and declares the
[Apache License 2.0](https://github.com/huangyh09/hilearn/blob/master/LICENSE).
Its source distribution is obtained through the package dependency; the
HiLearn package or its source archive is not vendored in this repository.
The [official plotting documentation](https://hilearn.readthedocs.io/en/latest/plotting.html)
describes its interfaces. The dependency's license remains its own and is not a
project-wide license declaration for this repository.

## GROPlant licensing

The original project code and associated documentation are licensed under MIT
(Copyright (c) 2026 GROPlant contributors). The project's pretrained weights
and original generated figures/data are licensed under CC BY 4.0. The complete
scope and attribution instructions are in [licensing.md](licensing.md).

Third-party reference DNA, original experiments and dependencies retain their
own terms. The upstream SeiPlant MIT notice is preserved in
`LICENSES/SeiPlant-MIT.txt`; the external HiLearn package remains Apache-2.0.
