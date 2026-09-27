# Phenotypic insufficiency of local protein representations

**Fibre-preserving rearrangements of human C2H2 zinc-finger arrays reorder DNA-recognition residues under an invariant 25-spectrum**

Zhengyi Chen<sup>1</sup> and Ruqing Chen<sup>2</sup>

<sup>1</sup> Guangxi Key Laboratory of Drug Discovery and Optimization, Guangxi Engineering Research Center for Pharmaceutical Molecular Screening and Druggability Evaluation, School of Pharmacy, Guilin Medical University, Guilin 541199, P. R. China (chenzhengyi@glmc.edu.cn)
<sup>2</sup> GUT Geoservice Inc., Montreal, Québec, Canada (ruqing@hotmail.com)

DOI: [10.5281/zenodo.22998488](https://doi.org/10.5281/zenodo.22998488)

## About

A protein representation built from local sequence statistics is constant on the *spectral fibre* of a sequence: all sequences with the same multiset of k-mers. For C2H2 zinc-finger proteins, whose fingers read DNA in order, this blind spot can matter for the binding phenotype.

The computations in this repository cover all 17 reviewed human C2H2 zinc-finger proteins with a non-trivial 25-fibre. For each, fibre-mates were drawn as random Eulerian trails of the de Bruijn graph of order 25. Fingers were located by the motif C-X(2-4)-C-X(12)-H-X(3-5)-H, and each finger was summarized by its recognition residues at helix positions −1, 2, 3 and 6.

- Of 1,178 distinct fibre-mates, 883 (75%) place different recognition residues at some array position.
- In all 1,178, the multiset of fingers is the same as in the wild type, and a window-25 convolutional encoder gives identical output.
- In ZNF91 (31 fingers), one fibre-mate changes the recognition residues at 17 array positions. In ZNF721 (23 fingers), one changes 17.
- KLF18 has a large fibre, but none of its 285 mates reorders its three fingers.

Binding was not measured. Statements about binding sites refer to finger-additive recognition models.

Companion papers: [10.5281/zenodo.22950095](https://doi.org/10.5281/zenodo.22950095), [10.5281/zenodo.22959211](https://doi.org/10.5281/zenodo.22959211), [10.5281/zenodo.22960439](https://doi.org/10.5281/zenodo.22960439), [10.5281/zenodo.22963882](https://doi.org/10.5281/zenodo.22963882), [10.5281/zenodo.22966020](https://doi.org/10.5281/zenodo.22966020).

## Repository layout

```
zinc-finger-order/
├── paper/      zinc_finger_order.tex, zinc_finger_order.pdf    the paper
├── code/       finger_order.py       fibre-mates, finger parsing, recognition residues, encoder check
│               make_figures.py       finger-order ladders for ZNF91 and ZNF721
│               spectral_fibres.py    UniProt loading and random Eulerian trails (from the companion repository)
├── data/       human_reviewed.tsv.gz   UniProt snapshot (reviewed, organism 9606)
├── results/    finger_order_output.txt   transcript of finger_order.py
│               finger_order.json         per-protein summary
│               finger_pairs.tsv          every fibre-mate: residues changed, array positions changed, recognition residues
│               case_ZNF91.json, case_ZNF721.json
│               fibre_results.json        N_k and κ per protein (from the companion paper)
├── figures/    fig_ladder.pdf
├── README.md, CITATION.cff, requirements.txt
└── LICENSE (code, MIT), LICENSE-CC-BY-4.0.md (text and figures)
```

## Reproducing

```
cd code
python -u finger_order.py     # about 4 minutes; needs torch for the encoder check
python make_figures.py        # needs matplotlib
cd ../paper
pdflatex zinc_finger_order
pdflatex zinc_finger_order
```

## How to cite

```bibtex
@misc{ChenChen2026ZincFingerOrder,
  author    = {Chen, Zhengyi and Chen, Ruqing},
  title     = {Phenotypic insufficiency of local protein representations: fibre-preserving
               rearrangements of human C2H2 zinc-finger arrays reorder DNA-recognition residues
               under an invariant 25-spectrum},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.22998488},
  url       = {https://doi.org/10.5281/zenodo.22998488}
}
```

## License

Text and figures: [CC BY 4.0](LICENSE-CC-BY-4.0.md). Code: [MIT](LICENSE). The UniProt data are distributed by UniProt under CC BY 4.0.
