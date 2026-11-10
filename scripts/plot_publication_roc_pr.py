"""Reproduce the selected repeat-retained v1 ROC/PR figure with HiLearn 0.2.2."""
from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.metadata
import importlib.util
import io
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt, font_manager as fm
from matplotlib.text import Text
from pypdf import PdfReader, PdfWriter, Transformation
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

ROOT = Path(__file__).resolve().parents[1]
SPECS = [("arabidopsis", "Arabidopsis", 45066, 14),
         ("wheat", "Wheat", 98656, 5), ("maize", "Maize", 49786, 3)]
HILEARN_SHA256 = "d5d591457564af8567eaa8cbed135fd6dea001854070b6b29fd34434bf6ade01"


def official_curves():
    """Load the installed, unchanged official module without unrelated optional imports."""
    if importlib.metadata.version("hilearn") != "0.2.2":
        raise RuntimeError("Install hilearn==0.2.2 to reproduce this selected figure.")
    package = importlib.util.find_spec("hilearn")
    source = Path(next(iter(package.submodule_search_locations))) / "plot" / "curve_plot.py"
    if hashlib.sha256(source.read_bytes()).hexdigest() != HILEARN_SHA256:
        raise RuntimeError("The installed HiLearn curve module differs from the frozen source.")
    spec = importlib.util.spec_from_file_location("publication_hilearn_curves", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_metrics(truth, score):
    """All-score-threshold metrics, distinct from the native HiLearn plot areas."""
    positive = int(truth.sum()); negative = len(truth) - positive
    auc = (rankdata(score)[truth].sum() - positive * (positive + 1) / 2) / (positive * negative)
    order = np.argsort(-score, kind="stable")
    sorted_score = score[order]; sorted_truth = truth[order]
    ends = np.r_[np.flatnonzero(sorted_score[:-1] != sorted_score[1:]), len(score) - 1]
    tp = np.cumsum(sorted_truth, dtype=np.int64)[ends]
    recall = tp / positive; precision = tp / (ends + 1)
    ap = np.dot(np.diff(np.r_[0., recall]), precision)
    return float(auc), float(ap)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "publication/data/roc")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "publication/generated/roc_pr")
    parser.add_argument("--font-family", default="Arial")
    parser.add_argument("--font-file", type=Path, help="Optional locally licensed TTF; never redistributed.")
    args = parser.parse_args()
    output = args.output_dir.resolve(); output.mkdir(parents=True, exist_ok=True)
    if args.font_file:
        fm.fontManager.addfont(str(args.font_file))
        family = fm.FontProperties(fname=str(args.font_file)).get_name()
    else:
        regular_path = fm.findfont(fm.FontProperties(family=args.font_family))
        family = fm.FontProperties(fname=regular_path).get_name()
    regular_path = fm.findfont(fm.FontProperties(family=family))
    bold_path = fm.findfont(fm.FontProperties(family=family, weight="bold"))
    pdfmetrics.registerFont(TTFont("Publication", regular_path))
    pdfmetrics.registerFont(TTFont("Publication-Bold", bold_path))
    plt.rcParams.update({"font.family": family, "font.size": 7, "axes.labelsize": 7,
                         "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
                         "axes.linewidth": .5, "pdf.fonttype": 42, "svg.fonttype": "none",
                         "path.simplify": False})
    curves = official_curves()
    width, height = 183., 75.
    fig = plt.figure(figsize=(width / 25.4, height / 25.4), facecolor="white")
    records = []; hashes = {}
    for index, (key, name, n, srx) in enumerate(SPECS):
        path = args.data_dir / f"{key}.full_input.npz"
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path, allow_pickle=False) as data:
            predicted = data["predicted"]; labels = data["labels_full"]
            overlap = data["overlap_bp"]
            assert predicted.shape == labels.shape == overlap.shape == (n, 2)
            assert np.array_equal(labels, overlap > 0)
            assert np.all(data["end0"] - data["start0"] == 1024)
            truth = labels.ravel(order="C").astype(bool)
            score = predicted.ravel(order="C").astype(np.float64)
            np.testing.assert_array_equal(truth, data["truth"])
            np.testing.assert_array_equal(score, data["score"])
        x = 11 + index * 62
        ax = fig.add_axes([x / width, 21 / height, 47 / width, 42 / height])
        plt.sca(ax)
        state = np.random.get_state()
        try:
            np.random.seed(42)
            fpr, tpr, roc_thresholds, roc_auc = curves.ROC_plot(
                truth, score, threshold=None, color="#E64B35", legend_on=False,
                base_line=True, linewidth=.85)
            roc_line, baseline = ax.lines[-2:]
            np.random.seed(42)
            recall, precision, pr_thresholds, pr_auc = curves.PR_curve(
                truth, score, threshold=None, color="#3C5488", legend_on=False,
                base_line=False, linewidth=.85)
            pr_line = ax.lines[-1]
        finally:
            np.random.set_state(state)
        roc_line.set_label(f"ROC curve (AUROC = {roc_auc:.3f})")
        pr_line.set(linestyle="--", label=f"PR curve (PR AUC = {pr_auc:.3f})")
        baseline.set(color="#909090", linewidth=.5, label="Random ROC")
        ax.set(xlim=(0, 1), ylim=(0, 1.05), xlabel="FPR (ROC) / Recall (PR)",
               ylabel="TPR (ROC) / Precision (PR)")
        ticks = [0, .2, .4, .6, .8, 1]
        ax.set_xticks(ticks); ax.set_yticks(ticks)
        ax.get_xticklabels()[0].set_ha("left")
        ax.get_xticklabels()[-1].set_ha("right")
        ax.get_yticklabels()[0].set_va("bottom")
        ax.tick_params(pad=1.5, width=.5, length=2, labelsize=7)
        ax.grid(True, color="#E3E3E3", linewidth=.3)
        legend = ax.legend(handles=[roc_line, pr_line, baseline], loc="lower right",
                           fontsize=7, framealpha=.94, handlelength=1.3,
                           handletextpad=.35, borderpad=.32, labelspacing=.25)
        legend.get_frame().set_linewidth(.4)
        fig.text((x + 23.5) / width, 71.8 / height, f"{name} HiLearn ROC and PR",
                 ha="center", fontsize=7, fontweight="bold")
        fig.text((x + 23.5) / width, 67.9 / height,
                 f"{n:,} test windows; two strands pooled", ha="center", fontsize=7)
        fig.text((x - 10.5) / width, 71.8 / height, "abc"[index],
                 fontsize=8, fontweight="bold")
        full_auc, full_ap = exact_metrics(truth, score)
        records.append(dict(species=key, n_windows=n, n_observations=len(truth),
                            n_SRX=srx, positive_n=int(truth.sum()), prevalence=float(truth.mean()),
                            HiLearn_AUROC=float(roc_auc), HiLearn_PR_AUC=float(pr_auc),
                            full_exact_AUROC=full_auc, full_exact_AP=full_ap,
                            n_ROC_thresholds=len(roc_thresholds), n_PR_thresholds=len(pr_thresholds)))
        np.savez_compressed(output / f"{key}.native_curves.npz", roc_fpr=fpr, roc_tpr=tpr,
                            roc_thresholds=roc_thresholds, pr_recall=recall, pr_precision=precision,
                            pr_thresholds=pr_thresholds)
    fig.text(0, 5 / height,
             "HiLearn 0.2.2; seed 42; all test records retained. Native PR AUC is not AP.", fontsize=7)
    fig.canvas.draw()
    for text in fig.findobj(match=Text):
        if text.get_visible() and text.get_text():
            assert text.get_fontsize() in [7., 8.]
    stem = output / "three_species_ROC_PR_183mm"
    for extension in ["pdf", "svg", "png"]:
        fig.savefig(stem.with_suffix("." + extension), dpi=600, facecolor="white")
    plt.close(fig)
    with (output / "metrics.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(records)
    mm = 72 / 25.4
    buffer = io.BytesIO(); page = canvas.Canvas(buffer, pagesize=(210 * mm, 297 * mm))
    page.setFont("Publication-Bold", 7)
    page.drawString(13.5 * mm, 278 * mm, "Repeat-retained strand_seiplant_v1: 14 / 5 / 3 SRX")
    page.setFont("Publication", 7)
    caption = [
        "All formal held-out test windows are included; the two strands are pooled in matching order.",
        "Positive truth means at least one base overlaps a same-strand experimental consensus interval.",
        "The native, unmodified HiLearn 0.2.2 functions are called, with seed 42 reset before each call.",
        "HiLearn samples up to 2,000 score thresholds, never test observations. No score filter is used.",
        "Its PR curve includes synthetic display endpoints; the reported PR AUC is trapezoidal, not AP.",
        "Full-threshold exact AUROC and AP are stored separately in metrics.tsv.",
        "Grey diagonal: ROC random baseline. Native PR diagonal is disabled; PR chance is prevalence.",
        "Ordinary text: 7 pt; panel letters: 8 pt. The 183 x 75 mm figure is placed without scaling.",
        f"Actual installed font family: {family}. No proprietary font files are distributed."]
    for row, sentence in enumerate(caption):
        page.drawString(13.5 * mm, (175 - row * 5) * mm, sentence)
    page.showPage(); page.save(); buffer.seek(0)
    a4 = PdfReader(buffer).pages[0]
    a4.merge_transformed_page(PdfReader(stem.with_suffix(".pdf")).pages[0],
                              Transformation().translate(tx=13.5 * mm, ty=192 * mm))
    writer = PdfWriter(); writer.add_page(a4)
    with (output / "three_species_ROC_PR_A4.pdf").open("wb") as handle:
        writer.write(handle)
    (output / "reproduction.json").write_text(json.dumps({
        "model_version": "strand_seiplant_v1", "repeat_sequences_retained": True,
        "n_SRX": [14, 5, 3], "font_family": family, "text_pt": 7, "panel_pt": 8,
        "hilearn_version": "0.2.2", "official_curve_source_sha256": HILEARN_SHA256,
        "input_sha256": hashes, "records": records,
        "score_filter": None, "HiLearn_PR_AUC_is_not_AP": True}, indent=2) + "\n", encoding="utf-8")
    print(f"Created ROC/PR figure and metrics in {output}")


if __name__ == "__main__":
    main()
