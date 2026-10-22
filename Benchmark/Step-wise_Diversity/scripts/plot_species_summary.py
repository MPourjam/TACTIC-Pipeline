"""Render the archived top-species measurements without rerunning clustering."""
import csv
from pathlib import Path


def plot(summary, output, top_n=25, vertical=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    with Path(summary).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        rows = [r for r in reader if not r["taxonomy_label"].startswith("__")][:top_n]
    cluster_fields = [f for f in fields if f.startswith("Number of Clusters at ")]
    series = ["unique_v1v9_sequences", *cluster_fields]
    labels = []
    for row in rows:
        taxon = row["taxonomy_label"].split(";")[-1]
        labels.append(taxon[:1].upper() + taxon[1:])
    colors = ["#CC79A7", "#56B4E9", "#009E73", "#E69F00"]
    positions = np.arange(len(rows))
    width = 0.8 / len(series)
    fig, ax = plt.subplots(figsize=(max(8, len(rows) * .55), 6) if vertical else (10, max(5, len(rows) * .42)))
    for index, field in enumerate(series):
        values = [int(row[field].replace(",", "")) for row in rows]
        offset = (index - (len(series) - 1) / 2) * width
        legend = "Unique V1-V9 sequences" if index == 0 else field.replace("Number of Clusters at", "Clusters at")
        if vertical:
            ax.bar(positions + offset, values, width=width, color=colors[index], label=legend)
        else:
            ax.barh(positions + offset, values, height=width, color=colors[index], label=legend)
    if vertical:
        ax.set_xticks(positions, labels, rotation=45, ha="right", fontstyle="italic")
        ax.set_ylabel("Feature count")
    else:
        ax.set_yticks(positions, labels, fontstyle="italic")
        ax.invert_yaxis()
        ax.set_xlabel("Feature count")
    ax.set_title(f"Top {len(rows)} species: unique sequences and clusters")
    ax.legend(frameon=False)
    fig.tight_layout()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".svg", ".png"):
        fig.savefig(output.with_suffix(suffix), bbox_inches="tight", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--top-n", type=int, default=25)
    parser.add_argument("--vertical", action="store_true")
    args = parser.parse_args()
    plot(args.summary, args.output, args.top_n, args.vertical)
