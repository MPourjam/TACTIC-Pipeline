import sys
from pathlib import Path

import pytest
from ete3 import Tree


# analysis.py imports processing_helper as a top-level module.
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tactic_pipeline"))
from tactic_pipeline.analysis import analysis  # noqa: E402


@pytest.mark.parametrize(
    ("method", "suffix"),
    [
        (["rapidnj"], "-Tree-nj.tre"),
        (["FastTree"], "-Tree-FastTree.tre"),
    ],
)
def test_plant_tree_removes_redundant_quotes_without_changing_lengths(
        tmp_path, monkeypatch, method, suffix):
    alignment = tmp_path / "aligned.fasta"
    alignment.write_text(">Zotu1\nACGT\n>Zotu2\nACGT\n", encoding="utf-8")
    monkeypatch.chdir(tmp_path)

    def fake_system_sub(command):
        Path(command[-1]).write_text(
            "('Zotu1':0.123,'OTU2':0.456,'sample one':0.7);\n",
            encoding="utf-8",
        )

    monkeypatch.setattr(analysis, "system_sub", fake_system_sub)
    analysis.plant_tree(alignment, method=method, output_file_name_root="IDs")

    tree_text = (tmp_path / f"IDs{suffix}").read_text(encoding="utf-8")
    assert tree_text == "(Zotu1:0.123,OTU2:0.456,'sample one':0.7);\n"


@pytest.mark.parametrize("quoted", [True, False])
def test_extract_and_rename_subtree_accepts_old_and_new_labels(
        tmp_path, quoted):
    source = tmp_path / "ZOTUs-Tree-nj.tre"
    label = lambda name: f"'{name}'" if quoted else name
    source.write_text(
        f"(({label('Zotu1')}:0.1,{label('Zotu2')}:0.2):0.3,"
        f"{label('Zotu3')}:0.4);\n",
        encoding="utf-8",
    )
    output = tmp_path / "SOTUs-Tree-nj-TIC.tre"

    analysis.extract_and_rename_subtree(
        source,
        subset=["Zotu1", "Zotu3"],
        rename_map={"Zotu1": "SOTU1", "Zotu3": "SOTU3"},
        output_file=output,
    )

    assert set(Tree(str(output)).get_leaf_names()) == {"SOTU1", "SOTU3"}
    assert "'" not in output.read_text(encoding="utf-8")
