import html
import re
from pathlib import Path
import pytest
from tactic_pipeline.analysis.analysis_helper import (
    get_pipeline_readme_content,
    write_pipeline_readme,
    write_all_pipeline_readmes,
    prepare_picrust2_input,
    remove_intermediate_files,
)
from tactic_pipeline.analysis.readme_html import parse_readme_text, render_readme_html

PIPELINE_TYPES = ["root", "otu", "zotu", "tac", "tic", "picrust2", "shared"]


def test_get_pipeline_readme_content_otu():
    content = get_pipeline_readme_content("otu")
    assert "OTU Pipeline Output Files" in content
    assert "OTUs-Table.tab" in content
    assert "OTUs-Seqs.fasta" in content
    assert "OTUs-Tree-nj.tre" in content
    assert "OTUs-Krona.html" in content
    assert "SampleMinCount_Normalized-OTUs-Table.tab" in content
    assert "Spike_Normalized-OTUs-Table.tab" in content
    assert "*-8Levels.*" in content
    assert "What it is:" in content
    assert "What it is used for:" in content


def test_get_pipeline_readme_content_zotu():
    content = get_pipeline_readme_content("zotu")
    assert "ZOTU Pipeline Output Files" in content
    assert "ZOTUs-Table.tab" in content
    assert "ZOTUs-Seqs.fasta" in content
    assert "ZOTUs-Tree-nj.tre" in content
    assert "ZOTUs-Krona.html" in content
    assert "SampleMinCount_Normalized-ZOTUs-Table.tab" in content
    assert "Spike_Normalized-ZOTUs-Table.tab" in content
    assert "*-8Levels.*" in content


def test_get_pipeline_readme_content_tac():
    content = get_pipeline_readme_content("tac")
    assert "TAC Pipeline Output Files" in content
    assert "OTUs-Table-TAC.tab" in content
    assert "OTUs-Seqs-TAC.fasta" in content
    assert "Map-ZOTUs-OTUs-TAC.tab" in content
    assert "OTUs-Tree-nj-TAC.tre" in content
    assert "OTUs-Krona.html" in content


def test_get_pipeline_readme_content_tic():
    content = get_pipeline_readme_content("tic")
    assert "TIC Pipeline Output Files" in content
    assert "SOTUs-Table-TIC.tab" in content
    assert "SOTUs-Seqs-TIC.fasta" in content
    assert "Map-SOTU-ZOTU.tab" in content
    assert "Map-GOTU-SOTU.tab" in content
    assert "Map-FOTU-GOTU.tab" in content
    assert "Invalid-Tax-Seqs-TIC.fasta" in content
    assert "ZOTU-Table-OnlyBacteria-FullTaxonomy.tab" in content
    assert "ZOTU-Table-All-FullTaxonomy.tab" in content


def test_get_pipeline_readme_content_picrust2():
    content = get_pipeline_readme_content("picrust2")
    assert "PICRUSt2 Input Files" in content
    assert "picrust2_pipeline.py" in content
    assert "-For-PICRUSt2.tab" in content
    assert "-For-PICRUSt2.fasta" in content


def test_get_pipeline_readme_content_shared():
    content = get_pipeline_readme_content("shared")
    assert "Shared Intermediate Sequence Files" in content
    assert "sorted.fasta" in content
    assert "derep.fasta" in content
    assert "relabel.tab" in content


def test_get_pipeline_readme_content_root():
    content = get_pipeline_readme_content("root")
    assert "Analysis Results Directory Overview" in content
    assert "Output files explanation.html" in content
    assert "Output files explanation.txt" not in content
    assert "OTU_Pipeline/" in content
    assert "ZOTU_Pipeline/" in content
    assert "TAC_Pipeline/" in content
    assert "TIC_Pipeline/" in content
    assert "PICRUSt2-Inputs/" in content
    assert "Pipeline_log.txt" in content
    assert "Analysis_log.txt" in content
    assert "spike_mapping_file.csv" in content


def test_write_pipeline_readme(tmp_path):
    readme_path_str = write_pipeline_readme(str(tmp_path), "otu")
    readme_path = Path(readme_path_str)
    assert readme_path.is_file()
    assert readme_path.name == "Output files explanation.html"
    assert not (tmp_path / "Output files explanation.txt").exists()
    page = readme_path.read_text(encoding="utf-8")
    assert "<title>OTU Pipeline Output Files</title>" in page


def test_write_all_pipeline_readmes_with_active_pipelines(tmp_path):
    active = ["otu", "tac"]
    created = write_all_pipeline_readmes(str(tmp_path), active_pipelines=active)
    assert len(created) == 3
    assert all(Path(p).is_file() for p in created)

    assert (tmp_path / "Output files explanation.html").is_file()
    assert (tmp_path / "OTU_Pipeline" / "Output files explanation.html").is_file()
    assert (tmp_path / "TAC_Pipeline" / "Output files explanation.html").is_file()
    # Inactive and empty should not have Output files explanation.html
    assert not (tmp_path / "TIC_Pipeline" / "Output files explanation.html").exists()

    otu_text = (tmp_path / "OTU_Pipeline" / "Output files explanation.html").read_text(
        encoding="utf-8"
    )
    assert "OTU Pipeline Output Files" in otu_text
    assert "OTUs-Table.tab" in otu_text

    tac_text = (tmp_path / "TAC_Pipeline" / "Output files explanation.html").read_text(
        encoding="utf-8"
    )
    assert "TAC Pipeline Output Files" in tac_text
    assert "OTUs-Table-TAC.tab" in tac_text


def test_write_all_pipeline_readmes_auto_discovers_existing_dirs_with_files(tmp_path):
    tic_dir = tmp_path / "TIC_Pipeline"
    tic_dir.mkdir(parents=True)
    (tic_dir / "SOTUs-Table-TIC.tab").write_text("SOTU\tsample1\ttax\n", encoding="utf-8")

    shared_dir = tmp_path / "Shared"
    shared_dir.mkdir(parents=True)
    (shared_dir / "sorted.fasta").write_text(">seq1\nACGT\n", encoding="utf-8")

    # Empty dir with no files should be skipped
    empty_zotu_dir = tmp_path / "ZOTU_Pipeline"
    empty_zotu_dir.mkdir(parents=True)

    created = write_all_pipeline_readmes(str(tmp_path))
    assert len(created) == 3
    assert all(Path(p).is_file() for p in created)

    assert (tmp_path / "Output files explanation.html").is_file()
    assert (tic_dir / "Output files explanation.html").is_file()
    assert (shared_dir / "Output files explanation.html").is_file()
    assert not (empty_zotu_dir / "Output files explanation.html").exists()

    tic_text = (tic_dir / "Output files explanation.html").read_text(encoding="utf-8")
    assert "TIC Pipeline Output Files" in tic_text


def test_prepare_picrust2_input_creates_picrust2_readme(tmp_path):
    zotu_dir = tmp_path / "ZOTU_Pipeline"
    zotu_dir.mkdir(parents=True)
    (zotu_dir / "ZOTUs-Table.tab").write_text(
        "#OTU ID\tsample1\ttaxonomy\nZOTU_1\t10\tBacteria;\n",
        encoding="utf-8"
    )
    (zotu_dir / "ZOTUs-Seqs.fasta").write_text(
        ">ZOTU_1 tax=Bacteria;\nACGTACGT\n",
        encoding="utf-8"
    )

    prepare_picrust2_input(str(tmp_path))

    picrust_dir = tmp_path / "PICRUSt2-Inputs"
    assert picrust_dir.is_dir()
    readme_file = picrust_dir / "Output files explanation.html"
    assert readme_file.is_file()
    assert not (picrust_dir / "Output files explanation.txt").exists()
    text = readme_file.read_text(encoding="utf-8")
    assert "PICRUSt2 Input Files" in text
    assert "picrust2_pipeline.py" in text


@pytest.mark.parametrize("pipeline_type", PIPELINE_TYPES)
def test_html_readme_documents_every_file_of_text_readme(pipeline_type):
    content = get_pipeline_readme_content(pipeline_type)
    page = render_readme_html(content)
    file_names = [line[2:].split(" (")[0].strip() for line in content.splitlines() if line.startswith("• ")]
    assert file_names
    for name in file_names:
        assert html.escape(name, quote=False) in page, f"{name} missing in HTML of {pipeline_type}"
    title = content.splitlines()[1].split(" — ")[-1]
    assert f"<title>{html.escape(title)}</title>" in page


@pytest.mark.parametrize("pipeline_type", PIPELINE_TYPES)
def test_parse_readme_text_keeps_every_description(pipeline_type):
    content = get_pipeline_readme_content(pipeline_type)
    doc = parse_readme_text(content)
    descriptions = [
        field["text"]
        for section in doc["sections"] for block in section["blocks"] if block["kind"] == "files"
        for field in block["fields"]
    ]
    # Every "- ..." description line of the text version ends up in a parsed field
    expected = [line for line in content.splitlines() if re.match(r"^  - ", line)]
    assert len(descriptions) == len(expected)
    assert all(text for text in descriptions)


def test_write_pipeline_readme_writes_html_with_links_to_existing_files(tmp_path):
    (tmp_path / "ZOTUs-Table.tab").write_text("#OTU ID\ts1\ttaxonomy\n", encoding="utf-8")
    write_pipeline_readme(str(tmp_path), "zotu")
    html_path = tmp_path / "Output files explanation.html"
    assert html_path.is_file()
    page = html_path.read_text(encoding="utf-8")
    assert page.startswith("<!DOCTYPE html>")
    assert '<a class="open" href="ZOTUs-Table.tab">Open</a>' in page
    # Files not in the folder are documented but not linked
    assert "ZOTUs-Tree-nj.tre" in page
    assert 'href="ZOTUs-Tree-nj.tre"' not in page


def test_write_pipeline_readme_falls_back_to_text(tmp_path, monkeypatch):
    def _failing_render(*args, **kwargs):
        raise RuntimeError("rendering failed")

    monkeypatch.setattr("tactic_pipeline.analysis.analysis_helper.render_readme_html", _failing_render)
    written = Path(write_pipeline_readme(str(tmp_path), "otu"))
    assert written.name == "Output files explanation.txt"
    assert "TACTIC Pipeline — OTU Pipeline Output Files" in written.read_text(encoding="utf-8")
    assert not (tmp_path / "Output files explanation.html").exists()


def test_write_pipeline_readme_replaces_outdated_text_fallback(tmp_path):
    (tmp_path / "Output files explanation.txt").write_text("old", encoding="utf-8")
    written = Path(write_pipeline_readme(str(tmp_path), "otu"))
    assert written.name == "Output files explanation.html"
    assert not (tmp_path / "Output files explanation.txt").exists()


def test_root_html_readme_links_to_pipeline_pages(tmp_path):
    write_all_pipeline_readmes(str(tmp_path), active_pipelines=["zotu"])
    root_page = (tmp_path / "Output files explanation.html").read_text(encoding="utf-8")
    assert (tmp_path / "ZOTU_Pipeline" / "Output files explanation.html").is_file()
    assert 'href="ZOTU_Pipeline/Output%20files%20explanation.html"' in root_page
    # Inactive pipeline folders without content are not linked
    assert 'href="TIC_Pipeline/' not in root_page


def test_picrust2_html_readme_lists_matching_files(tmp_path):
    picrust_dir = tmp_path / "PICRUSt2-Inputs"
    picrust_dir.mkdir()
    (picrust_dir / "ZOTU_Pipeline__ZOTUs-Table-For-PICRUSt2.tab").write_text("x", encoding="utf-8")
    write_pipeline_readme(str(picrust_dir), "picrust2")
    page = (picrust_dir / "Output files explanation.html").read_text(encoding="utf-8")
    assert 'href="ZOTU_Pipeline__ZOTUs-Table-For-PICRUSt2.tab"' in page


def test_render_readme_html_escapes_text():
    content = (
        "====\nTACTIC Pipeline — Test <Page>\n====\n\n"
        "Files\n-----\n"
        "• <Pipeline>__a&b.tab\n"
        "  - What it is: uses <script>alert(1)</script> & `--input`\n"
        "====\n"
    )
    page = render_readme_html(content)
    assert "<script>alert(1)</script>" not in page
    assert "&lt;Pipeline&gt;__a&amp;b.tab" in page
    assert "<code>--input</code>" in page


def test_remove_intermediate_files(tmp_path):
    for name in ["test_ZOTUs-Seqs.fasta", "test_ZOTUs-Seqs.csv", "shortened_algn_test_ZOTUs-Seqs.fasta",
                 "ZOTUs-Seqs.fasta", "ZOTUs-Table.tab"]:
        (tmp_path / name).write_text("x", encoding="utf-8")
    removed = sorted(Path(p).name for p in remove_intermediate_files(str(tmp_path)))
    assert removed == ["shortened_algn_test_ZOTUs-Seqs.fasta", "test_ZOTUs-Seqs.csv", "test_ZOTUs-Seqs.fasta"]
    assert sorted(p.name for p in tmp_path.iterdir()) == ["ZOTUs-Seqs.fasta", "ZOTUs-Table.tab"]
    assert remove_intermediate_files(str(tmp_path / "missing")) == []


def test_prepare_picrust2_input_skips_intermediate_files(tmp_path):
    zotu_dir = tmp_path / "ZOTU_Pipeline"
    zotu_dir.mkdir()
    for name in ["ZOTUs-Seqs.fasta", "test_ZOTUs-Seqs.fasta", "shortened_algn_test_ZOTUs-Seqs.fasta"]:
        (zotu_dir / name).write_text(">ZOTU_1 tax=Bacteria;\nACGT\n", encoding="utf-8")
    prepare_picrust2_input(str(tmp_path))
    fastas = sorted(p.name for p in (tmp_path / "PICRUSt2-Inputs").glob("*.fasta"))
    assert fastas == ["ZOTU_Pipeline__ZOTUs-Seqs-For-PICRUSt2.fasta"]
