from pathlib import Path
from tactic_pipeline.analysis.analysis_helper import (
    get_pipeline_readme_content,
    write_pipeline_readme,
    write_all_pipeline_readmes,
    prepare_picrust2_input,
)


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
    assert readme_path.name == "README.txt"
    text = readme_path.read_text(encoding="utf-8")
    assert "TACTIC Pipeline — OTU Pipeline Output Files" in text


def test_write_all_pipeline_readmes_with_active_pipelines(tmp_path):
    active = ["otu", "tac"]
    created = write_all_pipeline_readmes(str(tmp_path), active_pipelines=active)
    assert len(created) == 3
    assert all(Path(p).is_file() for p in created)

    assert (tmp_path / "README.txt").is_file()
    assert (tmp_path / "OTU_Pipeline" / "README.txt").is_file()
    assert (tmp_path / "TAC_Pipeline" / "README.txt").is_file()
    # Inactive and empty should not have README.txt
    assert not (tmp_path / "TIC_Pipeline" / "README.txt").exists()

    otu_text = (tmp_path / "OTU_Pipeline" / "README.txt").read_text(encoding="utf-8")
    assert "OTU Pipeline Output Files" in otu_text
    assert "OTUs-Table.tab" in otu_text

    tac_text = (tmp_path / "TAC_Pipeline" / "README.txt").read_text(encoding="utf-8")
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

    assert (tmp_path / "README.txt").is_file()
    assert (tic_dir / "README.txt").is_file()
    assert (shared_dir / "README.txt").is_file()
    assert not (empty_zotu_dir / "README.txt").exists()

    tic_text = (tic_dir / "README.txt").read_text(encoding="utf-8")
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
    readme_file = picrust_dir / "README.txt"
    assert readme_file.is_file()
    text = readme_file.read_text(encoding="utf-8")
    assert "PICRUSt2 Input Files" in text
    assert "picrust2_pipeline.py" in text
