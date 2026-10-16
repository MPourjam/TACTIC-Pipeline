from pathlib import Path
from ticlust import TICAnalysis

from tactic_pipeline.analysis.analysis import (
    addTax_to_table,
    cleanup as analysis_cleanup,
)
from tactic_pipeline.analysis.analysis_helper import (
    build_cluster_id_kingdom_map,
    build_kingdom_lookup_from_sina_csv,
    detect_silva_version,
    ensure_7level_semicolons_in_fasta,
    ensure_7level_semicolons_in_table,
    format_taxonomy_levels,
    generate_kingdom_enriched_pipeline_outputs,
    is_silva_144_or_newer,
    parse_sina_taxonomy,
    prepare_picrust2_input,
    restore_kingdom_in_taxonomy,
    write_taxonomy_levels_readme,
)


def test_detect_silva_version_and_is_silva_144_or_newer(tmp_path):
    assert detect_silva_version("/databases/SILVA/144") == (144,)
    assert detect_silva_version("/databases/SILVA/144/SILVA.arb") == (144,)
    assert detect_silva_version("/databases/SILVA/138.2") == (138, 2)
    assert detect_silva_version("/tmp/custom_db") is None

    assert is_silva_144_or_newer("/databases/SILVA/144") is True
    assert is_silva_144_or_newer("/databases/SILVA/138.2") is False

    # Fallback content-based detection when path has no version
    v144_samples = [
        "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;metagenome;",
    ]
    v144_other_samples = [
        "Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Odoribacteraceae;Odoribacter;",
    ]
    v138_samples = [
        "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;",
    ]
    assert is_silva_144_or_newer("/tmp/custom_db", sample_taxonomies=v144_samples) is True
    assert is_silva_144_or_newer("/tmp/custom_db", sample_taxonomies=v144_other_samples) is True
    assert is_silva_144_or_newer("/tmp/custom_db", sample_taxonomies=v138_samples) is False


def test_parse_sina_taxonomy_pre144_and_144_plus():
    # Pre-144 (6 ranks from SINA -> 7 levels, 7 semicolons)
    raw_138 = "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;"
    tax_138, dom_138, kng_138 = parse_sina_taxonomy(raw_138, is_v144_plus=False)
    assert tax_138 == "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;;"
    assert tax_138.count(";") == 7
    assert dom_138 == "Bacteria"
    assert kng_138 == ""

    # SILVA 144 with 8th source level (metagenome) and -ati kingdom
    raw_144 = (
        "Bacteria;Bacillati;Clostridiota;Clostridia;"
        "Acetivibrionales;Lachnospiraceae;Incertae Sedis;metagenome"
    )
    tax_144, dom_144, kng_144 = parse_sina_taxonomy(raw_144, is_v144_plus=True)
    assert tax_144 == "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;;"
    assert tax_144.count(";") == 7
    assert dom_144 == "Bacteria"
    assert kng_144 == "Bacillati"

    # SILVA 144 with --other kingdom and trailing semicolon
    raw_144_other = (
        "Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;"
        "Bacteroidales;Odoribacteraceae;Odoribacter;"
    )
    tax_other, dom_other, kng_other = parse_sina_taxonomy(raw_144_other, is_v144_plus=True)
    assert tax_other == "Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Odoribacteraceae;Odoribacter;;"
    assert tax_other.count(";") == 7
    assert dom_other == "Bacteria"
    assert kng_other == "Pseudomonadati--other"

    # Empty / unclassified
    empty_tax, _, _ = parse_sina_taxonomy("", is_v144_plus=True)
    assert empty_tax == ";;;;;;;"
    assert empty_tax.count(";") == 7


def test_restore_kingdom_in_7level_and_8level_outputs(tmp_path):
    sina_csv = tmp_path / "test_ZOTUs-Seqs.csv"
    sina_csv.write_text(
        "name,full_name,align_quality_slv,align_ident_slv,align_score_slv,lca_tax_embl_ebi_ena,lca_tax_slv\n"
        "Zotu1,Zotu1,100,99.5,1.0,,Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;metagenome;\n"
        "Zotu2,Zotu2,100,99.1,1.0,,Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Odoribacteraceae;Odoribacter;\n"
        "Zotu3,Zotu3,100,95.0,1.0,,Bacteria;Fusobacteriati;;;;;;\n"
        "Zotu4,Zotu4,100,80.0,1.0,,Unclassified;\n",
        encoding="utf-8",
    )

    seq_k_map, phylum_k_map = build_kingdom_lookup_from_sina_csv(str(sina_csv))
    assert seq_k_map["Zotu1"] == "Bacillati"
    assert seq_k_map["Zotu2"] == "Pseudomonadati--other"
    assert seq_k_map["Zotu3"] == "Fusobacteriati"
    assert phylum_k_map[("Bacteria", "Clostridiota")] == "Bacillati"
    assert phylum_k_map[("Bacteria", "Bacteroidota")] == "Pseudomonadati--other"

    # 7-level restoration for ZOTU/OTU (Domain..Genus, 7 semicolons)
    zotu1_classic = "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;;"
    zotu1_7lvl = restore_kingdom_in_taxonomy(
        zotu1_classic,
        seq_id="Zotu1",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=7,
    )
    assert zotu1_7lvl == "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Incertae Sedis;"
    assert zotu1_7lvl.count(";") == 7

    # Kingdom-only classified ZOTU3 restored via seq_id map
    zotu3_classic = "Bacteria;;;;;;;"
    zotu3_7lvl = restore_kingdom_in_taxonomy(
        zotu3_classic,
        seq_id="Zotu3",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=7,
    )
    assert zotu3_7lvl == "Bacteria;Fusobacteriati;;;;;;"
    assert zotu3_7lvl.count(";") == 7

    # 8-level restoration for TIC SOTU (Domain..SOTU, 8 semicolons)
    sotu_map_file = tmp_path / "Map-SOTU-ZOTU.tab"
    sotu_map_file.write_text(
        "SOTU\tZOTU\n"
        "SOTU101\tZotu1\n"
        "SOTU201\tZotu2\n",
        encoding="utf-8",
    )
    tic_k_map = build_cluster_id_kingdom_map(str(sotu_map_file), seq_k_map)
    sotu1_classic = "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;GOTU101;SOTU101;"
    sotu1_8lvl = restore_kingdom_in_taxonomy(
        sotu1_classic,
        seq_id="SOTU101",
        seq_kingdom_map=tic_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=8,
    )
    assert sotu1_8lvl == (
        "Bacteria;Bacillati;Clostridiota;Clostridia;"
        "Acetivibrionales;Lachnospiraceae;GOTU101;SOTU101;"
    )
    assert sotu1_8lvl.count(";") == 8

    # FullTaxonomy NA-Kingdom expansion -> NA-Domain;NA-Kingdom;... (8 semicolons)
    na_classic = "NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species;"
    na_8lvl = restore_kingdom_in_taxonomy(
        na_classic,
        seq_id="Zotu4",
        seq_kingdom_map=tic_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=8,
    )
    assert na_8lvl == "NA-Domain;NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species;"
    assert na_8lvl.count(";") == 8


def test_tic_analysis_and_dual_output_generation_with_silva_144(tmp_path):
    zotu_dir = tmp_path / "ZOTU_Pipeline"
    tic_dir = tmp_path / "TIC_Pipeline"
    zotu_dir.mkdir()
    tic_dir.mkdir()

    # Simulate SILVA 144 SINA output and stripped 7-level ZOTU files
    sina_csv = zotu_dir / "test_ZOTUs-Seqs.csv"
    sina_csv.write_text(
        "name,full_name,align_quality_slv,align_ident_slv,align_score_slv,lca_tax_embl_ebi_ena,lca_tax_slv\n"
        "Zotu1,Zotu1,100,99.5,1.0,,Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;metagenome;\n"
        "Zotu2,Zotu2,100,99.1,1.0,,Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;\n",
        encoding="utf-8",
    )
    seq_k_map, phylum_k_map = build_kingdom_lookup_from_sina_csv(str(sina_csv))

    tax1, _, _ = parse_sina_taxonomy(
        "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;metagenome;",
        is_v144_plus=True,
    )
    tax2, _, _ = parse_sina_taxonomy(
        "Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;",
        is_v144_plus=True,
    )

    seq1 = "TACGTAGGGGGCAAGCGTTATCCGGATTTACTGGGTGTAAAGGGAGCGTAGACGG"
    seq2 = "TGAGGAAGGTCCCCCACATTGGAACTGAGACACGGTCCAAACTCCTACGGGAGGC"

    zotu_fasta = zotu_dir / "ZOTUs-Seqs.fasta"
    zotu_fasta.write_text(
        f">Zotu1 tax={tax1}\n{seq1}\n"
        f">Zotu2 tax={tax2}\n{seq2}\n",
        encoding="utf-8",
    )
    zotu_table = zotu_dir / "ZOTUs-Table.tab"
    zotu_table.write_text(
        "#OTU ID\tSampleA\tSampleB\n"
        "Zotu1\t100.0\t50.0\n"
        "Zotu2\t20.0\t80.0\n",
        encoding="utf-8",
    )
    addTax_to_table(str(zotu_table), str(zotu_fasta))

    # Generate 7-level outputs for ZOTU_Pipeline
    generate_kingdom_enriched_pipeline_outputs(
        str(zotu_dir),
        target_levels=7,
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
    )
    assert (zotu_dir / "ZOTUs-Table-7Levels.tab").is_file()
    assert (zotu_dir / "ZOTUs-Seqs-7Levels.fasta").is_file()
    zotu_7lvl_text = (zotu_dir / "ZOTUs-Table-7Levels.tab").read_text(encoding="utf-8")
    assert "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;" in zotu_7lvl_text
    assert "Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;" in zotu_7lvl_text

    # Stage into TIC_Pipeline and run TICAnalysis directly
    staged_fasta = tic_dir / "ZOTUs-Seqs.fasta"
    staged_table = tic_dir / "ZOTUs-Table.tab"
    staged_fasta.write_text(zotu_fasta.read_text(encoding="utf-8"), encoding="utf-8")
    staged_table.write_text(zotu_table.read_text(encoding="utf-8"), encoding="utf-8")

    tic_analysis = TICAnalysis(
        taxed_fasta_file_path=staged_fasta,
        zotu_table_file=staged_table,
    )
    tic_analysis.run(threads=1)

    for tab_f in sorted(tic_analysis.tic_wd.glob("*.tab")):
        if not tab_f.name.startswith("Map-"):
            ensure_7level_semicolons_in_table(str(tab_f))
    for fa_f in sorted(tic_analysis.tic_wd.glob("*.fasta")):
        ensure_7level_semicolons_in_fasta(str(fa_f))

    (tic_analysis.tic_wd / "Non-Bacteria-Sequences.fasta").rename(tic_dir / "Invalid-Tax-Seqs-TIC.fasta")
    (tic_analysis.tic_wd / "SOTUs-Seq.fasta").rename(tic_dir / "SOTUs-Seqs-TIC.fasta")
    (tic_analysis.tic_wd / "SOTUs-Table.tab").rename(tic_dir / "SOTUs-Table-TIC.tab")
    for item in list(tic_analysis.tic_wd.iterdir()):
        item.rename(tic_dir / item.name)
    tic_analysis.tic_wd.rmdir()

    tic_k_map = build_cluster_id_kingdom_map(
        str(tic_dir / "Map-SOTU-ZOTU.tab"),
        seq_k_map,
    )
    generate_kingdom_enriched_pipeline_outputs(
        str(tic_dir),
        target_levels=8,
        seq_kingdom_map=tic_k_map,
        phylum_kingdom_map=phylum_k_map,
    )

    # Staged ZOTU inputs inside TIC_Pipeline should not be duplicated as -8Levels
    assert not (tic_dir / "ZOTUs-Table-8Levels.tab").exists()
    assert not (tic_dir / "ZOTUs-Seqs-8Levels.fasta").exists()

    # TIC 7-level standard table has 7 semicolons per taxonomy
    sotu_std_lines = (tic_dir / "SOTUs-Table-TIC.tab").read_text(encoding="utf-8").strip().splitlines()[1:]
    assert len(sotu_std_lines) == 2
    for row in sotu_std_lines:
        tax_col = row.split("\t")[-1]
        assert tax_col.count(";") == 7

    # TIC 8-level enriched table has 8 semicolons per taxonomy and includes Kingdom + SOTU
    sotu_8lvl_text = (tic_dir / "SOTUs-Table-TIC-8Levels.tab").read_text(encoding="utf-8")
    sotu_8lvl_lines = sotu_8lvl_text.strip().splitlines()[1:]
    assert len(sotu_8lvl_lines) == 2
    for row in sotu_8lvl_lines:
        tax_col = row.split("\t")[-1]
        assert tax_col.count(";") == 8
    assert "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;SOTU" in sotu_8lvl_text
    assert "Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;SOTU" in sotu_8lvl_text

    # PICRUSt2 preparation ignores -7Levels and -8Levels duplicates
    prepare_picrust2_input(str(tmp_path))
    picrust_files = sorted(p.name for p in (tmp_path / "PICRUSt2-Inputs").iterdir())
    assert not any("7Levels" in name or "8Levels" in name for name in picrust_files)

    # Cleanup preserves both standard and -7Levels/-8Levels files
    analysis_cleanup(str(zotu_dir))
    analysis_cleanup(str(tic_dir))
    assert (zotu_dir / "ZOTUs-Table-7Levels.tab").is_file()
    assert (tic_dir / "SOTUs-Table-TIC-8Levels.tab").is_file()
    assert (tic_dir / "ZOTU-Table-All-FullTaxonomy-8Levels.tab").is_file()


def test_detect_silva_version_edge_cases(tmp_path):
    # None, empty, whitespace strings
    assert detect_silva_version(None) is None
    assert detect_silva_version("") is None
    assert detect_silva_version("   ") is None

    # Archive / file extensions
    assert detect_silva_version("/mnt/dbs/SILVA_144_SSURef_NR99_03_07_24_opt.arb.gz") == (144,)
    assert detect_silva_version("/mnt/dbs/silva.145.0.arb") == (145, 0)
    assert detect_silva_version("/databases/SILVA/144.1/SILVA.arb") == (144, 1)

    # Disambiguation: path with numbers in non-SILVA segments
    assert detect_silva_version("/data/experiment_144/ref_dbs/SILVA_138.2/db.arb") == (138, 2)
    assert detect_silva_version("/scratch_2024/batch_144/SILVA/138.1/") == (138, 1)

    # Parent directory containing SILVA/<version>/SILVA.arb structure
    custom_root = tmp_path / "custom_dbs"
    (custom_root / "SILVA" / "144").mkdir(parents=True)
    (custom_root / "SILVA" / "144" / "SILVA.arb").touch()
    assert detect_silva_version(str(custom_root)) == (144,)

    # Multiple subfolders in SILVA/ picks highest version
    (custom_root / "SILVA" / "138.2").mkdir(parents=True)
    (custom_root / "SILVA" / "138.2" / "SILVA.arb").touch()
    assert detect_silva_version(str(custom_root)) == (144,)

    # Non-SILVA reference paths
    assert detect_silva_version("/opt/ref_dbs/Greengenes2/2022.10/gg2.qza") is None
    assert detect_silva_version("/opt/databases/ncbi_16s/") is None


def test_is_silva_144_or_newer_fallback_edge_cases():
    # Empty, None, or whitespace-only inputs
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=[]) is False
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=[None, "", "   "]) is False
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=["Unclassified;", "No hit;"]) is False

    # Single-rank lines
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=["Bacteria;"]) is False
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=["Root;Bacteria"]) is False

    # Classic pre-144 taxonomy (rank 1 is phylum, does not end in -ati)
    classic = ["Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;"]
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=classic) is False

    # SILVA >= 144 Archaea with -ati kingdom
    archaea_144 = [
        "Archaea;Methanobacteriati;Methanobacteriota;Methanobacteria;Methanobacteriales;Methanobacteriaceae;Methanobrevibacter;"
    ]
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=archaea_144) is True

    # Case-insensitive -ati
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=["Bacteria;bacillati;Bacillota;..."]) is True
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=["Bacteria;BACILLATI;Bacillota;..."]) is True

    # Rank 1 is NOT kingdom, but later ranks contain -ati
    false_alarm = ["Bacteria;Firmicutes;Someclassati;Bacillales;Familyati;Bacillus;"]
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=false_alarm) is False

    # Mixed list where one line is unclassified and another has -ati
    mixed = ["Unclassified;", ";;;;;;", "Bacteria;Spirochaetati;Spirochaetota;Spirochaetia;"]
    assert is_silva_144_or_newer("/tmp/no_ver", sample_taxonomies=mixed) is True


def test_parse_sina_taxonomy_deep_edge_cases():
    # Null and blank strings
    for is_v144 in (False, True):
        tax, dom, kng = parse_sina_taxonomy(None, is_v144_plus=is_v144)
        assert tax == ";;;;;;;"
        assert tax.count(";") == 7
        assert dom == ""
        assert kng == ""

        tax, dom, kng = parse_sina_taxonomy("   ", is_v144_plus=is_v144)
        assert tax == ";;;;;;;"
        assert tax.count(";") == 7

        tax, dom, kng = parse_sina_taxonomy(";;;;", is_v144_plus=is_v144)
        assert tax == ";;;;;;;"
        assert tax.count(";") == 7

    # Deeply incomplete taxonomies with SILVA >= 144
    tax1, dom1, kng1 = parse_sina_taxonomy("Bacteria", is_v144_plus=True)
    assert tax1 == "Bacteria;;;;;;;"
    assert tax1.count(";") == 7
    assert dom1 == "Bacteria"
    assert kng1 == ""

    tax2, dom2, kng2 = parse_sina_taxonomy("Bacteria;Bacillati", is_v144_plus=True)
    assert tax2 == "Bacteria;;;;;;;"
    assert tax2.count(";") == 7
    assert dom2 == "Bacteria"
    assert kng2 == "Bacillati"

    tax3, dom3, kng3 = parse_sina_taxonomy("Bacteria;Bacillati;Firmicutes", is_v144_plus=True)
    assert tax3 == "Bacteria;Firmicutes;;;;;;"
    assert tax3.count(";") == 7
    assert dom3 == "Bacteria"
    assert kng3 == "Bacillati"

    tax4, dom4, kng4 = parse_sina_taxonomy("Bacteria;Bacillati;Firmicutes;Bacilli", is_v144_plus=True)
    assert tax4 == "Bacteria;Firmicutes;Bacilli;;;;;"
    assert tax4.count(";") == 7
    assert dom4 == "Bacteria"
    assert kng4 == "Bacillati"

    # Whitespace surrounding each token
    messy_144 = " Bacteria ;  Bacillati ; Firmicutes ; Bacilli ; Bacillales ; Bacillaceae ; Bacillus ; metagenome "
    clean_tax, clean_dom, clean_kng = parse_sina_taxonomy(messy_144, is_v144_plus=True)
    assert clean_tax == "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;"
    assert clean_tax.count(";") == 7
    assert clean_dom == "Bacteria"
    assert clean_kng == "Bacillati"

    # Excessive tokens (>8 ranks in SILVA >= 144)
    long_144 = (
        "Bacteria;Bacillati;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;"
        "subsp_subtilis;metagenome;extra_flag_1;extra_flag_2"
    )
    long_tax, long_dom, long_kng = parse_sina_taxonomy(long_144, is_v144_plus=True)
    assert long_tax == "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;"
    assert long_tax.count(";") == 7
    assert long_dom == "Bacteria"
    assert long_kng == "Bacillati"

    # Excessive tokens in Pre-144 (>6 ranks)
    long_138 = "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;subsp_subtilis;extra_flag"
    p138_tax, p138_dom, p138_kng = parse_sina_taxonomy(long_138, is_v144_plus=False)
    assert p138_tax == "Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;"
    assert p138_tax.count(";") == 7
    assert p138_dom == "Bacteria"
    assert p138_kng == ""


def test_build_kingdom_lookup_from_sina_csv_edge_cases(tmp_path):
    # Non-existent file
    missing = tmp_path / "non_existent.csv"
    assert build_kingdom_lookup_from_sina_csv(str(missing)) == ({}, {})

    # File with empty lines, comments, short lines, and missing taxonomy
    sina_csv = tmp_path / "edge_cases.csv"
    sina_csv.write_text(
        "name,full_name,align_quality_slv,align_ident_slv,align_score_slv,lca_tax_embl_ebi_ena,lca_tax_slv\n"
        "# comment row\n"
        "\n"
        "Zotu1,\"Zotu1, annotated sequence\",100,99.2,1.0,,Bacteria;Bacillati;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;metagenome;\n"
        "Zotu2,Zotu2,100,95.0,1.0,,Bacteria;\n"
        "Zotu3,Zotu3,100,90.0,1.0,,Unclassified;\n"
        "Zotu4,Zotu4,100,85.0,1.0,,\n"
        "truncated_line,100,90.0\n",
        encoding="utf-8",
    )
    seq_k_map, phylum_k_map = build_kingdom_lookup_from_sina_csv(str(sina_csv))
    assert seq_k_map == {"Zotu1": "Bacillati"}
    assert phylum_k_map == {("Bacteria", "Firmicutes"): "Bacillati"}
    assert "Zotu2" not in seq_k_map
    assert "Zotu3" not in seq_k_map
    assert "Zotu4" not in seq_k_map


def test_build_cluster_id_kingdom_map_edge_cases(tmp_path):
    zotu_k_map = {"Zotu1": "Bacillati", "Zotu2": "Pseudomonadati--other"}

    # Missing map file returns copy of zotu map
    missing_map = tmp_path / "no_map.tab"
    assert build_cluster_id_kingdom_map(str(missing_map), zotu_k_map) == zotu_k_map

    # Map file with comments, various headers, and unknown ZOTU entries
    map_file = tmp_path / "Map-Test.tab"
    map_file.write_text(
        "# Cluster mapping file\n"
        "OTU\tZOTU\n"
        "\n"
        "OTU1\tZotu1\n"
        "OTU2\tZotu2\n"
        "OTU3\tZotu_Unclassified\n"
        "OTU1\tZotu99_Duplicate\n",
        encoding="utf-8",
    )
    cluster_k_map = build_cluster_id_kingdom_map(str(map_file), zotu_k_map)
    assert cluster_k_map["OTU1"] == "Bacillati"
    assert cluster_k_map["OTU2"] == "Pseudomonadati--other"
    assert "OTU3" not in cluster_k_map  # Zotu_Unclassified had no kingdom


def test_restore_kingdom_in_taxonomy_edge_cases():
    seq_k_map = {"Zotu1": "Bacillati"}
    phylum_k_map = {("Bacteria", "Clostridiota"): "Bacillati"}

    # Fallback to phylum_kingdom_map when seq_id is not in seq_kingdom_map
    classic_tax = "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;;"
    restored_by_phylum = restore_kingdom_in_taxonomy(
        classic_tax,
        seq_id="Zotu_Unseen",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=7,
    )
    assert restored_by_phylum == "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;"
    assert restored_by_phylum.count(";") == 7

    # Handles seq_id with leading '>' character
    restored_by_fasta_id = restore_kingdom_in_taxonomy(
        classic_tax,
        seq_id=">Zotu1",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=7,
    )
    assert restored_by_fasta_id == "Bacteria;Bacillati;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;"

    # Neither seq_id nor phylum is known: preserves slot as empty with exact semicolon count
    unknown_tax = "Bacteria;UnknownPhylum;UnknownClass;UnknownOrder;UnknownFamily;UnknownGenus;;"
    restored_unknown_7 = restore_kingdom_in_taxonomy(
        unknown_tax,
        seq_id="Zotu_Unknown",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=7,
    )
    assert restored_unknown_7 == "Bacteria;;UnknownPhylum;UnknownClass;UnknownOrder;UnknownFamily;UnknownGenus;"
    assert restored_unknown_7.count(";") == 7

    restored_unknown_8 = restore_kingdom_in_taxonomy(
        unknown_tax,
        seq_id="Zotu_Unknown",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=8,
    )
    assert restored_unknown_8 == "Bacteria;;UnknownPhylum;UnknownClass;UnknownOrder;UnknownFamily;UnknownGenus;;"
    assert restored_unknown_8.count(";") == 8

    # Empty taxonomy input
    assert restore_kingdom_in_taxonomy("", target_levels=7) == ";;;;;;;"
    assert restore_kingdom_in_taxonomy("", target_levels=8) == ";;;;;;;;"

    # Incomplete taxonomy (padded up to target levels with kingdom inserted at index 1)
    short_tax = "Bacteria;Clostridiota;"
    restored_short_7 = restore_kingdom_in_taxonomy(
        short_tax,
        seq_id="Zotu1",
        seq_kingdom_map=seq_k_map,
        phylum_kingdom_map=phylum_k_map,
        target_levels=7,
    )
    assert restored_short_7 == "Bacteria;Bacillati;Clostridiota;;;;;"
    assert restored_short_7.count(";") == 7

    # NA-taxonomies: sequence with known kingdom in seq_k_map
    na_tax = "NA-Kingdom;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species;"
    restored_na_known = restore_kingdom_in_taxonomy(
        na_tax,
        seq_id="Zotu1",
        seq_kingdom_map=seq_k_map,
        target_levels=8,
    )
    assert restored_na_known == "NA-Domain;Bacillati;NA-Phylum;NA-Class;NA-Order;NA-Family;NA-Genus;NA-Species;"
    assert restored_na_known.count(";") == 8


def test_ensure_7level_semicolons_in_table_and_fasta_edge_cases(tmp_path):
    # Non-existent files do not raise exceptions
    ensure_7level_semicolons_in_table(str(tmp_path / "no_table.tab"))
    ensure_7level_semicolons_in_fasta(str(tmp_path / "no_fasta.fasta"))

    # Table without Taxonomy column (e.g. sample count matrix) is unmodified
    count_table = tmp_path / "counts.tab"
    original_counts = "#OTU ID\tSampleA\tSampleB\nOTU1\t10\t20\n"
    count_table.write_text(original_counts, encoding="utf-8")
    ensure_7level_semicolons_in_table(str(count_table))
    assert count_table.read_text(encoding="utf-8") == original_counts

    # Table with uppercase TAXONOMY column and 6-semicolon entries converted to 7 semicolons
    tax_table = tmp_path / "tax_table.tab"
    tax_table.write_text(
        "#OTU ID\tSampleA\tTAXONOMY\n"
        "# Comment row inside table\n"
        "Zotu1\t50\tBacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;SpeciesA\n",
        encoding="utf-8",
    )
    ensure_7level_semicolons_in_table(str(tax_table))
    updated_lines = tax_table.read_text(encoding="utf-8").strip().splitlines()
    assert updated_lines[2].split("\t")[2] == "Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;SpeciesA;"
    assert updated_lines[2].split("\t")[2].count(";") == 7

    # Running ensure_7level_semicolons_in_table again is idempotent
    ensure_7level_semicolons_in_table(str(tax_table))
    idempotent_lines = tax_table.read_text(encoding="utf-8").strip().splitlines()
    assert idempotent_lines[2].split("\t")[2].count(";") == 7

    # FASTA with mixed headers (some with tax=, some without, some with extra attributes)
    fasta_file = tmp_path / "test_mixed.fasta"
    fasta_file.write_text(
        ">Seq_no_tax\nATCG\n"
        ">Zotu1 tax=Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;\nATCG\n"
        ">Zotu2 tax=Bacteria;Firmicutes; size=42;\nATCG\n",
        encoding="utf-8",
    )
    ensure_7level_semicolons_in_fasta(str(fasta_file))
    fasta_lines = fasta_file.read_text(encoding="utf-8").strip().splitlines()
    assert fasta_lines[0] == ">Seq_no_tax"
    assert "tax=Bacteria;Clostridiota;Clostridia;Acetivibrionales;Lachnospiraceae;Roseburia;; " in fasta_lines[2] or fasta_lines[2].endswith(";;")
    assert fasta_lines[2].count(";") == 7
    assert "tax=Bacteria;Firmicutes;;;;;; size=42;" in fasta_lines[4]


def test_tac_pipeline_full_enrichment_flow(tmp_path):
    tac_dir = tmp_path / "TAC_Pipeline"
    tac_dir.mkdir()

    seq_k_map = {"Zotu1": "Bacillati", "Zotu2": "Pseudomonadati--other"}
    phylum_k_map = {("Bacteria", "Firmicutes"): "Bacillati"}

    # Simulate Map-ZOTUs-OTUs-TAC.tab
    map_file = tac_dir / "Map-ZOTUs-OTUs-TAC.tab"
    map_file.write_text(
        "ZOTU\tOTU\n"
        "Zotu1\tOTU1\n"
        "Zotu2\tOTU2\n",
        encoding="utf-8",
    )
    tac_k_map = build_cluster_id_kingdom_map(str(map_file), seq_k_map)
    assert tac_k_map["OTU1"] == "Bacillati"
    assert tac_k_map["OTU2"] == "Pseudomonadati--other"

    # Simulate standard TAC outputs (7 levels, 7 semicolons)
    tac_table = tac_dir / "OTUs-Table-TAC.tab"
    tac_table.write_text(
        "#OTU ID\tSample1\tTaxonomy\n"
        "OTU1\t100\tBacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;\n"
        "OTU2\t50\tBacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;;\n",
        encoding="utf-8",
    )
    tac_fasta = tac_dir / "OTUs-Seqs-TAC.fasta"
    tac_fasta.write_text(
        ">OTU1 tax=Bacteria;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;\nATGC\n"
        ">OTU2 tax=Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;;\nGCAT\n",
        encoding="utf-8",
    )
    krona_html = tac_dir / "OTUs-Krona.html"
    krona_html.write_text("<html>Krona 7-level</html>", encoding="utf-8")

    # Staged ZOTU files from preprocessing inside TAC_Pipeline
    staged_zotu_table = tac_dir / "ZOTUs-Table.tab"
    staged_zotu_table.write_text(tac_table.read_text(encoding="utf-8"), encoding="utf-8")
    staged_zotu_fasta = tac_dir / "ZOTUs-Seqs.fasta"
    staged_zotu_fasta.write_text(tac_fasta.read_text(encoding="utf-8"), encoding="utf-8")

    # Mock Krona creator function
    created_krona_calls = []
    def mock_create_krona(tool, table_path, out_html):
        created_krona_calls.append((tool, table_path, out_html))
        (tac_dir / out_html).write_text("<html>Krona 8-level</html>", encoding="utf-8")

    generate_kingdom_enriched_pipeline_outputs(
        str(tac_dir),
        target_levels=8,
        seq_kingdom_map=tac_k_map,
        phylum_kingdom_map=phylum_k_map,
        krona_tool="ktImportText",
        create_krona_fn=mock_create_krona,
    )

    # Verify generated files
    assert (tac_dir / "OTUs-Table-TAC-8Levels.tab").is_file()
    assert (tac_dir / "OTUs-Seqs-TAC-8Levels.fasta").is_file()
    assert (tac_dir / "OTUs-Krona-8Levels.html").is_file()
    assert len(created_krona_calls) == 1
    assert created_krona_calls[0][0] == "ktImportText"
    assert created_krona_calls[0][2] == "OTUs-Krona-8Levels.html"

    # Verify staged ZOTU inputs and Map files were NOT enriched as -8Levels
    assert not (tac_dir / "ZOTUs-Table-8Levels.tab").exists()
    assert not (tac_dir / "ZOTUs-Seqs-8Levels.fasta").exists()
    assert not (tac_dir / "Map-ZOTUs-OTUs-TAC-8Levels.tab").exists()

    # Check 8-level taxonomies and semicolon counts in table
    enriched_tab_text = (tac_dir / "OTUs-Table-TAC-8Levels.tab").read_text(encoding="utf-8")
    for row in enriched_tab_text.strip().splitlines()[1:]:
        tax = row.split("\t")[-1]
        assert tax.count(";") == 8
    assert "Bacteria;Bacillati;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;" in enriched_tab_text
    assert "Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;;" in enriched_tab_text

    # Check 8-level taxonomies and semicolon counts in FASTA
    enriched_fa_text = (tac_dir / "OTUs-Seqs-TAC-8Levels.fasta").read_text(encoding="utf-8")
    assert ">OTU1 tax=Bacteria;Bacillati;Firmicutes;Bacilli;Bacillales;Bacillaceae;Bacillus;;" in enriched_fa_text
    assert ">OTU2 tax=Bacteria;Pseudomonadati--other;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Bacteroides;;" in enriched_fa_text


def test_write_taxonomy_levels_readme(tmp_path):
    readme_path_str = write_taxonomy_levels_readme(str(tmp_path))
    readme_path = Path(readme_path_str)
    assert readme_path.is_file()
    assert readme_path.name == "README_Taxonomy_Levels.txt"
    text = readme_path.read_text(encoding="utf-8")
    assert "SILVA reference database version >= 144" in text
    assert "Standard Outputs" in text
    assert "*-8Levels.*" in text
    assert "*-7Levels.*" in text
    assert "Which file should I use?" in text
    assert "PICRUSt2" in text

