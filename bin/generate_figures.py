#!/usr/bin/env python3

"""
GBS-Sentinel Integrated Surveillance Visualization

This module generates publication-quality figures exclusively from
actual outputs produced by the GBS-Sentinel workflow.

Integrated data sources:
    1. AMRFinderPlus
    2. Virulence-factor screening
    3. MLST
    4. GBS-SBG capsular serotyping
    5. Core-genome phylogeny

No synthetic biological observations are created.

When a particular figure cannot be scientifically generated because
the required real data are unavailable, the figure is skipped and
a corresponding report is written under Reports/.

Outputs:
    PNG
    SVG
    PDF
    TSV summary and matrix files
"""

from __future__ import annotations

import argparse
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, Optional

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

try:
    from ete3 import Tree, TreeStyle, NodeStyle, TextFace

    HAS_ETE3 = True
except ImportError:
    HAS_ETE3 = False


# =============================================================================
# Logging
# =============================================================================

def log(message: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{timestamp}] {message}", flush=True)


def write_skip_report(
    outdir: Path,
    figure_name: str,
    reason: str,
) -> None:
    reports_dir = outdir / "Reports"
    reports_dir.mkdir(parents=True, exist_ok=True)

    report_path = reports_dir / f"{figure_name}_SKIPPED.txt"

    report_path.write_text(
        f"Figure: {figure_name}\n"
        f"Status: NOT GENERATED\n"
        f"Reason: {reason}\n"
        f"Timestamp: {datetime.now().isoformat()}\n"
    )

    log(f"Skipped {figure_name}: {reason}")


# =============================================================================
# Figure output
# =============================================================================

def save_figure(
    fig,
    outdir: Path,
    basename: str,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    for extension in ("png", "svg", "pdf"):
        fig.savefig(
            outdir / f"{basename}.{extension}",
            bbox_inches="tight",
            dpi=300,
        )

    plt.close(fig)

    log(
        f"Saved {basename}.{{png,svg,pdf}}"
    )


# =============================================================================
# Utility functions
# =============================================================================

def safe_text(value: object) -> str:
    if value is None:
        return ""

    if pd.isna(value):
        return ""

    return str(value).strip()


def unique_join(
    values: Iterable[object],
    separator: str = ";",
) -> str:

    cleaned = sorted(
        {
            safe_text(value)
            for value in values
            if safe_text(value)
        }
    )

    return separator.join(cleaned)


# =============================================================================
# AMRFinderPlus
# =============================================================================

def parse_amrfinder_file(
    path: Path,
) -> pd.DataFrame:
    """
    Parse an AMRFinderPlus TSV.

    Expected columns from AMRFinderPlus include:

        Gene symbol
        Sequence name
        Scope
        Element type
        Class
        Subclass
        Method
        % Coverage of reference
        % Identity to reference
        Accession of closest sequence
    """

    columns = [
        "sample",
        "gene",
        "sequence",
        "scope",
        "element_type",
        "class",
        "subclass",
        "method",
        "coverage",
        "identity",
        "accession",
    ]

    empty = pd.DataFrame(columns=columns)

    if not path.exists():
        return empty

    if path.stat().st_size == 0:
        return empty

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
        )
    except Exception as exc:
        log(
            f"WARNING: Failed to read AMRFinderPlus "
            f"file {path}: {exc}"
        )
        return empty

    if df.empty:
        return empty

    normalized = {
        safe_text(column).lower(): column
        for column in df.columns
    }


    def get_column(
        *expected_names: str,
    ) -> pd.Series:

        for expected in expected_names:

            actual = normalized.get(
                expected.lower()
            )

            if actual is not None:
                return (
                    df[actual]
                    .fillna("")
                    .astype(str)
                )

        return pd.Series(
            [""] * len(df),
            index=df.index,
        )

    sample = path.stem

    if sample.endswith("_amrfinder"):
        sample = sample[:-10]

    result = pd.DataFrame(
    {
        "sample": sample,
        "gene": get_column(
            "Element symbol",
            "Gene symbol",
        ),
        "sequence": get_column(
            "Contig id",
            "Sequence name",
        ),
        "scope": get_column("Scope"),
        "element_type": get_column(
            "Type",
            "Element type",
        ),
        "class": get_column("Class"),
        "subclass": get_column("Subclass"),
        "method": get_column("Method"),
        "coverage": get_column(
            "% Coverage of reference"
        ),
        "identity": get_column(
            "% Identity to reference"
        ),
        "accession": get_column(
            "Closest reference accession",
            "Accession of closest sequence",
        ),
    }
)
    
    result = result[
        result["gene"].str.strip() != ""
    ].copy()

    return result


def load_amrfinder(
    directory: Path,
) -> pd.DataFrame:

    files = sorted(
        directory.glob("*_amrfinder.tsv")
    )

    frames = []

    for file_path in files:
        frame = parse_amrfinder_file(
            file_path
        )

        if not frame.empty:
            frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "sample",
                "gene",
                "sequence",
                "scope",
                "element_type",
                "class",
                "subclass",
                "method",
                "coverage",
                "identity",
                "accession",
            ]
        )

    return pd.concat(
        frames,
        ignore_index=True,
    )


# =============================================================================
# Virulence
# =============================================================================

def parse_virulence_file(
    path: Path,
) -> pd.DataFrame:

    empty = pd.DataFrame(
        columns=[
            "sample",
            "gene",
            "product",
            "database",
            "accession",
            "identity",
            "coverage",
        ]
    )

    if not path.exists():
        return empty

    if path.stat().st_size == 0:
        return empty

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
        )
    except Exception as exc:
        log(
            f"WARNING: Failed to read virulence "
            f"file {path}: {exc}"
        )
        return empty

    if df.empty:
        return empty

    columns = {
        safe_text(column).upper(): column
        for column in df.columns
    }

    gene_col = columns.get("GENE")
    product_col = columns.get("PRODUCT")
    database_col = columns.get("DATABASE")
    accession_col = columns.get("ACCESSION")
    identity_col = columns.get("%IDENTITY")
    coverage_col = columns.get("%COVERAGE")

    if gene_col is None:
        return empty

    sample = path.stem

    if sample.endswith("_vf"):
        sample = sample[:-3]

    result = pd.DataFrame(
        {
            "sample": sample,
            "gene": (
                df[gene_col]
                .fillna("")
                .astype(str)
            ),
            "product": (
                df[product_col]
                .fillna("")
                .astype(str)
                if product_col is not None
                else ""
            ),
            "database": (
                df[database_col]
                .fillna("")
                .astype(str)
                if database_col is not None
                else ""
            ),
            "accession": (
                df[accession_col]
                .fillna("")
                .astype(str)
                if accession_col is not None
                else ""
            ),
            "identity": (
                df[identity_col]
                .fillna("")
                .astype(str)
                if identity_col is not None
                else ""
            ),
            "coverage": (
                df[coverage_col]
                .fillna("")
                .astype(str)
                if coverage_col is not None
                else ""
            ),
        }
    )

    result = result[
        result["gene"].str.strip() != ""
    ].copy()

    return result


def load_virulence(
    directory: Path,
) -> pd.DataFrame:

    files = sorted(
        directory.glob("*_vf.tsv")
    )

    frames = []

    for file_path in files:
        frame = parse_virulence_file(
            file_path
        )

        if not frame.empty:
            frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "sample",
                "gene",
                "product",
                "database",
                "accession",
                "identity",
                "coverage",
            ]
        )

    return pd.concat(
        frames,
        ignore_index=True,
    )


# =============================================================================
# MLST
# =============================================================================

def normalize_sample_name(
    value: str,
) -> str:

    name = Path(value).name

    for suffix in (
        ".fasta",
        ".fa",
        ".fna",
        ".fsa",
    ):
        if name.endswith(suffix):
            name = name[: -len(suffix)]

    return name


def parse_mlst_file(
    path: Path,
) -> pd.DataFrame:

    columns = [
        "sample",
        "scheme",
        "ST",
        "alleles",
    ]

    rows = []

    if not path.exists():
        return pd.DataFrame(columns=columns)

    if path.stat().st_size == 0:
        return pd.DataFrame(columns=columns)

    try:
        lines = path.read_text(
            encoding="utf-8",
            errors="replace",
        ).splitlines()
    except Exception as exc:
        log(
            f"WARNING: Failed to read MLST "
            f"file {path}: {exc}"
        )
        return pd.DataFrame(columns=columns)

    for line in lines:

        line = line.strip()

        if not line:
            continue

        if line.startswith("#"):
            continue

        parts = line.split("\t")

        if len(parts) < 3:
            continue

        input_name = parts[0]

        sample = normalize_sample_name(
            input_name
        )

        scheme = parts[1].strip()
        st = parts[2].strip()

        alleles = (
            "\t".join(parts[3:])
            if len(parts) > 3
            else ""
        )

        rows.append(
            {
                "sample": sample,
                "scheme": scheme,
                "ST": st,
                "alleles": alleles,
            }
        )

    return pd.DataFrame(
        rows,
        columns=columns,
    )


def load_mlst(
    directory: Path,
) -> pd.DataFrame:

    files = sorted(
        directory.glob("*_mlst.tsv")
    )

    frames = []

    for file_path in files:

        frame = parse_mlst_file(
            file_path
        )

        if not frame.empty:
            frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "sample",
                "scheme",
                "ST",
                "alleles",
            ]
        )

    result = pd.concat(
        frames,
        ignore_index=True,
    )

    result = result.drop_duplicates(
        subset=["sample"],
        keep="last",
    )

    return result


# =============================================================================
# GBS-SBG Serotyping
# =============================================================================

def parse_serotype_file(
    path: Path,
) -> pd.DataFrame:

    columns = [
        "sample",
        "serotype",
        "best_match",
        "confidence",
        "status",
    ]

    empty = pd.DataFrame(
        columns=columns
    )

    if not path.exists():
        return empty

    if path.stat().st_size == 0:
        return empty

    try:
        df = pd.read_csv(
            path,
            sep="\t",
            dtype=str,
        )
    except Exception as exc:
        log(
            f"WARNING: Failed to read serotype "
            f"file {path}: {exc}"
        )
        return empty

    if df.empty:
        return empty

    required = [
        "sample_id",
        "serotype",
        "best_match",
        "confidence",
        "status",
    ]

    if not all(
        column in df.columns
        for column in required
    ):
        log(
            f"WARNING: Serotype file {path} "
            f"does not contain the expected schema."
        )
        return empty

    result = df[required].copy()

    result.columns = columns

    for column in columns:
        result[column] = (
            result[column]
            .fillna("NA")
            .astype(str)
            .str.strip()
        )

    return result


def load_serotyping(
    directory: Path,
) -> pd.DataFrame:

    files = sorted(
        directory.glob("*_serotype.tsv")
    )

    frames = []

    for file_path in files:

        frame = parse_serotype_file(
            file_path
        )

        if not frame.empty:
            frames.append(frame)

    if not frames:
        return pd.DataFrame(
            columns=[
                "sample",
                "serotype",
                "best_match",
                "confidence",
                "status",
            ]
        )

    result = pd.concat(
        frames,
        ignore_index=True,
    )

    result = result.drop_duplicates(
        subset=["sample"],
        keep="last",
    )

    return result


# =============================================================================
# Matrix builders
# =============================================================================

def build_gene_matrix(
    frame: pd.DataFrame,
) -> pd.DataFrame:

    if frame.empty:
        return pd.DataFrame()

    if "sample" not in frame.columns:
        return pd.DataFrame()

    if "gene" not in frame.columns:
        return pd.DataFrame()

    subset = frame[
        ["sample", "gene"]
    ].copy()

    subset["sample"] = (
        subset["sample"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    subset["gene"] = (
        subset["gene"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    subset = subset[
        (subset["sample"] != "")
        &
        (subset["gene"] != "")
    ]

    if subset.empty:
        return pd.DataFrame()

    subset = subset.drop_duplicates()

    matrix = pd.crosstab(
        subset["sample"],
        subset["gene"],
    )

    return matrix.astype(int)


# =============================================================================
# AMR class matrix
# =============================================================================

def build_amr_class_matrix(
    amr: pd.DataFrame,
) -> pd.DataFrame:

    if amr.empty:
        return pd.DataFrame()

    subset = amr[
        ["sample", "class"]
    ].copy()

    subset["sample"] = (
        subset["sample"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    subset["class"] = (
        subset["class"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    subset = subset[
        (subset["sample"] != "")
        &
        (subset["class"] != "")
        &
        (subset["class"].str.lower() != "na")
    ]

    if subset.empty:
        return pd.DataFrame()

    subset = subset.drop_duplicates()

    return pd.crosstab(
        subset["sample"],
        subset["class"],
    ).astype(int)


# =============================================================================
# Plotting
# =============================================================================

def plot_heatmap(
    matrix: pd.DataFrame,
    title: str,
    outdir: Path,
    basename: str,
    cmap: str = "YlOrRd",
    colorbar_label: str = "Presence",
) -> bool:

    if matrix.empty:
        return False

    fig_width = max(
        8,
        min(24, matrix.shape[1] * 0.45),
    )

    fig_height = max(
        4,
        min(24, matrix.shape[0] * 0.55),
    )

    fig, ax = plt.subplots(
        figsize=(fig_width, fig_height)
    )

    sns.heatmap(
        matrix,
        cmap=cmap,
        linewidths=0.4,
        linecolor="white",
        annot=(
            matrix.shape[0] <= 25
            and matrix.shape[1] <= 30
        ),
        fmt="g",
        cbar_kws={
            "label": colorbar_label
        },
        ax=ax,
    )

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel("")
    ax.set_ylabel("Sample")

    save_figure(
        fig,
        outdir,
        basename,
    )

    matrix.to_csv(
        outdir / f"{basename}_matrix.tsv",
        sep="\t",
    )

    return True


def plot_frequency(
    matrix: pd.DataFrame,
    title: str,
    outdir: Path,
    basename: str,
) -> bool:

    if matrix.empty:
        return False

    frequency = (
        matrix.sum(axis=0)
        .sort_values(ascending=False)
    )

    counts = frequency.reset_index()
    counts.columns = [
        "feature",
        "sample_count",
    ]

    counts.to_csv(
        outdir / f"{basename}_counts.tsv",
        sep="\t",
        index=False,
    )

    width = max(
        8,
        min(24, len(frequency) * 0.45),
    )

    fig, ax = plt.subplots(
        figsize=(width, 5)
    )

    frequency.plot(
        kind="bar",
        ax=ax,
        edgecolor="black",
    )

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel("")
    ax.set_ylabel("Number of samples")
    ax.tick_params(
        axis="x",
        rotation=60,
    )

    save_figure(
        fig,
        outdir,
        basename,
    )

    return True


def plot_sample_counts(
    matrix: pd.DataFrame,
    title: str,
    outdir: Path,
    basename: str,
) -> bool:

    if matrix.empty:
        return False

    counts = (
        matrix.sum(axis=1)
        .sort_values(ascending=False)
    )

    counts.to_csv(
        outdir / f"{basename}_counts.tsv",
        sep="\t",
        header=["feature_count"],
    )

    fig, ax = plt.subplots(
        figsize=(
            max(7, len(counts) * 0.65),
            5,
        )
    )

    counts.plot(
        kind="bar",
        ax=ax,
        edgecolor="black",
    )

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel("Sample")
    ax.set_ylabel("Number of detected features")
    ax.tick_params(
        axis="x",
        rotation=45,
    )

    save_figure(
        fig,
        outdir,
        basename,
    )

    return True


def plot_cooccurrence(
    matrix: pd.DataFrame,
    title: str,
    outdir: Path,
    basename: str,
) -> bool:

    if matrix.empty:
        return False

    qualifying = [
        gene
        for gene in matrix.columns
        if matrix[gene].sum() >= 2
    ]

    if len(qualifying) < 2:
        return False

    pairs = []

    for gene_1, gene_2 in combinations(
        qualifying,
        2,
    ):

        co_occurrence = int(
            (
                (matrix[gene_1] == 1)
                &
                (matrix[gene_2] == 1)
            ).sum()
        )

        if co_occurrence > 0:

            pairs.append(
                {
                    "gene_1": gene_1,
                    "gene_2": gene_2,
                    "co_occurrence": co_occurrence,
                }
            )

    if not pairs:
        return False

    pair_df = pd.DataFrame(
        pairs
    )

    pair_df.to_csv(
        outdir / f"{basename}_pairs.tsv",
        sep="\t",
        index=False,
    )

    genes = sorted(
        set(pair_df["gene_1"])
        |
        set(pair_df["gene_2"])
    )

    matrix_out = pd.DataFrame(
        0,
        index=genes,
        columns=genes,
    )

    for _, row in pair_df.iterrows():

        g1 = row["gene_1"]
        g2 = row["gene_2"]
        value = int(
            row["co_occurrence"]
        )

        matrix_out.loc[g1, g2] = value
        matrix_out.loc[g2, g1] = value

    fig, ax = plt.subplots(
        figsize=(
            max(7, len(genes) * 0.5),
            max(7, len(genes) * 0.5),
        )
    )

    sns.heatmap(
        matrix_out,
        annot=True,
        fmt="d",
        cmap="Blues",
        linewidths=0.3,
        ax=ax,
    )

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
    )

    save_figure(
        fig,
        outdir,
        basename,
    )

    return True


def plot_category_distribution(
    data: pd.DataFrame,
    column: str,
    title: str,
    outdir: Path,
    basename: str,
) -> bool:

    if data.empty:
        return False

    if column not in data.columns:
        return False

    values = (
        data[column]
        .fillna("NA")
        .astype(str)
        .str.strip()
    )

    values = values[
        values != ""
    ]

    if values.empty:
        return False

    counts = (
        values
        .value_counts()
        .sort_values(ascending=False)
    )

    counts.to_csv(
        outdir / f"{basename}_counts.tsv",
        sep="\t",
        header=["sample_count"],
    )

    fig, ax = plt.subplots(
        figsize=(
            max(7, len(counts) * 0.7),
            5,
        )
    )

    counts.plot(
        kind="bar",
        ax=ax,
        edgecolor="black",
    )

    ax.set_title(
        title,
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel(column)
    ax.set_ylabel("Number of samples")
    ax.tick_params(
        axis="x",
        rotation=45,
    )

    save_figure(
        fig,
        outdir,
        basename,
    )

    return True


def plot_cross_tab(
    data: pd.DataFrame,
    row_column: str,
    column_column: str,
    title: str,
    outdir: Path,
    basename: str,
) -> bool:

    if data.empty:
        return False

    if (
        row_column not in data.columns
        or column_column not in data.columns
    ):
        return False

    subset = data[
        [row_column, column_column]
    ].copy()

    subset[row_column] = (
        subset[row_column]
        .fillna("NA")
        .astype(str)
    )

    subset[column_column] = (
        subset[column_column]
        .fillna("NA")
        .astype(str)
    )

    matrix = pd.crosstab(
        subset[row_column],
        subset[column_column],
    )

    if matrix.empty:
        return False

    return plot_heatmap(
        matrix,
        title,
        outdir,
        basename,
        cmap="YlGnBu",
        colorbar_label="Number of samples",
    )


# =============================================================================
# Surveillance summary
# =============================================================================

def build_surveillance_summary(
    mlst: pd.DataFrame,
    serotyping: pd.DataFrame,
    amr: pd.DataFrame,
    virulence: pd.DataFrame,
) -> pd.DataFrame:

    sample_names = set()

    for data in (
        mlst,
        serotyping,
        amr,
        virulence,
    ):

        if (
            not data.empty
            and "sample" in data.columns
        ):

            sample_names.update(
                data["sample"]
                .dropna()
                .astype(str)
                .tolist()
            )

    if not sample_names:
        return pd.DataFrame()

    summary = pd.DataFrame(
        {
            "sample": sorted(
                sample_names
            )
        }
    )

    # -----------------------------------------------------------------
    # MLST
    # -----------------------------------------------------------------

    if not mlst.empty:

        mlst_simple = (
            mlst[
                [
                    "sample",
                    "scheme",
                    "ST",
                ]
            ]
            .drop_duplicates(
                subset=["sample"]
            )
        )

        summary = summary.merge(
            mlst_simple,
            on="sample",
            how="left",
        )

    else:

        summary["scheme"] = "NA"
        summary["ST"] = "NA"

    # -----------------------------------------------------------------
    # Serotyping
    # -----------------------------------------------------------------

    if not serotyping.empty:

        serotype_simple = (
            serotyping[
                [
                    "sample",
                    "serotype",
                    "best_match",
                    "confidence",
                    "status",
                ]
            ]
            .drop_duplicates(
                subset=["sample"]
            )
        )

        summary = summary.merge(
            serotype_simple,
            on="sample",
            how="left",
        )

    else:

        summary["serotype"] = "NA"
        summary["best_match"] = "NA"
        summary["confidence"] = "NA"
        summary["status"] = "NA"

    # -----------------------------------------------------------------
    # AMR
    # -----------------------------------------------------------------

    if not amr.empty:

        amr_summary = (
            amr
            .groupby("sample")
            .agg(
                amr_gene_count=(
                    "gene",
                    "nunique",
                ),
                amr_genes=(
                    "gene",
                    lambda x: unique_join(x),
                ),
                amr_classes=(
                    "class",
                    lambda x: unique_join(x),
                ),
                amr_subclasses=(
                    "subclass",
                    lambda x: unique_join(x),
                ),
                amr_methods=(
                    "method",
                    lambda x: unique_join(x),
                ),
            )
            .reset_index()
        )

        summary = summary.merge(
            amr_summary,
            on="sample",
            how="left",
        )

    else:

        summary["amr_gene_count"] = 0
        summary["amr_genes"] = ""
        summary["amr_classes"] = ""
        summary["amr_subclasses"] = ""
        summary["amr_methods"] = ""

    # -----------------------------------------------------------------
    # Virulence
    # -----------------------------------------------------------------

    if not virulence.empty:

        vf_summary = (
            virulence
            .groupby("sample")
            .agg(
                virulence_gene_count=(
                    "gene",
                    "nunique",
                ),
                virulence_genes=(
                    "gene",
                    lambda x: unique_join(x),
                ),
            )
            .reset_index()
        )

        summary = summary.merge(
            vf_summary,
            on="sample",
            how="left",
        )

    else:

        summary["virulence_gene_count"] = 0
        summary["virulence_genes"] = ""

    # -----------------------------------------------------------------
    # Normalize
    # -----------------------------------------------------------------

    text_columns = [
        "scheme",
        "ST",
        "serotype",
        "best_match",
        "confidence",
        "status",
        "amr_genes",
        "amr_classes",
        "amr_subclasses",
        "amr_methods",
        "virulence_genes",
    ]

    for column in text_columns:

        summary[column] = (
            summary[column]
            .fillna("NA")
            .astype(str)
        )

    count_columns = [
        "amr_gene_count",
        "virulence_gene_count",
    ]

    for column in count_columns:

        summary[column] = (
            pd.to_numeric(
                summary[column],
                errors="coerce",
            )
            .fillna(0)
            .astype(int)
        )

    # -----------------------------------------------------------------
    # Human-readable derived fields
    # -----------------------------------------------------------------

    summary["serotype_call_valid"] = (
        summary["status"] == "OK"
    )

    summary["total_surveillance_features"] = (
        summary["amr_gene_count"]
        +
        summary["virulence_gene_count"]
    )

    return summary


# =============================================================================
# Integrated burden matrix
# =============================================================================

def plot_integrated_burden(
    summary: pd.DataFrame,
    outdir: Path,
) -> bool:

    if summary.empty:
        return False

    matrix = summary[
        [
            "sample",
            "amr_gene_count",
            "virulence_gene_count",
        ]
    ].copy()

    matrix = matrix.set_index(
        "sample"
    )

    matrix.columns = [
        "AMR determinants",
        "Virulence determinants",
    ]

    return plot_heatmap(
        matrix,
        "Integrated AMR and Virulence Burden",
        outdir,
        "integrated_amr_virulence_burden",
        cmap="YlOrRd",
        colorbar_label="Number of detected determinants",
    )


# =============================================================================
# AMR class distribution
# =============================================================================

def plot_amr_class_distribution(
    amr: pd.DataFrame,
    outdir: Path,
) -> bool:

    if amr.empty:
        return False

    classes = (
        amr["class"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    classes = classes[
        (classes != "")
        &
        (classes.str.lower() != "na")
    ]

    if classes.empty:
        return False

    counts = (
        classes
        .value_counts()
        .sort_values(ascending=False)
    )

    counts.to_csv(
        outdir / "amr_class_counts.tsv",
        sep="\t",
        header=["detection_count"],
    )

    fig, ax = plt.subplots(
        figsize=(
            max(7, len(counts) * 0.7),
            5,
        )
    )

    counts.plot(
        kind="bar",
        ax=ax,
        edgecolor="black",
    )

    ax.set_title(
        "AMRFinderPlus Antimicrobial Classes",
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel("AMR class")
    ax.set_ylabel("Detection count")
    ax.tick_params(
        axis="x",
        rotation=45,
    )

    save_figure(
        fig,
        outdir,
        "amr_class_distribution",
    )

    return True


# =============================================================================
# AMR method distribution
# =============================================================================

def plot_amr_method_distribution(
    amr: pd.DataFrame,
    outdir: Path,
) -> bool:

    if amr.empty:
        return False

    methods = (
        amr["method"]
        .fillna("")
        .astype(str)
        .str.strip()
    )

    methods = methods[
        methods != ""
    ]

    if methods.empty:
        return False

    counts = (
        methods
        .value_counts()
        .sort_values(ascending=False)
    )

    counts.to_csv(
        outdir / "amr_method_counts.tsv",
        sep="\t",
        header=["detection_count"],
    )

    fig, ax = plt.subplots(
        figsize=(
            max(7, len(counts) * 0.7),
            5,
        )
    )

    counts.plot(
        kind="bar",
        ax=ax,
        edgecolor="black",
    )

    ax.set_title(
        "AMRFinderPlus Detection Method Distribution",
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel("Detection method")
    ax.set_ylabel("Detection count")
    ax.tick_params(
        axis="x",
        rotation=45,
    )

    save_figure(
        fig,
        outdir,
        "amr_method_distribution",
    )

    return True


# =============================================================================
# Phylogenetic visualization
# =============================================================================

def clean_tree_sample_name(
    name: str,
) -> str:

    value = str(name)

    for suffix in (
        ".fna",
        ".fasta",
        ".fa",
        ".fsa",
        ".fna.gz",
        ".fasta.gz",
    ):
        value = value.replace(
            suffix,
            "",
        )

    return value


def render_tree(
    tree_path: Path,
    annotations: Dict[str, str],
    title: str,
    outdir: Path,
    basename: str,
) -> bool:

    if not HAS_ETE3:
        return False

    if not tree_path.exists():
        return False

    if tree_path.stat().st_size == 0:
        return False

    try:

        tree = Tree(
            str(tree_path),
            format=1,
        )

    except Exception as exc:

        log(
            f"WARNING: unable to parse tree "
            f"for {basename}: {exc}"
        )

        return False

    for leaf in tree.iter_leaves():

        sample = clean_tree_sample_name(
            leaf.name
        )

        style = NodeStyle()

        style["size"] = 7

        leaf.set_style(style)

        if sample in annotations:

            label = annotations[
                sample
            ]

            if label:

                leaf.add_face(
                    TextFace(
                        f" [{label}]",
                        fsize=8,
                    ),
                    column=0,
                    position="branch-right",
                )

    tree_style = TreeStyle()

    tree_style.show_leaf_name = True
    tree_style.mode = "r"

    tree_style.title.add_face(
        TextFace(
            title,
            fsize=14,
        ),
        column=0,
    )

    for extension in (
        "png",
        "svg",
        "pdf",
    ):

        tree.render(
            str(
                outdir
                / f"{basename}.{extension}"
            ),
            tree_style=tree_style,
            w=180,
            units="mm",
        )

    log(
        f"Saved {basename}.{{png,svg,pdf}}"
    )

    return True


# =============================================================================
# Dashboard
# =============================================================================

def plot_dashboard(
    summary: pd.DataFrame,
    outdir: Path,
) -> bool:

    if summary.empty:
        return False

    summary.to_csv(
        outdir / "sample_summary_dashboard.tsv",
        sep="\t",
        index=False,
    )

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(16, 10),
    )

    # AMR burden
    summary.plot(
        x="sample",
        y="amr_gene_count",
        kind="bar",
        legend=False,
        edgecolor="black",
        ax=axes[0, 0],
    )

    axes[0, 0].set_title(
        "AMR determinants per sample"
    )

    axes[0, 0].set_xlabel("")
    axes[0, 0].set_ylabel(
        "AMRFinderPlus detections"
    )

    axes[0, 0].tick_params(
        axis="x",
        rotation=45,
    )

    # Virulence burden
    summary.plot(
        x="sample",
        y="virulence_gene_count",
        kind="bar",
        legend=False,
        edgecolor="black",
        ax=axes[0, 1],
    )

    axes[0, 1].set_title(
        "Virulence determinants per sample"
    )

    axes[0, 1].set_xlabel("")
    axes[0, 1].set_ylabel(
        "Virulence detections"
    )

    axes[0, 1].tick_params(
        axis="x",
        rotation=45,
    )

    # ST distribution
    st_counts = (
        summary["ST"]
        .value_counts()
    )

    st_counts.plot(
        kind="bar",
        edgecolor="black",
        ax=axes[1, 0],
    )

    axes[1, 0].set_title(
        "MLST sequence type distribution"
    )

    axes[1, 0].set_xlabel("Sequence type")
    axes[1, 0].set_ylabel(
        "Number of samples"
    )

    axes[1, 0].tick_params(
        axis="x",
        rotation=45,
    )

    # Serotype distribution
    serotype_counts = (
        summary[
            summary["serotype_call_valid"]
        ]["serotype"]
        .value_counts()
    )

    if serotype_counts.empty:

        axes[1, 1].text(
            0.5,
            0.5,
            "No successful GBS-SBG calls",
            ha="center",
            va="center",
            transform=axes[1, 1].transAxes,
        )

    else:

        serotype_counts.plot(
            kind="bar",
            edgecolor="black",
            ax=axes[1, 1],
        )

    axes[1, 1].set_title(
        "GBS capsular serotype distribution"
    )

    axes[1, 1].set_xlabel("Serotype")
    axes[1, 1].set_ylabel(
        "Number of samples"
    )

    axes[1, 1].tick_params(
        axis="x",
        rotation=45,
    )

    fig.suptitle(
        "GBS-Sentinel Surveillance Dashboard",
        fontsize=17,
        fontweight="bold",
    )

    fig.tight_layout()

    save_figure(
        fig,
        outdir,
        "sample_summary_dashboard",
    )

    return True


# =============================================================================
# Main
# =============================================================================

def main() -> None:

    parser = argparse.ArgumentParser(
        description=(
            "GBS-Sentinel integrated surveillance "
            "figure generation"
        )
    )

    parser.add_argument(
        "--amr-dir",
        required=True,
        help="Directory containing AMRFinderPlus *_amrfinder.tsv files",
    )

    parser.add_argument(
        "--vf-dir",
        required=True,
        help="Directory containing *_vf.tsv files",
    )

    parser.add_argument(
        "--mlst-dir",
        required=True,
        help="Directory containing *_mlst.tsv files",
    )

    parser.add_argument(
        "--serotype-dir",
        required=True,
        help="Directory containing *_serotype.tsv files",
    )

    parser.add_argument(
        "--tree",
        required=True,
        help="Core-genome Newick tree",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )

    args = parser.parse_args()

    amr_dir = Path(args.amr_dir)
    vf_dir = Path(args.vf_dir)
    mlst_dir = Path(args.mlst_dir)
    serotype_dir = Path(args.serotype_dir)
    tree_path = Path(args.tree)
    outdir = Path(args.outdir)

    reports_dir = (
        outdir / "Reports"
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    reports_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    log(
        "=================================================="
    )

    log(
        "GBS-Sentinel Integrated Visualization"
    )

    log(
        "=================================================="
    )

    log(
        f"AMRFinderPlus directory: {amr_dir}"
    )

    log(
        f"Virulence directory:    {vf_dir}"
    )

    log(
        f"MLST directory:         {mlst_dir}"
    )

    log(
        f"Serotyping directory:   {serotype_dir}"
    )

    log(
        f"Phylogenetic tree:      {tree_path}"
    )

    # -------------------------------------------------------------------------
    # Load actual pipeline results
    # -------------------------------------------------------------------------

    amr = load_amrfinder(
        amr_dir
    )

    virulence = load_virulence(
        vf_dir
    )

    mlst = load_mlst(
        mlst_dir
    )

    serotyping = load_serotyping(
        serotype_dir
    )

    log(
        f"Loaded AMRFinderPlus hits: {len(amr)}"
    )

    log(
        f"Loaded virulence hits:     {len(virulence)}"
    )

    log(
        f"Loaded MLST records:       {len(mlst)}"
    )

    log(
        f"Loaded serotype records:   {len(serotyping)}"
    )

    # -------------------------------------------------------------------------
    # Build matrices
    # -------------------------------------------------------------------------

    amr_matrix = build_gene_matrix(
        amr
    )

    vf_matrix = build_gene_matrix(
        virulence
    )

    amr_class_matrix = build_amr_class_matrix(
        amr
    )

    # -------------------------------------------------------------------------
    # Canonical surveillance table
    # -------------------------------------------------------------------------

    summary = build_surveillance_summary(
        mlst=mlst,
        serotyping=serotyping,
        amr=amr,
        virulence=virulence,
    )

    if summary.empty:

        raise RuntimeError(
            "No sample-level surveillance data "
            "were available for integration."
        )

    summary_path = (
        outdir
        / "surveillance_summary.tsv"
    )

    summary.to_csv(
        summary_path,
        sep="\t",
        index=False,
    )

    log(
        f"Canonical surveillance table: {summary_path}"
    )

    # -------------------------------------------------------------------------
    # Save individual standardized data tables
    # -------------------------------------------------------------------------

    amr.to_csv(
        outdir / "amrfinder_summary.tsv",
        sep="\t",
        index=False,
    )

    virulence.to_csv(
        outdir / "virulence_summary.tsv",
        sep="\t",
        index=False,
    )

    mlst.to_csv(
        outdir / "mlst_summary.tsv",
        sep="\t",
        index=False,
    )

    serotyping.to_csv(
        outdir / "serotype_summary.tsv",
        sep="\t",
        index=False,
    )

    # -------------------------------------------------------------------------
    # AMR figures
    # -------------------------------------------------------------------------

    if amr_matrix.empty:

        write_skip_report(
            outdir,
            "amr_heatmap",
            "No AMRFinderPlus gene detections were available.",
        )

    else:

        plot_heatmap(
            amr_matrix,
            "AMRFinderPlus Resistance Determinant Presence/Absence",
            outdir,
            "amr_heatmap",
            cmap="YlOrRd",
        )

        plot_frequency(
            amr_matrix,
            "AMRFinderPlus Resistance Determinant Frequency",
            outdir,
            "amr_frequency",
        )

        plot_sample_counts(
            amr_matrix,
            "AMRFinderPlus Resistance Determinants per Sample",
            outdir,
            "resistance_gene_counts",
        )

        plot_cooccurrence(
            amr_matrix,
            "AMRFinderPlus Resistance Determinant Co-occurrence",
            outdir,
            "amr_cooccurrence",
        )


    plot_amr_class_distribution(
        amr,
        outdir,
    )

    plot_amr_method_distribution(
        amr,
        outdir,
    )

    if not amr_class_matrix.empty:

        plot_heatmap(
            amr_class_matrix,
            "AMRFinderPlus AMR Class Distribution by Sample",
            outdir,
            "amr_class_heatmap",
            cmap="YlGnBu",
        )

    # -------------------------------------------------------------------------
    # Virulence figures
    # -------------------------------------------------------------------------

    if vf_matrix.empty:

        write_skip_report(
            outdir,
            "virulence_heatmap",
            "No virulence-factor detections were available.",
        )

    else:

        plot_heatmap(
            vf_matrix,
            "Virulence Factor Presence/Absence",
            outdir,
            "virulence_heatmap",
            cmap="PuBu",
        )

        plot_frequency(
            vf_matrix,
            "Virulence Factor Frequency",
            outdir,
            "virulence_frequency",
        )

        plot_sample_counts(
            vf_matrix,
            "Virulence Determinants per Sample",
            outdir,
            "virulence_gene_counts",
        )

        plot_cooccurrence(
            vf_matrix,
            "Virulence Factor Co-occurrence",
            outdir,
            "virulence_cooccurrence",
        )

    # -------------------------------------------------------------------------
    # Combined AMR and virulence
    # -------------------------------------------------------------------------

    pieces = []

    if not amr_matrix.empty:
        pieces.append(
            amr_matrix.add_prefix(
                "AMR:"
            )
        )

    if not vf_matrix.empty:
        pieces.append(
            vf_matrix.add_prefix(
                "VF:"
            )
        )

    if pieces:

        combined = pd.concat(
            pieces,
            axis=1,
        ).fillna(0)

        combined = combined.astype(int)

        plot_heatmap(
            combined,
            "Integrated AMR and Virulence Determinants",
            outdir,
            "combined_amr_virulence_heatmap",
            cmap="RdYlBu_r",
        )

    # -------------------------------------------------------------------------
    # MLST
    # -------------------------------------------------------------------------

    plot_category_distribution(
        summary,
        "ST",
        "MLST Sequence Type Distribution",
        outdir,
        "mlst_distribution",
    )

    # -------------------------------------------------------------------------
    # Serotyping
    # -------------------------------------------------------------------------

    successful_serotypes = summary[
        summary["serotype_call_valid"]
    ].copy()

    if successful_serotypes.empty:

        write_skip_report(
            outdir,
            "serotype_distribution",
            "No successful GBS-SBG serotype calls were available.",
        )

    else:

        plot_category_distribution(
            successful_serotypes,
            "serotype",
            "GBS Capsular Serotype Distribution",
            outdir,
            "serotype_distribution",
        )

    # -------------------------------------------------------------------------
    # ST × Serotype
    # -------------------------------------------------------------------------

    if not successful_serotypes.empty:

        plot_cross_tab(
            successful_serotypes,
            "ST",
            "serotype",
            "MLST Sequence Type × GBS Capsular Serotype",
            outdir,
            "st_vs_serotype",
        )

    # -------------------------------------------------------------------------
    # ST × AMR burden
    # -------------------------------------------------------------------------

    st_amr = summary[
        ["ST", "amr_gene_count"]
    ].copy()

    if not st_amr.empty:

        st_amr_matrix = (
            st_amr
            .groupby("ST")
            ["amr_gene_count"]
            .agg(
                [
                    "count",
                    "mean",
                    "median",
                    "max",
                ]
            )
        )

        st_amr_matrix.to_csv(
            outdir / "st_amr_burden_summary.tsv",
            sep="\t",
        )

        plot_heatmap(
            st_amr_matrix,
            "AMR Burden Statistics by MLST Sequence Type",
            outdir,
            "st_amr_burden",
            cmap="YlOrRd",
            colorbar_label="AMR burden statistic",
        )

    # -------------------------------------------------------------------------
    # Serotype × AMR burden
    # -------------------------------------------------------------------------

    if not successful_serotypes.empty:

        serotype_amr = (
            successful_serotypes
            .groupby("serotype")
            ["amr_gene_count"]
            .agg(
                [
                    "count",
                    "mean",
                    "median",
                    "max",
                ]
            )
        )

        serotype_amr.to_csv(
            outdir / "serotype_amr_burden_summary.tsv",
            sep="\t",
        )

        plot_heatmap(
            serotype_amr,
            "AMR Burden Statistics by GBS Serotype",
            outdir,
            "serotype_amr_burden",
            cmap="YlOrRd",
            colorbar_label="AMR burden statistic",
        )

    # -------------------------------------------------------------------------
    # ST × Virulence burden
    # -------------------------------------------------------------------------

    st_vf = (
        summary
        .groupby("ST")
        ["virulence_gene_count"]
        .agg(
            [
                "count",
                "mean",
                "median",
                "max",
            ]
        )
    )

    if not st_vf.empty:

        st_vf.to_csv(
            outdir / "st_virulence_burden_summary.tsv",
            sep="\t",
        )

        plot_heatmap(
            st_vf,
            "Virulence Burden Statistics by MLST Sequence Type",
            outdir,
            "st_virulence_burden",
            cmap="PuBu",
            colorbar_label="Virulence burden statistic",
        )

    # -------------------------------------------------------------------------
    # Serotype × Virulence burden
    # -------------------------------------------------------------------------

    if not successful_serotypes.empty:

        serotype_vf = (
            successful_serotypes
            .groupby("serotype")
            ["virulence_gene_count"]
            .agg(
                [
                    "count",
                    "mean",
                    "median",
                    "max",
                ]
            )
        )

        serotype_vf.to_csv(
            outdir / "serotype_virulence_burden_summary.tsv",
            sep="\t",
        )

        plot_heatmap(
            serotype_vf,
            "Virulence Burden Statistics by GBS Serotype",
            outdir,
            "serotype_virulence_burden",
            cmap="PuBu",
            colorbar_label="Virulence burden statistic",
        )

    # -------------------------------------------------------------------------
    # Integrated AMR/Virulence burden
    # -------------------------------------------------------------------------

    plot_integrated_burden(
        summary,
        outdir,
    )

    # -------------------------------------------------------------------------
    # AMR burden versus virulence burden
    # -------------------------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(8, 6)
    )

    ax.scatter(
        summary["amr_gene_count"],
        summary["virulence_gene_count"],
        s=70,
        edgecolor="black",
    )

    for _, row in summary.iterrows():

        ax.annotate(
            row["sample"],
            (
                row["amr_gene_count"],
                row["virulence_gene_count"],
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8,
        )

    ax.set_title(
        "AMR Burden versus Virulence Burden",
        fontsize=15,
        fontweight="bold",
    )

    ax.set_xlabel(
        "AMRFinderPlus resistance determinants"
    )

    ax.set_ylabel(
        "Virulence determinants"
    )

    save_figure(
        fig,
        outdir,
        "amr_vs_virulence_burden",
    )

    # -------------------------------------------------------------------------
    # Surveillance dashboard
    # -------------------------------------------------------------------------

    plot_dashboard(
        summary,
        outdir,
    )

    # -------------------------------------------------------------------------
    # Phylogenetic annotations
    # -------------------------------------------------------------------------

    if tree_path.exists():

        st_annotations = {
            row["sample"]:
                f"ST:{row['ST']}"
            for _, row in summary.iterrows()
        }

        serotype_annotations = {
            row["sample"]:
                f"Serotype:{row['serotype']}"
            for _, row in summary.iterrows()
        }

        amr_annotations = {
            row["sample"]:
                f"AMR:{row['amr_gene_count']}"
            for _, row in summary.iterrows()
        }

        virulence_annotations = {
            row["sample"]:
                f"VF:{row['virulence_gene_count']}"
            for _, row in summary.iterrows()
        }

        integrated_annotations = {
            row["sample"]:
                (
                    f"ST:{row['ST']} | "
                    f"Serotype:{row['serotype']} | "
                    f"AMR:{row['amr_gene_count']} | "
                    f"VF:{row['virulence_gene_count']}"
                )
            for _, row in summary.iterrows()
        }

        render_tree(
            tree_path,
            st_annotations,
            "Core-Genome Phylogeny with MLST Sequence Types",
            outdir,
            "tree_with_st_labels",
        )

        render_tree(
            tree_path,
            serotype_annotations,
            "Core-Genome Phylogeny with GBS Capsular Serotypes",
            outdir,
            "tree_with_serotype_annotation",
        )

        render_tree(
            tree_path,
            amr_annotations,
            "Core-Genome Phylogeny with AMR Burden",
            outdir,
            "tree_with_amr_annotation",
        )

        render_tree(
            tree_path,
            virulence_annotations,
            "Core-Genome Phylogeny with Virulence Burden",
            outdir,
            "tree_with_virulence_annotation",
        )

        render_tree(
            tree_path,
            integrated_annotations,
            "Integrated GBS-Sentinel Surveillance Phylogeny",
            outdir,
            "tree_integrated_surveillance",
        )

    else:

        write_skip_report(
            outdir,
            "phylogeny_visualizations",
            "Core-genome Newick tree not found.",
        )

    # -------------------------------------------------------------------------
    # Final report
    # -------------------------------------------------------------------------

    figure_files = sorted(
        path.name
        for path in outdir.iterdir()
        if path.suffix.lower()
        in {
            ".png",
            ".svg",
            ".pdf",
        }
    )

    successful_serotype_count = int(
        summary[
            "serotype_call_valid"
        ].sum()
    )

    report = reports_dir / "visualization_summary.txt"

    report.write_text(
        "GBS-Sentinel Integrated Visualization Report\n"
        "============================================\n\n"
        f"Generated: {datetime.now().isoformat()}\n"
        f"Samples integrated: {len(summary)}\n"
        f"AMRFinderPlus detections: {len(amr)}\n"
        f"Virulence detections: {len(virulence)}\n"
        f"MLST records: {len(mlst)}\n"
        f"Serotyping records: {len(serotyping)}\n"
        f"Successful serotype calls: {successful_serotype_count}\n\n"
        "Integrated biological dimensions:\n"
        "  - MLST sequence type\n"
        "  - GBS capsular serotype\n"
        "  - AMRFinderPlus resistance determinants\n"
        "  - AMR class/subclass information\n"
        "  - Virulence determinants\n"
        "  - Core-genome phylogeny\n\n"
        "Canonical table:\n"
        "  - surveillance_summary.tsv\n\n"
        "Generated figures:\n"
        + "\n".join(
            f"  - {name}"
            for name in figure_files
        )
        + "\n"
    )

    log(
        "=================================================="
    )

    log(
        "GBS-Sentinel visualization completed"
    )

    log(
        f"Samples integrated: {len(summary)}"
    )

    log(
        f"Successful serotype calls: "
        f"{successful_serotype_count}"
    )

    log(
        f"Figures generated: {len(figure_files)}"
    )

    log(
        f"Canonical table: {summary_path}"
    )

    log(
        "=================================================="
    )


if __name__ == "__main__":
    main()

    