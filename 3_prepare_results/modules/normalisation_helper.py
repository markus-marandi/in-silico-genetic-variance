from collections import defaultdict
import gzip
import re
from pathlib import Path
from typing import Any, Iterable

import pandas as pd
import polars as pl
from polars.exceptions import ComputeError


def strip_ensembl_version(gene_id: str) -> str:
    """strip ensembl version suffix.

    args:
        gene_id (str): ensembl id with or without version.

    returns:
        str: ensembl id without version.
    """
    if not isinstance(gene_id, str):
        return gene_id
    return gene_id.split('.')[0]


def normalize_gene_id_column(df: pd.DataFrame, col: str) -> pd.Series:
    """return version-stripped ensembl ids from a dataframe column.

    args:
        df (pd.DataFrame): input frame.
        col (str): column with ensembl ids.

    returns:
        pd.Series: normalized ids.
    """
    return df[col].astype(str).map(strip_ensembl_version)


def _parse_gtf_attributes(attr: str) -> dict[str, str]:
    # minimal robust parser for gtf attributes
    out: dict[str, str] = {}
    for item in attr.strip().strip(';').split(';'):
        item = item.strip()
        if not item:
            continue
        if ' ' not in item:
            continue
        key, val = item.split(' ', 1)
        val = val.strip().strip('"')
        out[key] = val
    return out


def parse_gtf_genes(gtf_path: Path) -> pd.DataFrame:
    """parse gtf gene records with gene_id, gene_name, gene_type, coords, tss, length.

    args:
        gtf_path (Path): path to gencode gtf.

    returns:
        pd.DataFrame: gene metadata indexed by gene_id.
    """
    rows = []
    with gtf_path.open('rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature != 'gene':
                continue
            a = _parse_gtf_attributes(attrs)
            gene_id = strip_ensembl_version(a.get('gene_id', ''))
            gene_name = a.get('gene_name', None)
            gene_type = a.get('gene_type', a.get('gene_biotype', None))
            start_i = int(start)
            end_i = int(end)
            if strand == '+':
                tss = start_i
            else:
                tss = end_i
            rows.append(
                {
                    'gene_id': gene_id,
                    'gene_symbol': gene_name,
                    'gene_type': gene_type,
                    'chrom': chrom,
                    'start': start_i,
                    'end': end_i,
                    'strand': strand,
                    'tss': tss,
                    'gene_length': end_i - start_i + 1,
                }
            )
    df = pd.DataFrame(rows)
    return df


def _normalize_label(s: str) -> str:
    # normalize tissue label strings for matching
    return re.sub(r'[^a-z0-9]+', '_', s.lower()).strip('_')


def find_muscle_tpm_column(columns: Iterable[str]) -> str:
    """find muscle skeletal column name in gtex median tpm.

    args:
        columns (Iterable[str]): dataframe columns from gct file.

    returns:
        str: column name for muscle skeletal.
    """
    candidates = [c for c in columns if c not in {'Name', 'Description', 'gene_id', 'gene_name'}]
    # try exact known labels first
    for exact in ['Muscle_Skeletal', 'Muscle - Skeletal']:
        if exact in candidates:
            return exact
    # normalized contains both words
    for c in candidates:
        norm = _normalize_label(c)
        if 'muscle' in norm and 'skeletal' in norm:
            return c
    raise ValueError('muscle skeletal column not found in tpm gct')


def load_tpm_muscle(tpm_path: Path) -> pd.DataFrame:
    """load median tpm for muscle skeletal.

    args:
        tpm_path (Path): path to gct with median tpm.

    returns:
        pd.DataFrame: columns: gene_id, gene_symbol(optional), tpm_muscle.
    """
    gct = pd.read_csv(tpm_path, sep='\t', skiprows=2, low_memory=False)
    gene_col = 'Name' if 'Name' in gct.columns else (
        'gene_id' if 'gene_id' in gct.columns else None
    )
    if gene_col is None:
        raise ValueError('gene id column not found in tpm gct')
    muscle_col = find_muscle_tpm_column(gct.columns)
    out = gct[[gene_col, 'Description', muscle_col]].rename(
        columns={gene_col: 'gene_id', 'Description': 'gene_symbol', muscle_col: 'tpm_muscle'}
    )
    out['gene_id'] = out['gene_id'].map(strip_ensembl_version)
    return out


_ag_vid_re = re.compile(r'^(?:chr)?(?P<chrom>[0-9XYM]+):(?P<pos>\d+):(?P<ref>[ACGTN]+)>(?P<alt>[ACGTN]+)$')
_gtex_vid_re = re.compile(r'^(?P<chrom>chr[0-9XYM]+)_(?P<pos>\d+)_(?P<ref>[ACGTN]+)_(?P<alt>[ACGTN]+)(?:_b\d+)?$')


def ag_variant_to_canonical(variant_id: Any) -> str | None:
    """convert ag variant string to canonical chr_pos_ref_alt.

    args:
        variant_id (Any): input id like chr20:52181783:T>C.

    returns:
        str | None: canonical id or None when unparseable.
    """
    if not isinstance(variant_id, str):
        return None
    m = _ag_vid_re.match(variant_id)
    if not m:
        return None
    chrom = m.group('chrom')
    if not chrom.startswith('chr'):
        chrom = f'chr{chrom}'
    return f"{chrom}_{m.group('pos')}_{m.group('ref')}_{m.group('alt')}"


def gtex_variant_to_canonical(variant_id: Any) -> str | None:
    """convert gtex variant string to canonical chr_pos_ref_alt.

    args:
        variant_id (Any): input id like chr20_52181783_T_C_b38.

    returns:
        str | None: canonical id or None when unparseable.
    """
    if not isinstance(variant_id, str):
        return None
    m = _gtex_vid_re.match(variant_id)
    if not m:
        return None
    return f"{m.group('chrom')}_{m.group('pos')}_{m.group('ref')}_{m.group('alt')}"


_METHOD_MAP: dict[str, str] = {
    # normalize to the requested friendly labels
    'genemasklfcscorer': 'gene_exonmask_delta_log2',
    'genemask_lfc': 'gene_exonmask_delta_log2',
    'gene_exonmask_delta_log2': 'gene_exonmask_delta_log2',
    'center2.001kb_diff_mean': 'center2.001kb_diff_mean',
    'center10.001kb_diff_mean': 'center10.001kb_diff_mean',
    'center100.001kb_diff_mean': 'center100.001kb_diff_mean',
    'center2.001kb_l2_log1p': 'center2.001kb_l2_log1p',
    'center10.001kb_l2_log1p': 'center10.001kb_l2_log1p',
    'center100.001kb_l2_log1p': 'center100.001kb_l2_log1p',
}


def friendly_method_name(name: Any) -> str:
    """map raw scorer/method names to clean labels.

    args:
        name (Any): raw method name.

    returns:
        str: cleaned method label.
    """
    if not isinstance(name, str):
        return 'unknown'
    key = name.strip().lower()
    return _METHOD_MAP.get(key, name)


def choose_column(df: pd.DataFrame, candidates: list[str]) -> str:
    """select first existing column from candidates.

    args:
        df (pd.DataFrame): input frame.
        candidates (list[str]): ordered column names to try.

    returns:
        str: chosen column name.
    """
    for c in candidates:
        if c in df.columns:
            return c
    # try case-insensitive
    lower = {c.lower(): c for c in df.columns}
    for c in candidates:
        if c.lower() in lower:
            return lower[c.lower()]
    raise KeyError(f'none of the columns found: {candidates}')


# center-scorer label normalization to strict canonical names used in plots
_cm_re = re.compile(r'CenterMaskScorer.*width=(\d+).*aggregation_type=([A-Z0-9_]+)')


def _normalize_center_label(s: object) -> str:
    """normalize raw center-scorer names to canonical labels.

    args:
        s (object): raw label (e.g., CenterMaskScorer(...)).

    returns:
        str: canonical label like 'center2.001kb_diff_mean'.
    """
    if not isinstance(s, str):
        return 'unknown'
    if s in {
        'gene_exonmask_delta_log2',
        'center2.001kb_diff_mean', 'center10.001kb_diff_mean', 'center100.001kb_diff_mean',
        'center2.001kb_l2_log1p', 'center10.001kb_l2_log1p', 'center100.001kb_l2_log1p',
    }:
        return s
    m = _cm_re.search(s)
    if m:
        w = int(m.group(1)); agg = m.group(2)
        width_map = {2001: 'center2.001kb', 10001: 'center10.001kb', 100001: 'center100.001kb'}
        agg_map = {'DIFF_MEAN': 'diff_mean', 'L2_DIFF_LOG1P': 'l2_log1p'}
        return f"{width_map.get(w, f'center{w}')}_{agg_map.get(agg, agg.lower())}"
    if 'GeneMaskLFCScorer' in s:
        return 'gene_exonmask_delta_log2'
    return s

def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """merge overlapping/adjacent intervals."""
    if not intervals:
        return []
    intervals = sorted(intervals, key=lambda x: x[0])
    merged: list[tuple[int, int]] = []
    cur_start, cur_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_end + 1:  # overlapping or directly adjacent
            cur_end = max(cur_end, e)
        else:
            merged.append((cur_start, cur_end))
            cur_start, cur_end = s, e
    merged.append((cur_start, cur_end))
    return merged


def _interval_total_length(intervals: list[tuple[int, int]]) -> int:
    """total length of merged intervals in bp."""
    if not intervals:
        return 0
    merged = _merge_intervals(intervals)
    return sum(e - s + 1 for s, e in merged)


def parse_gtf_exon_cds_lengths(gtf_path: Path) -> pd.DataFrame:
    """compute exonic and coding (CDS) lengths per gene from a GTF.

    args:
        gtf_path (Path): path to gencode gtf.

    returns:
        pd.DataFrame: columns: gene_id, exonic_length, coding_length.
    """
    exon_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cds_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)

    with gtf_path.open('rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature not in {'exon', 'CDS'}:
                continue
            a = _parse_gtf_attributes(attrs)
            raw_gene_id = a.get('gene_id', '')
            if not raw_gene_id:
                continue
            gene_id = strip_ensembl_version(raw_gene_id)
            s = int(start)
            e = int(end)
            if feature == 'exon':
                exon_intervals[gene_id].append((s, e))
            elif feature == 'CDS':
                cds_intervals[gene_id].append((s, e))

    rows = []
    all_gene_ids = set(exon_intervals) | set(cds_intervals)
    for gid in sorted(all_gene_ids):
        exonic_length = _interval_total_length(exon_intervals.get(gid, []))
        coding_length = _interval_total_length(cds_intervals.get(gid, []))
        rows.append(
            {
                'gene_id': gid,
                'exonic_length': exonic_length,
                'coding_length': coding_length,
            }
        )
    return pd.DataFrame(rows)


def parse_gtf_gene_with_lengths(gtf_path: Path) -> pd.DataFrame:
    """combine gene metadata (including genomic length) with exonic and coding lengths.

    args:
        gtf_path (Path): path to gencode gtf.

    returns:
        pd.DataFrame: one row per gene_id with:
            gene_id, gene_symbol, gene_type, chrom, start, end, strand, tss,
            genomic_length, exonic_length, coding_length,
            intronic_length, utr_length
    """
    genes = parse_gtf_genes(gtf_path).copy()
    genes = genes.rename(columns={'gene_length': 'genomic_length'})

    exon_cds = parse_gtf_exon_cds_lengths(gtf_path)

    out = genes.merge(exon_cds, on='gene_id', how='left')
    out[['exonic_length', 'coding_length']] = (
        out[['exonic_length', 'coding_length']].fillna(0).astype(int)
    )

    # intragenic non-exonic sequence
    out['intronic_length'] = (
        (out['genomic_length'] - out['exonic_length'])
        .clip(lower=0)
        .astype(int)
    )

    # exonic non-coding sequence (5'UTR + 3'UTR)
    out['utr_length'] = (
        (out['exonic_length'] - out['coding_length'])
        .clip(lower=0)
        .astype(int)
    )

    return out



def load_mane_select(mane_path: Path) -> pd.DataFrame:
    #load MANE Select mapping gene_id -> manu transcript id."""
    mane = pd.read_csv(mane_path, sep='\t', dtype=str)
    if 'MANE_status' in mane.columns:
        mane = mane[mane['MANE_status'] == 'MANE Select'].copy()
    if 'Ensembl_Gene' not in mane.columns or 'Ensembl_nuc' not in mane.columns:
        raise ValueError('MANE summary missing Ensembl_Gene/Ensembl_nuc')
    mane['gene_id'] = mane['Ensembl_Gene'].map(strip_ensembl_version)
    mane = mane.rename(columns={'Ensembl_nuc': 'mane_transcript_id'})
    return mane[['gene_id', 'mane_transcript_id']].drop_duplicates('gene_id')


def parse_gtf_utr5_lengths(gtf_path: Path, mane: pd.DataFrame) -> pd.DataFrame:
    """compute 5'UTR length per gene using MANE Select transcripts."""
    if not {'gene_id', 'mane_transcript_id'} <= set(mane.columns):
        raise ValueError('mane must have gene_id and mane_transcript_id')

    gene_to_tx = dict(zip(mane['gene_id'], mane['mane_transcript_id']))
    tx_to_gene = {tx: gid for gid, tx in gene_to_tx.items()}
    mane_tx = set(tx_to_gene.keys())

    utr_by_tx: dict[str, list[tuple[int, int]]] = defaultdict(list)
    cds_by_tx: dict[str, list[tuple[int, int]]] = defaultdict(list)
    strand_by_tx: dict[str, str] = {}
    gene_by_tx: dict[str, str] = {}

    with gtf_path.open('rt') as fh:
        for line in fh:
            if not line or line.startswith('#'):
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) != 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = parts
            if feature not in {'UTR', 'CDS'}:
                continue
            a = _parse_gtf_attributes(attrs)
            tx_id = a.get('transcript_id')
            if not tx_id or tx_id not in mane_tx:
                continue
            raw_gid = a.get('gene_id', '')
            if not raw_gid:
                continue
            gid = strip_ensembl_version(raw_gid)

            s = int(start)
            e = int(end)
            strand_by_tx.setdefault(tx_id, strand)
            gene_by_tx.setdefault(tx_id, gid)
            if feature == 'UTR':
                utr_by_tx[tx_id].append((s, e))
            else:
                cds_by_tx[tx_id].append((s, e))

    rows: list[dict[str, int | str]] = []
    for tx_id, utr_ints in utr_by_tx.items():
        cds_ints = cds_by_tx.get(tx_id, [])
        if not cds_ints or not utr_ints:
            continue
        gid = gene_by_tx.get(tx_id)
        strand = strand_by_tx.get(tx_id, '+')
        cds_min = min(s for s, _ in cds_ints)
        cds_max = max(e for _, e in cds_ints)

        if strand == '+':
            utr5_ints = [(s, e) for s, e in utr_ints if e < cds_min]
        else:
            utr5_ints = [(s, e) for s, e in utr_ints if s > cds_max]

        utr5_len = _interval_total_length(utr5_ints)
        rows.append({'gene_id': gid, 'transcript_id': tx_id, 'utr5_length': utr5_len})

    if not rows:
        return pd.DataFrame(columns=['gene_id', 'utr5_length'])

    df = pd.DataFrame(rows)
    gene_utr5 = (
        df.groupby('gene_id', as_index=False)['utr5_length']
          .max()
          .astype({'utr5_length': int})
    )
    return gene_utr5


def load_vgh_metrics(metrics_path: Path) -> pd.DataFrame:
    """load vgh gene-level metrics."""
    metrics = [
        'ncRVIS',
        'loeuf_score',
        'ncGERP',
        'RVIS_score',
        'ncCADD',
        'pHaplo',
        'pTriplo',
        'Episcore',
        'pLI',
        'median_tpm',
        'num_enh',
        'num_super_enh',
        'tau',
        'vg_eqtl',
    ]

    read_kwargs = dict(
        separator='\t',
        has_header=True,
        null_values=['NA', 'NaN', ''],
    )

    try:
        lf = pl.read_csv(metrics_path, ignore_errors=False, **read_kwargs)
    except ComputeError as exc:
        # detect ragged lines and report explicitly
        expected_cols = None
        ragged_count = 0
        ragged_examples: list[str] = []
        opener = gzip.open if metrics_path.suffix in {'.gz', '.bgz'} else open

        with opener(metrics_path, 'rt', encoding='utf-8', errors='replace') as fh:
            for idx, line in enumerate(fh, start=1):
                if idx == 1:
                    expected_cols = len(line.rstrip('\n').split('\t'))
                    continue
                if not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if expected_cols is not None and len(parts) != expected_cols:
                    ragged_count += 1
                    if len(ragged_examples) < 5:
                        preview = '\t'.join(parts[:5])
                        ragged_examples.append(
                            f"line {idx}: fields={len(parts)} expected={expected_cols} preview={preview}"
                        )

        if expected_cols is None:
            expected_cols = 0

        print("warning: vgh parse error ->", exc)
        print(
            f"warning: vgh ragged rows dropped={ragged_count}, expected_cols={expected_cols}, "
            "examples (first 5):"
        )
        for example in ragged_examples:
            print(f"warning:   {example}")

        lf = pl.read_csv(
            metrics_path,
            ignore_errors=True,
            truncate_ragged_lines=True,
            **read_kwargs,
        )

    if 'vgh' not in lf.columns:
        raise ValueError('vgh column not found in VGH metrics file')

    lf = lf.rename({'vgh': 'gene_id'}).with_columns(
        pl.col('gene_id').cast(pl.Utf8).str.split('.').list.get(0)
    )

    cast_cols = [c for c in metrics if c in lf.columns]
    if cast_cols:
        lf = lf.with_columns(pl.col(cast_cols).cast(pl.Float64, strict=False))

    present_metrics = [c for c in metrics if c in lf.columns]
    cols = ['gene_id'] + present_metrics
    return lf.select(cols).unique('gene_id').to_pandas()