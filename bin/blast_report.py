#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import base64
from io import BytesIO
import datetime
import os
import sys
import argparse
from jinja2 import Environment, FileSystemLoader, select_autoescape

def parser_args(args=None):
    Description = "Generate HTML reports from filtered BLAST results."
    Epilog = """Example usage:
    python blast_report.py --blast_file sample.filter.blastn.txt --fasta_file sample.scaffolds.fa --sample_name sample --id ticket_id --output_html sample_blast_report.html
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)


    parser.add_argument(
        "-b",
        "--blast_file",
        required=True,
        type=str,
        help="BLAST results file (required)",
    )
    parser.add_argument(
        "-f",
        "--fasta_file",
        required=True,
        type=str,
        help="FASTA file (required)",
    )
    parser.add_argument(
        "-s",
        "--sample_name",
        required=True,
        type=str,
        help="Sample name (required)",
    )
    parser.add_argument(
        "-i",
        "--id",
        required=False,
        type=str,
        help="Run ID to be in the report (optional)",
    )
    parser.add_argument(
        "-ml",
        "--min_qlen",
        default=200,
        type=int,
        help="Minimum query length to consider a BLAST hit (optional)",
    )
    parser.add_argument(
        "-mc",
        "--min_coverage",
        default=50,
        type=int,
        help="Minimum coverage to consider a BLAST hit (optional)",
    )
    parser.add_argument(
        "-o",
        "--output_html",
        required=True,
        type=str,
        help="Output HTML report file (required)",
    )
    parser.add_argument(
        "-sr"
        "--suggest_min_rows",
        type=int,
        default=20,
        help='Minimum number of rows (hits) required to consider auto-suggestion (default: 20)'
    )
    parser.add_argument(
        "-si"
        "--suggest_min_identity",
        type=float,
        default=90.0,
        help='Minimum max %% identity required to consider auto-suggestion (default: 90)'
    )
    parser.add_argument(
        "-sb"
        "--suggest_min_bitscore",
        type=float,
        default=300,
        help='Minimum max bitscore required to consider auto-suggestion (default: 300)'
    )
    parser.add_argument(
        "-ns"
        "--no_suggest",
        action='store_true',
        help='Disable automated genotype suggestion'
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help='DPI for output plots. Default: 300'
    )

    return parser.parse_args(args)

def extract_contig_headers(fasta_content):
    headers = []
    for line in fasta_content.split('\n'):
        if line.startswith('>'):
            headers.append(line)
    return headers

def encode_plot_to_base64(fig, dpi=300):
    buffer = BytesIO()
    fig.savefig(buffer, format='png', dpi=dpi, bbox_inches='tight', facecolor='white', edgecolor='none')
    buffer.seek(0)
    img_base64 = base64.b64encode(buffer.read()).decode('utf-8')
    buffer.close()
    return f"data:image/png;base64,{img_base64}"

def generate_report_data(blast_file, fasta_file, sample_name, id, min_qlen=200, min_coverage=50, suggest_enabled=True,
                         suggest_min_rows=20, suggest_min_identity=90.0, suggest_min_bitscore=300, plot_dpi=300):

    df = pd.read_csv(blast_file, sep="\t", header=0, index_col=0)

    unique_contigs = df['qaccver'].unique()
    df[['contig','temp1']] = df['qaccver'].str.split('_length_', expand=True)
    df[['length','coverage']] = df['temp1'].str.split('_cov_', expand=True)
    df = df.drop(['qaccver', 'temp1', 'length'], axis=1)
    df['coverage'] = pd.to_numeric(df['coverage'])
    df = df[df['qlen'] > min_qlen]
    df = df[df['coverage'] > min_coverage]

    if df.empty or df['contig'].nunique() == 0:
        raise ValueError("No contigs meet filtering criteria")

    # Suggestion logic
    suggestion = "Automatic suggestion disabled."
    if suggest_enabled:
        try:
            max_pident = df.loc[df['pident'].idxmax()]
            max_bitscore = df.loc[df['bitscore'].idxmax()]
            grouped_pident_medians = df.groupby('sscinames')['pident'].median()
            highest_median_pident = grouped_pident_medians.idxmax()
            grouped_bitscore_medians = df.groupby('sscinames')['bitscore'].median()
            highest_median_bitscore = grouped_bitscore_medians.idxmax()
            species_counts = df['sscinames'].value_counts()
            top_species_count = species_counts.get(max_pident['sscinames'], 0)
            if (
                (max_pident['sscinames'] == max_bitscore['sscinames']) and
                (top_species_count >= suggest_min_rows) and
                (df['pident'].max() >= suggest_min_identity) and
                (df['bitscore'].max() >= suggest_min_bitscore) and
                (highest_median_pident == max_pident['sscinames']) and
                (highest_median_bitscore == max_bitscore['sscinames'])
            ):
                suggestion = str(max_pident['sscinames'])
            else:
                suggestion = "Please do manual assessment."
        except Exception as e:
            print(f"⚠️ Suggestion logic failed: {e}")
            suggestion = "Please do manual assessment."

    # --------------------------------------------------------------
    # Count number of contigs / genotypes and prepare plotting data
    n_contigs = df['contig'].nunique()
    u_contigs = df[['contig','coverage','sscinames']].drop_duplicates()
    n_sscinames = df['sscinames'].nunique()
    l_contigs = df[['contig','qlen','sscinames']].drop_duplicates()

    if n_contigs == 0 or n_sscinames == 0:
        raise ValueError(f"Invalid contig or genotype count for {file_name}")

    pident_medians = df.groupby('sscinames')['pident'].median().sort_values(ascending=True)
    bitscore_medians = df.groupby('sscinames')['bitscore'].median().sort_values(ascending=True)
    coverage_medians = df.groupby('contig')['coverage'].median().sort_values(ascending=True)
    length_medians = df.groupby('contig')['qlen'].median().sort_values(ascending=True)

    pident_order = pident_medians.index.tolist()
    bitscore_order = bitscore_medians.index.tolist()
    coverage_order = coverage_medians.index.tolist()
    length_order = length_medians.index.tolist()

    if n_sscinames+n_contigs == 2:
        fig_height = 3
    elif n_sscinames+n_contigs > 8:
        fig_height = (n_sscinames+n_contigs)*0.8
    else:
        fig_height = n_sscinames+n_contigs

    h_ratio = 1 if n_contigs == 1 else n_contigs * 0.7
    legend_status = n_contigs > 1
    fig_width = 14 if legend_status else 10

    fig = plt.figure(figsize=(fig_width, fig_height))
    gs = fig.add_gridspec(2, 2, height_ratios=[n_sscinames, h_ratio])

    # Plot 1: identity score
    ax = fig.add_subplot(gs[0 , 0])
    ax.set_title("Figur 1: Identitet per genotyp")
    sns.pointplot(
        data=df, x="pident", y="sscinames", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", zorder=10, order=pident_order, legend=False
    )
    sns.stripplot(
        data=df, x="pident", y="sscinames", hue="contig",
        dodge=True, alpha=.4, legend=False, jitter=0.3, order=pident_order
    )
    if n_sscinames > 1:
        for i in range(n_sscinames):
            if i % 2 == 1:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_sscinames-0.3)
    ax.grid(True, axis="x", linestyle="--")
    ax.set(xlabel='BLAST identitet (%)', ylabel='Genotyp')

    # Plot 2: bit score
    ax = fig.add_subplot(gs[0 , 1])
    ax.set_title("Figur 2: Bit score per genotyp")
    sns.pointplot(
        data=df, x="bitscore", y="sscinames", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", zorder=10, order=bitscore_order, legend=legend_status
    )
    sns.stripplot(
        data=df, x="bitscore", y="sscinames", hue="contig",
        dodge=True, alpha=.4, legend=False, jitter=0.3, order=bitscore_order
    )
    if n_sscinames > 1:
        for i in range(n_sscinames):
            if i % 2 == 1:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_sscinames-0.3)
    ax.grid(True, axis="x", linestyle="--")
    if legend_status:
        ax.legend().set_title("Contig")
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
    ax.set(xlabel='BLAST bit score', ylabel=' ')

    # Plot 3: coverage
    ax = fig.add_subplot(gs[1, 0])
    ax.set_title("Figur 3: Täckning per contig")
    sns.pointplot(
        data=u_contigs, x="coverage", y="contig", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", order=coverage_order
    )
    if n_contigs > 1:
        for i in range(n_contigs):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_contigs-0.3)
    ax.grid(True, axis="x", linestyle="--")
    ax.set(xlabel='Täckning (x)', ylabel='Contig')
    ax.set_xscale('log')
    ax.set_xlim(left=50, right=10*df['coverage'].max())

    # Plot 4: length
    ax = fig.add_subplot(gs[1, 1])
    ax.set_title("Figur 4: Längd per contig")
    sns.pointplot(
        data=l_contigs, x="qlen", y="contig", hue="contig",
        dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
        markers="d", markersize=6, linestyle="none", order=length_order, legend=False
    )
    if n_contigs > 1:
        for i in range(n_contigs):
            if i % 2 == 0:
                ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
    ax.set_ylim(-0.5, n_contigs-0.3)
    ax.grid(True, axis="x", linestyle="--")
    ax.set(xlabel='Längd (bp)', ylabel=' ')
    ax.set_xlim(left=200, right=1.1*df['qlen'].max())

    fig.subplots_adjust(hspace=20, wspace=20)
    fig.tight_layout()

    img_buffer = BytesIO()
    fig.savefig(img_buffer, format='png', dpi=plot_dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    img_buffer.seek(0)
    img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
    img_data_uri = f"data:image/png;base64,{img_base64}"

    # Filter FASTA to only selected contigs
    filtered_fasta = filter_fasta_by_contigs(fasta_file, unique_contigs)

    contigs_content_html = filtered_fasta.replace('\n', '<br>')
    contig_headers = extract_contig_headers(filtered_fasta)
    contigs_summary = '<br>'.join(contig_headers)
    contig_count = len(contig_headers)

    # Warning
    warningtext = ''
    if df['pident'].max() < 90:
        warningtext = "&#9888; OBS! Highest identity < 90%"

    return {
        'sample_name': sample_name,
        'id': id,
        'time_stamp': datetime.datetime.now().strftime("%Y-%m-%d"),
        'is_error_report': False,
        'warningtext': warningtext,
        'img_data_uri': img_data_uri,
        'contigs_content': contigs_content_html,
        'contigs_summary': contigs_summary,
        'contig_count': contig_count,
        'suggestion': suggestion
    }


def generate_error_report_data(blast_file, fasta_file, sample_name, id, plot_dpi=300):
    file_name = os.path.basename(blast_file)

    img_data_uri = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="
    contigs_content = ""
    contigs_summary = ""
    contig_count = 0

    try:
        df = pd.read_csv(blast_file, sep="\t", header=0, index_col=0)
        unique_contigs = df['qaccver'].unique()

        df[['contig','temp1']] = df['qaccver'].str.split('_length_',expand=True)
        df[['length','coverage']] = df['temp1'].str.split('_cov_',expand=True)
        df = df.drop(['qaccver', 'temp1', 'length'], axis=1)
        df["coverage"] = pd.to_numeric(df["coverage"])

        if df.empty or df['contig'].nunique() == 0:
            print(f"Warning: No valid data found in BLAST file {file_name}")
            return img_data_uri, contigs_content, contigs_summary, contig_count

        n_contigs = df['contig'].nunique()
        u_contigs = df[['contig','coverage','sscinames']].drop_duplicates()
        n_sscinames = df['sscinames'].nunique()
        l_contigs = df[['contig','qlen','sscinames']].drop_duplicates()

        if n_contigs == 0 or n_sscinames == 0:
            print(f"Warning: Invalid contig or genotype count for {file_name}")
            return img_data_uri, contigs_content, contigs_summary, contig_count

        pident_medians = df.groupby('sscinames')['pident'].median().sort_values(ascending=True)
        bitscore_medians = df.groupby('sscinames')['bitscore'].median().sort_values(ascending=True)
        coverage_medians = df.groupby('contig')['coverage'].median().sort_values(ascending=True)
        length_medians = df.groupby('contig')['qlen'].median().sort_values(ascending=True)

        pident_order = pident_medians.index.tolist()
        bitscore_order = bitscore_medians.index.tolist()
        coverage_order = coverage_medians.index.tolist()
        length_order = length_medians.index.tolist()

        if n_sscinames+n_contigs == 2:
            fig_height = 3
        elif n_sscinames+n_contigs > 8:
            fig_height = (n_sscinames+n_contigs)*0.8
        else:
            fig_height = n_sscinames+n_contigs

        h_ratio = 1 if n_contigs == 1 else n_contigs * 0.7
        legend_status = n_contigs > 1
        fig_width = 14 if legend_status else 10

        fig = plt.figure(figsize=(fig_width, fig_height))
        gs = fig.add_gridspec(2, 2, height_ratios=[n_sscinames, h_ratio])

        # (plots identical to normal path)
        ax = fig.add_subplot(gs[0 , 0])
        ax.set_title("Figur 1: Identitet per genotyp")
        sns.pointplot(
            data=df, x="pident", y="sscinames", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", zorder=10, order=pident_order, legend=False
        )
        sns.stripplot(
            data=df, x="pident", y="sscinames", hue="contig",
            dodge=True, alpha=.4, legend=False, jitter=0.3, order=pident_order
        )
        if n_sscinames > 1:
            for i in range(n_sscinames):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_sscinames-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='BLAST identitet (%)', ylabel='Genotyp')

        ax = fig.add_subplot(gs[0 , 1])
        ax.set_title("Figur 2: Bit score per genotyp")
        sns.pointplot(
            data=df, x="bitscore", y="sscinames", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", zorder=10, order=bitscore_order, legend=legend_status
        )
        sns.stripplot(
            data=df, x="bitscore", y="sscinames", hue="contig",
            dodge=True, alpha=.4, legend=False, jitter=0.3, order=bitscore_order
        )
        if n_sscinames > 1:
            for i in range(n_sscinames):
                if i % 2 == 1:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_sscinames-0.3)
        ax.grid(True, axis="x", linestyle="--")
        if legend_status:
            ax.legend().set_title("Contig")
            sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))
        ax.set(xlabel='BLAST bit score', ylabel=' ')

        ax = fig.add_subplot(gs[1, 0])
        ax.set_title("Figur 3: Täckning per contig")
        sns.pointplot(
            data=u_contigs, x="coverage", y="contig", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", order=coverage_order
        )
        if n_contigs > 1:
            for i in range(n_contigs):
                if i % 2 == 0:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_contigs-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='Täckning (x)', ylabel='Contig')
        ax.set_xscale('log')
        max_coverage = df['coverage'].max()
        if max_coverage > 0:
            ax.set_xlim(left=max(0.1, df['coverage'].min()/2), right=10*max_coverage)
        else:
            ax.set_xlim(left=0.1, right=100)

        ax = fig.add_subplot(gs[1, 1])
        ax.set_title("Figur 4: Längd per contig")
        sns.pointplot(
            data=l_contigs, x="qlen", y="contig", hue="contig",
            dodge=.8 - .8 / n_contigs, palette="dark", errorbar=None,
            markers="d", markersize=6, linestyle="none", order=length_order, legend=False
        )
        if n_contigs > 1:
            for i in range(n_contigs):
                if i % 2 == 0:
                    ax.axhspan(i - 0.5, i + 0.5, facecolor='gray', alpha=0.2, zorder=-1)
        ax.set_ylim(-0.5, n_contigs-0.3)
        ax.grid(True, axis="x", linestyle="--")
        ax.set(xlabel='Längd (bp)', ylabel=' ')
        max_length = df['qlen'].max()
        if max_length > 0:
            ax.set_xlim(left=max(1, df['qlen'].min()/2), right=1.1*max_length)
        else:
            ax.set_xlim(left=1, right=1000)

        fig.subplots_adjust(hspace = 20, wspace=20)
        fig.tight_layout()

        img_buffer = BytesIO()
        fig.savefig(img_buffer, format='png', dpi=plot_dpi, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        img_buffer.seek(0)
        img_base64 = base64.b64encode(img_buffer.read()).decode('utf-8')
        img_data_uri = f"data:image/png;base64,{img_base64}"

    except Exception as e:
        print(f"Error generating plots for {file_name}: {e}")
        img_data_uri = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNkYPhfDwAChwGA60e6kgAAAABJRU5ErkJggg=="

    # Read fasta file (unfiltered original)
    try:
        filtered_fasta = filter_fasta_by_contigs(fasta_file, unique_contigs)
        contigs_content_html = filtered_fasta.replace('\n', '<br>')
        if contigs_content_html:
                contig_headers = extract_contig_headers(filtered_fasta)
                contig_count = len(contig_headers)
                contigs_summary = '<br>'.join(contig_headers)
        else:
            print(f"Warning: FASTA file not found: {seq_file_original}")
    except Exception as e:
        print(f"Error reading FASTA file for {file_name}: {e}")

    return img_data_uri, contigs_content, contigs_summary, contig_count

    return {
        'sample_name': sample_name,
        'id': id,
        'time_stamp': datetime.datetime.now().strftime("%Y-%m-%d"),
        'is_error_report': True,
        'img_data_uri': img_data_uri,
        'contigs_content': contigs_content_html,
        'contigs_summary': contigs_summary,
        'contig_count': contig_count
    }


def render_report(output_file, template, data, css_content):

    html_content = template.render(
        css_content=css_content,
        data = data
    )

    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

def filter_fasta_by_contigs(fasta_path, contigs_to_keep):
    """Read a FASTA file and return only sequences whose headers match the contigs_to_keep list."""
    filtered_seqs = []
    keep = False
    header = None
    seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # If we were collecting the previous sequence, save it
                if header and keep:
                    filtered_seqs.append((header, "".join(seq)))

                header = line
                seq = []

                # Check if contig name appears in header
                # Example header: >NODE_1_length_7412_cov_139.22
                contig_name = header[1:]

                keep = contig_name in contigs_to_keep
            else:
                seq.append(line)

        # Save last seq
        if header and keep:
            filtered_seqs.append((header, "".join(seq)))

    # Convert into FASTA format
    filtered_fasta = "\n".join(f"{h}\n{s}" for h, s in filtered_seqs)
    return filtered_fasta

def main(args=None, is_error=False):
    args = parser_args(args)

    # Build paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    asset_path = os.path.join(script_dir, "../assets")
    css_path = os.path.join(asset_path, "blast_report_template.css")

    # --- Load CSS content
    with open(css_path, "r", encoding="utf-8") as css_file:
        css_content = css_file.read()

    # --- Load Jinja2 environment from the template's folder
    template_dir = os.path.abspath(asset_path)

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(['html', 'xml'])
    )

    template = env.get_template("blast_report_template.html")

    if is_error:
        data = generate_error_report_data(args.blast_file, args.fasta_file, args.sample_name, args.id, plot_dpi=args.dpi)
    else:
        data = generate_report_data(args.blast_file, args.fasta_file, args.sample_name, args.id, args.min_qlen, args.min_coverage, plot_dpi=args.dpi)

    render_report(args.output_html, template, data, css_content)

    print(f"📁 Report saved to: {args.output_html}")

if __name__ == "__main__":
    main()
