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
from jinja2 import Environment, FileSystemLoader

# Path to template folder
TEMPLATE_DIR = os.path.dirname(os.path.abspath(__file__))

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
        help='Minimum max %% identity required to consider auto-suggestion (default: 90)')
    parser.add_argument(
        "-sb"
        "--suggest_min_bitscore",
        type=float,
        default=300,
        help='Minimum max bitscore required to consider auto-suggestion (default: 300)')
    parser.add_argument(
        "-ns"
        "--no_suggest",
        action='store_true',
        help='Disable automated genotype suggestion')
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

def generate_report_data(blast_file, output_base, suggest_enabled=True,
                         suggest_min_rows=20, suggest_min_identity=90.0, suggest_min_bitscore=300, dpi=300):

    df = pd.read_csv(blast_file, header=0)
    df[['contig','temp1']] = df['qseqid'].str.split('_length_', expand=True)
    df[['length','coverage']] = df['temp1'].str.split('_cov_', expand=True)
    df = df.drop(['qseqid', 'temp1', 'length'], axis=1)
    df['coverage'] = pd.to_numeric(df['coverage'])
    df = df[df['qlen'] > 200]
    df = df[df['coverage'] > 50]

    if df.empty or df['contig'].nunique() == 0:
        raise ValueError("No contigs meet filtering criteria")

    # Suggestion logic
    suggestion = "Automatic suggestion disabled."
    if suggest_enabled:
        try:
            max_pident = df.loc[df['pident'].idxmax()]
            max_bitscore = df.loc[df['bitscore'].idxmax()]
            species_counts = df['sscinames'].value_counts()
            top_species_count = species_counts.get(max_pident['sscinames'], 0)
            if (max_pident['sscinames'] == max_bitscore['sscinames'] and
                top_species_count >= suggest_min_rows and
                df['pident'].max() >= suggest_min_identity and
                df['bitscore'].max() >= suggest_min_bitscore):
                suggestion = str(max_pident['sscinames'])
            else:
                suggestion = "Please do manual assessment."
        except:
            suggestion = "Please do manual assessment."

    # Plotting example (identity per genotype)

    fig, ax = plt.subplots(figsize=(10,6))
    sns.pointplot(data=df, x="pident", y="scomname", hue="contig", dodge=True, errorbar=None, ax=ax)
    ax.set_title("Identity per genotype")
    img_data_uri = encode_plot_to_base64(fig, dpi=dpi)
    plt.close(fig)

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
        'seq_name': os.path.basename(blast_file).replace('.blast', ''),
        'time_stamp': datetime.datetime.now().strftime("%Y-%m-%d"),
        'is_error_report': False,
        'warningtext': warningtext,
        'img_data_uri': img_data_uri,
        'contigs_content': contigs_content_html,
        'contigs_summary': contigs_summary,
        'contig_count': contig_count,
        'suggestion': suggestion
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

def main(args=None):
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

    data = generate_report_data(args.blast_file, args.fasta_file, args.sample_name, args.id)

    render_report(args.output_html, template, data, css_content)

    print(f"Report saved to: {args.output_html}")

if __name__ == "__main__":
    main()
