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
    python blast_report.py --output_html <sample>_blast_report.html
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
            species_counts = df['scomname'].value_counts()
            top_species_count = species_counts.get(max_pident['scomname'], 0)
            if (max_pident['scomname'] == max_bitscore['scomname'] and
                top_species_count >= suggest_min_rows and
                df['pident'].max() >= suggest_min_identity and
                df['bitscore'].max() >= suggest_min_bitscore):
                suggestion = str(max_pident['scomname'])
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

    # Contigs
    seq_file = f"{output_base}/ev_contig/{os.path.basename(blast_file).replace('.blast', '_200bp_minCov50.fasta')}"
    with open(seq_file, 'r') as f:
        contigs_content = f.read()
    contig_headers = extract_contig_headers(contigs_content)
    contig_count = len(contig_headers)
    contigs_summary = '<br>'.join(contig_headers)
    contigs_content_html = contigs_content.replace('\n', '<br>')

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

def render_report(template_file, output_file, data):
    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template(template_file)
    html_content = template.render(**data)
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)

def main(args=None):
    args = parser_args(args)
    data = generate_report_data(args.blast_file, "./")
    render_report(args.output_html, os.path.join("./", f"{data['seq_name']}.html"), data)
    print(f"Report saved to: {args.output_html}")

if __name__ == "__main__":
    main()
