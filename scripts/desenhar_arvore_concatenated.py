#!/usr/bin/env python3
"""
Script to generate phylograms and cladograms in PDF with Bootstrap support (RAxML)
and Posterior Probability (MrBayes).
Usage: python3 script.py <raxml_file> <mrbayes_file> <output_name>
"""
import os
import sys
import re
import math
import toytree
import toyplot
import toyplot.pdf
from pathlib import Path

# ============================================================
# CONFIGURATION
# ============================================================
OUTGROUP_PREFIX = "Gamszarea"
BOLD_GROUP = "Simplicillium sp"

COL_PALETTE = {
    'very_high': '#1a4d0a',    # ≥95%
    'high': '#2a581d',         # 81-94%
    'medium': '#7bb142',       # 51-80%
    'low': '#e4f2cb'           # <51%
}

def read_tree(file):
    """Reads the RAxML tree and roots it at the outgroup."""
    try:
        print(f"     Reading RAxML tree: {Path(file).name}...")
        tre = toytree.tree(file)
        tre = tre.ladderize(True)

        full_name = [tip for tip in tre.get_tip_labels() if OUTGROUP_PREFIX in tip]
        if full_name:
            tre = tre.root(full_name[0])
            print(f"     Tree rooted at: {full_name[0]}")
        return tre
    except Exception as e:
        print(f"     Error reading tree: {e}")
        sys.exit(1)

def remove_underscores_labels(tree):
    """Removes underscores from taxon names."""
    tips = tree.get_tip_labels()
    new_map = {}
    for tip in tips:
        new_name = tip.replace("_", " ")
        new_map[tip] = new_name

    for node in tree.traverse():
        if node.name in new_map:
            node.name = new_map[node.name]

    return tree

def extract_mrbayes_support_nexus(mrbayes_file):
    """Extracts PP values from MrBayes NEXUS tree."""
    print(f"    Extracting PP from MrBayes: {Path(mrbayes_file).name}...")
    try:
        with open(mrbayes_file, 'r') as f:
            content = f.read()
        tree_match = re.search(r'tree\s+.*?\s*=\s*(.*?);', content, re.DOTALL | re.IGNORECASE)
        if not tree_match:
            print("      No tree found in MrBayes file")
            return []
        tree_string = tree_match.group(1)
        # FIXED: removed literal newlines inside the raw string
        probs = re.findall(r'\[&prob=([\d.e\-\+]+)', tree_string)
        return [f"{float(p):.2f}" if float(p) <= 1 else f"{float(p)/100:.2f}" for p in probs]
    except Exception as e:
        print(f"     Error extracting PP: {e}")
        return []

def combine_support(tree_raxml, probs_mrbayes):
    """Combines Bootstrap and Posterior Probability."""
    supports = tree_raxml.get_node_data("support")
    labels = []

    for i in range(tree_raxml.nnodes):
        bs = supports.get(i, None)
        label = ""

        if bs is not None and str(bs).lower() != 'nan' and str(bs) != "":
            try:
                bs_val = str(int(float(bs)))
                if not tree_raxml[i].is_leaf() and i < len(probs_mrbayes):
                    label = f"{bs_val}/{probs_mrbayes[i]}"
                else:
                    label = bs_val
            except:
                label = ""

        labels.append(label)

    return labels

def draw_tree_pdf(tre, labels, output_file, cladogram=False):
    """Generates a phylogram or cladogram PDF with aligned support values."""
    tree_type = "cladogram" if cladogram else "phylogram"
    print(f"     Generating {tree_type} PDF...")

    num_taxa = len(tre.get_tip_labels())
    height = max(1200, num_taxa * 80)
    width = 2000

    canvas = toyplot.Canvas(width=width, height=height, style={"background-color": "white"})
    axes = canvas.cartesian()
    axes.show = False

    axes.rectangle(
        -15, 4,
        14.75, 25.2,
        opacity=0.25,
        color='#EC94F4',
    )

    axes.rectangle(
        -15, 4,
        -0.5, 14.75,
        opacity=0.25,
        color='#94E3F4',
    )

    axes.rectangle(
        -15, 4,
        25.2, 74.3,
        opacity=0.25,
        color='#94E3F4',
    )

    supports = tre.get_node_data("support")
    node_colors = []
    for i in range(tre.nnodes):
        val = supports.get(i, 0)
        try:
            v = float(val)
            if math.isnan(v): v = 0
            if v >= 95: node_colors.append(COL_PALETTE['very_high'])
            elif v >= 81: node_colors.append(COL_PALETTE['high'])
            elif v >= 51: node_colors.append(COL_PALETTE['medium'])
            else: node_colors.append(COL_PALETTE['low'])
        except:
            node_colors.append(COL_PALETTE['low'])

    labels_clean = [label if label else "" for label in labels]

    for node in tre.traverse():
        if node.name and "**" in node.name:
            node.name = node.name.replace("**", "")

    tre.draw(
        axes=axes,
        node_labels=labels_clean,
        node_markers="o",
        node_sizes=28,
        node_colors=node_colors,
        tip_labels_style={
            "font-size": "10px",
            "font-weight": "normal"
        },
        node_labels_style={
            "font-size": "8px",
            "fill": "red",
            "font-weight": "bold"
        },
        use_edge_lengths=not cladogram,
        edge_widths=1.5,
        tip_labels=True
    )

    try:
        toyplot.pdf.render(canvas, output_file)
        print(f"    ✓ {tree_type.capitalize()} PDF generated: {output_file}")
    except Exception as e:
        print(f"    Error rendering PDF: {e}")
        svg_output = output_file.replace('.pdf', '.svg')
        toyplot.svg.render(canvas, svg_output)
        print(f"    ✓ SVG generated successfully: {svg_output}")

def main():
    """Main function: MrBayes -> RAxML -> Drawing."""
    if len(sys.argv) < 3:
        print("\n Error: Not enough arguments!")
        print("Usage: python3 script.py <RAxML_file> <MrBayes_file> <output_name>")
        sys.exit(1)

    raxml_f = sys.argv[1]
    mrbayes_f = sys.argv[2]

    if len(sys.argv) >= 4:
        output_f = sys.argv[3].strip()
    else:
        output_f = "final_tree"
    if output_f.endswith(".pdf"):
        output_f = output_f[:-4]

    print("\n" + "="*60)
    print(" PHYLOGENETIC PROCESSOR")
    print("="*60)

    pp_values = extract_mrbayes_support_nexus(mrbayes_f)
    print("    ✓ Posterior Probabilities extracted.")

    tree = read_tree(raxml_f)
    tree = remove_underscores_labels(tree)
    print("    ✓ Species names cleaned.")

    labels = combine_support(tree, pp_values)

    phylo_file = f"{output_f}.pdf"
    draw_tree_pdf(tree, labels, phylo_file, cladogram=False)

    clado_file = f"{output_f}_cladogram.pdf"
    draw_tree_pdf(tree, labels, clado_file, cladogram=True)

    print(f"\n Completed! Files saved:")
    print(f"   → Phylogram:  {phylo_file}")
    print(f"   → Cladogram: {clado_file}")
    print("="*60 + "\n")

if __name__ == "__main__":
    main()
