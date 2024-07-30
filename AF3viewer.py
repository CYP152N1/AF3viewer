import json
import matplotlib.pyplot as plt
import numpy as np
import os
import argparse
from Bio import PDB
import csv
from PIL import Image

def plot_plddt(plddt_data, folder_name, file_index):
    plt.figure()
    plt.plot(plddt_data)
    plt.xlabel('Residue Index')
    plt.ylabel('pLDDT')
    plt.title(f'pLDDT Scores for {folder_name}_model_{file_index}')
    plt.grid(True)
    output_path = os.path.join(folder_name, f'{folder_name}_plddt_{file_index}.png')
    plt.savefig(output_path)
    plt.close()

def plot_sequence_frequencies(sequences, folder_name):
    amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
    frequencies = {aa: 0 for aa in amino_acids}

    for seq in sequences:
        for aa in seq:
            if aa in frequencies:
                frequencies[aa] += 1

    plt.figure()
    plt.bar(frequencies.keys(), frequencies.values())
    plt.xlabel('Amino Acid')
    plt.ylabel('Frequency')
    plt.title(f'Amino Acid Frequencies in {folder_name}')
    plt.grid(True)
    output_path = os.path.join(folder_name, f'{folder_name}_amino_acid_frequencies.png')
    plt.savefig(output_path)
    plt.close()

def process_job_request(folder_name):
    job_request_path = os.path.join(folder_name, f'{folder_name}_job_request.json')
    with open(job_request_path, 'r') as f:
        job_data = json.load(f)

    sequences = []
    for chain in job_data[0]['sequences']:
        if 'proteinChain' in chain:
            sequences.append(chain['proteinChain']['sequence'])

    plot_sequence_frequencies(sequences, folder_name)

def plot_contact_probs(contact_probs, folder_name, file_index):
    plt.figure()
    plt.imshow(contact_probs, cmap='hot', interpolation='nearest')
    plt.colorbar()
    plt.title(f'Contact Probability Map for {folder_name}_model_{file_index}')
    output_path = os.path.join(folder_name, f'{folder_name}_contact_probs_{file_index}.png')
    plt.savefig(output_path)
    plt.close()

    # Contact probabilities data to CSV
    csv_output_path = os.path.join(folder_name, f'{folder_name}_contact_probs_{file_index}.csv')
    with open(csv_output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(contact_probs)

def plot_pae(pae_data, folder_name, file_index):
    plt.figure()
    plt.imshow(pae_data, cmap='Greens_r', interpolation='nearest', vmin=0, vmax=32)
    plt.colorbar()
    plt.clim(0, 32)
    plt.title(f'PAE of {folder_name}_model_{file_index}')
    output_path = os.path.join(folder_name, f'{folder_name}_pae_{file_index}.png')
    plt.savefig(output_path)
    plt.close()

    # PAE data to CSV
    csv_output_path = os.path.join(folder_name, f'{folder_name}_pae_{file_index}.csv')
    with open(csv_output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(pae_data)

def plot_matrix(matrix_data, title, output_path, cmap='Greens', labels=None):
    # Ensure matrix_data is a numpy array of floats
    try:
        matrix_data = np.array(matrix_data, dtype=float)
    except ValueError:
        print(f"Error: Non-numeric data found in matrix_data for {title}")
        return

    plt.figure()
    plt.imshow(matrix_data, cmap=cmap, interpolation='nearest', vmin=0, vmax=1)
    plt.colorbar()
    plt.clim(0, 1)
    plt.title(title)
    if labels:
        plt.xticks(np.arange(len(labels)), labels)
        plt.yticks(np.arange(len(labels)), labels)
    plt.savefig(output_path)
    plt.close()

def save_matrix_to_csv(matrix_data, output_path):
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(matrix_data)

def plot_and_save_chain_pair_matrices(data, folder_name, file_index, labels):
    chain_pair_iptm = data.get('chain_pair_iptm', [])
    chain_pair_pae_min = data.get('chain_pair_pae_min', [])

    if chain_pair_iptm:
        # Plot and save chain_pair_iptm
        iptm_title = f'Chain Pair iPTM for {folder_name}_model_{file_index}'
        iptm_output_path = os.path.join(folder_name, f'{folder_name}_chain_pair_iptm_{file_index}.png')
        plot_matrix(chain_pair_iptm, iptm_title, iptm_output_path, cmap='Greens', labels=labels)

        iptm_csv_path = os.path.join(folder_name, f'{folder_name}_chain_pair_iptm_{file_index}.csv')
        save_matrix_to_csv(chain_pair_iptm, iptm_csv_path)

    if chain_pair_pae_min:
        # Save chain_pair_pae_min to CSV (but do not plot PNG)
        pae_csv_path = os.path.join(folder_name, f'{folder_name}_chain_pair_pae_min_{file_index}.csv')
        save_matrix_to_csv(chain_pair_pae_min, pae_csv_path)

def plot_ca_structure_from_cif(file_path, folder_name, file_index):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure(f'{folder_name}_model_{file_index}', file_path)
    chains = list(structure.get_chains())
    
    color_map = plt.get_cmap('coolwarm')
    marker_styles = ['o', '^', 's', 'p', '*', 'X', 'D', '<', '>', 'v', 'h', 'H', '+', 'x', '|', '_']

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for chain_index, chain in enumerate(chains):
        atoms = [atom for atom in chain.get_atoms() if atom.get_name() in ['CA', 'P']]
        if not atoms:
            continue
        coords = np.array([atom.get_coord() for atom in atoms])
        res_nums = np.array([atom.get_parent().id[1] for atom in atoms])
        if res_nums.max() - res_nums.min() == 0:
            norm_res_nums = np.zeros_like(res_nums)
        else:
            norm_res_nums = (res_nums - res_nums.min()) / (res_nums.max() - res_nums.min())
        colors = color_map(norm_res_nums)

        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=colors, s=20, label=f'Chain {chain.id}', marker=marker_styles[chain_index % len(marker_styles)])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Cα Atoms of {folder_name}_model_{file_index}')
    ax.legend()
    output_path = os.path.join(folder_name, f'{folder_name}_structure_rank_{file_index}.png')
    plt.savefig(output_path)
    plt.close()


def plot_chain_colored_structure(file_path, folder_name, file_index):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure(f'{folder_name}_model_{file_index}', file_path)
    chains = list(structure.get_chains())

    color_map = plt.get_cmap('rainbow')

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    for chain_index, chain in enumerate(chains):
        atoms = [atom for atom in chain.get_atoms() if atom.get_name() in ['CA', 'P']]
        if not atoms:
            continue
        coords = np.array([atom.get_coord() for atom in atoms])
        color = color_map(chain_index / len(chains))

        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=[color], s=20, label=f'Chain {chain.id}')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'Chain-Colored Structure of {folder_name}_model_{file_index}_rank')
    ax.legend()
    output_path = os.path.join(folder_name, f'{folder_name}_chain_rank_{file_index}.png')
    plt.savefig(output_path)
    plt.close()

    # Plot with 180 degrees rotation around the Z-axis
    fig_rotated = plt.figure()
    ax_rotated = fig_rotated.add_subplot(111, projection='3d')

    for chain_index, chain in enumerate(chains):
        atoms = [atom for atom in chain.get_atoms() if atom.get_name() in ['CA', 'P']]
        if not atoms:
            continue
        coords = np.array([atom.get_coord() for atom in atoms])
        color = color_map(chain_index / len(chains))

        ax_rotated.scatter(-coords[:, 0], -coords[:, 1], coords[:, 2], c=[color], s=20, label=f'Chain {chain.id}')

    ax_rotated.set_xlabel('-X')
    ax_rotated.set_ylabel('-Y')
    ax_rotated.set_zlabel('Z')
    ax_rotated.set_title(f'Chain-Colored Structure of {folder_name}_model_{file_index}_rank')
    ax_rotated.legend()
    output_path_rotated = os.path.join(folder_name, f'{folder_name}_rotated_chain_rank_{file_index}.png')
    plt.savefig(output_path_rotated)
    plt.close()

def plot_plddt_from_cif(file_path, folder_name, file_index):
    parser = PDB.MMCIFParser()
    structure = parser.get_structure(f'{folder_name}_model_{file_index}', file_path)
    residues = []

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.id in ['CA', 'P']:  # CAまたはP原子を使用
                        residues.append({
                            'label_seq_id': residue.id[1],
                            'B_iso_or_equiv': atom.bfactor,
                            'chain_id': chain.id
                        })

    if residues:
        seq_ids = [res['label_seq_id'] for res in residues]
        b_factors = [res['B_iso_or_equiv'] for res in residues]
        chain_ids = [res['chain_id'] for res in residues]
        
        unique_chains = sorted(set(chain_ids))
        chain_colors = {chain: plt.get_cmap('rainbow')(i / len(unique_chains)) for i, chain in enumerate(unique_chains)}
        chain_images = []

        # Plot all chains with different colors
        plt.figure()
        for chain in unique_chains:
            chain_indices = [i for i, x in enumerate(chain_ids) if x == chain]
            plt.scatter([seq_ids[i] for i in chain_indices], [b_factors[i] for i in chain_indices],
                        color=chain_colors[chain], label=f'Chain {chain}', s=5)
        
        plt.xlabel('Residue Index')
        plt.ylabel('pLDDT')
        plt.ylim(0, 100)
        plt.title(f'Cα-pLDDT Scores of {folder_name}_model_{file_index}')
        plt.legend()
        output_path = os.path.join(folder_name, f'{folder_name}_cif_plddt_{file_index}.png')
        plt.savefig(output_path)
        plt.close()

        # Save data to CSV
        csv_output_path = os.path.join(folder_name, f'{folder_name}_cif_plddt_{file_index}.csv')
        with open(csv_output_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(['Residue Index', 'pLDDT', 'Chain'])
            writer.writerows(zip(seq_ids, b_factors, chain_ids))

        # Plot each chain separately and save images
        for chain in unique_chains:
            chain_indices = [i for i, x in enumerate(chain_ids) if x == chain]
            plt.figure()
            plt.scatter([seq_ids[i] for i in chain_indices], [b_factors[i] for i in chain_indices],
                        color=chain_colors[chain], s=5)
            plt.xlabel('Residue Index')
            plt.ylabel('pLDDT')
            plt.ylim(0, 100)
            plt.title(f'Chain {chain} pLDDT of {folder_name}_model_{file_index}')
            output_path = os.path.join(folder_name, f'{folder_name}_chain_{chain}_plddt_{file_index}.png')
            plt.savefig(output_path)
            plt.close()
            chain_images.append(output_path)
        
        # Create GIF for chain pLDDT plots
        create_gif(chain_images, os.path.join(folder_name, f'chain_{file_index}_plddt.gif'))

def load_confidence_data(folder_name):
    confidence_data = []
    for i in range(5):
        file_path = os.path.join(folder_name, f'{folder_name}_summary_confidences_{i}.json')
        with open(file_path, 'r') as f:
            data = json.load(f)
        data['model_index'] = i
        confidence_data.append(data)
    return confidence_data

def rank_and_save_cif(confidence_data, folder_name, rank_by='ranking_score'):
    if rank_by not in ['fraction_disordered', 'iptm', 'ptm', 'ranking_score']:
        raise ValueError(f"Invalid ranking criteria: {rank_by}")

    sorted_data = sorted(confidence_data, key=lambda x: x[rank_by], reverse=True)
    for rank, data in enumerate(sorted_data, start=1):
        original_cif_path = os.path.join(folder_name, f'{folder_name}_model_{data["model_index"]}.cif')
        new_cif_path = os.path.join(folder_name, f'{folder_name}_model_{data["model_index"]}_rank_{rank}.cif')
        os.system(f'cp {original_cif_path} {new_cif_path}')
    return sorted_data

def align_models(sorted_data, folder_name, align_chains=None):
    parser = PDB.MMCIFParser()
    superimposer = PDB.Superimposer()
    ref_structure = parser.get_structure('ref', os.path.join(folder_name, f'{folder_name}_model_{sorted_data[0]["model_index"]}_rank_1.cif'))

    for rank, data in enumerate(sorted_data, start=1):
        model_path = os.path.join(folder_name, f'{folder_name}_model_{data["model_index"]}_rank_{rank}.cif')
        model_structure = parser.get_structure(f'model_{rank}', model_path)
        ref_atoms = []
        model_atoms = []

        for ref_chain, model_chain in zip(ref_structure[0], model_structure[0]):
            if align_chains and ref_chain.id not in align_chains:
                continue
            ref_atoms.extend([atom for atom in ref_chain.get_atoms() if atom.get_name() in ['CA', 'P']])
            model_atoms.extend([atom for atom in model_chain.get_atoms() if atom.get_name() in ['CA', 'P']])

        superimposer.set_atoms(ref_atoms, model_atoms)
        superimposer.apply(model_structure.get_atoms())

        align_cif_path = os.path.join(folder_name, f'{folder_name}_model_{data["model_index"]}_rank_{rank}_align.cif')
        io = PDB.MMCIFIO()
        io.set_structure(model_structure)
        io.save(align_cif_path)

def create_gif(image_files, output_path):
    images = [Image.open(image_file) for image_file in image_files]
    images[0].save(output_path, save_all=True, append_images=images[1:], loop=0, duration=500)

def generate_html_report(folder_name):
    html_content = "<html><head><title>Analysis Report</title></head><body>"
    html_content += "<h1>Analysis Report</h1>"

    # Load and display job request details
    job_request_path = os.path.join(folder_name, f'{folder_name}_job_request.json')
    with open(job_request_path, 'r') as f:
        job_data = json.load(f)
    
    html_content += "<h2>Job Request Details</h2>"
    for item in job_data[0]['sequences']:
        if 'proteinChain' in item:
            protein_chain = item['proteinChain']
            sequence = protein_chain['sequence']
            count = protein_chain['count']
            length = len(sequence)
            html_content += f"<h3>Protein Chain (Count: {count}, Length: {length})</h3>"
            html_content += f"<textarea readonly style='width:100%; height:40px;'>{sequence}</textarea><br>"
            fasta_content = f">Protein Chain\n{sequence}"
            html_content += f"<a href='data:text/plain;charset=utf-8,{fasta_content}' download='protein_chain.fasta'><button>Save as Fasta</button></a><br><br>"

        if 'ligand' in item:
            ligand = item['ligand']
            ligand_name = ligand['ligand']
            count = ligand['count']
            html_content += f"<h3>Ligand (Count: {count})</h3>"
            html_content += f"<p>Name: {ligand_name}</p><br>"

        if 'ion' in item:
            ion = item['ion']
            ion_name = ion['ion']
            count = ion['count']
            html_content += f"<h3>Ion (Count: {count})</h3>"
            html_content += f"<p>Name: {ion_name}</p><br>"

        if 'dnaSequence' in item:
            dna_sequence = item['dnaSequence']
            sequence = dna_sequence['sequence']
            count = dna_sequence['count']
            length = len(sequence)
            html_content += f"<h3>DNA Sequence (Count: {count}, Length: {length})</h3>"
            html_content += f"<textarea readonly style='width:100%; height:40px;'>{sequence}</textarea><br>"
            fasta_content = f">DNA Sequence\n{sequence}"
            html_content += f"<a href='data:text/plain;charset=utf-8,{fasta_content}' download='dna_sequence.fasta'><button>Save as Fasta</button></a><br><br>"

        if 'rnaSequence' in item:
            rna_sequence = item['rnaSequence']
            sequence = rna_sequence['sequence']
            count = rna_sequence['count']
            length = len(sequence)
            html_content += f"<h3>RNA Sequence (Count: {count}, Length: {length})</h3>"
            html_content += f"<textarea readonly style='width:100%; height:40px;'>{sequence}</textarea><br>"
            fasta_content = f">RNA Sequence\n{sequence}"
            html_content += f"<a href='data:text/plain;charset=utf-8,{fasta_content}' download='rna_sequence.fasta'><button>Save as Fasta</button></a><br><br>"

    # 評価指標のセクションを追加
    html_content += "<h2>Model Evaluations</h2>"
    html_content += "<table border='1'><tr><th>Model</th><th>Fraction Disordered</th><th>Has Clash</th><th>iPTM</th><th>PTM</th><th>Ranking Score</th></tr>"

    for i in range(5):
        summary_confidences_path = os.path.join(folder_name, f'{folder_name}_summary_confidences_{i}.json')
        with open(summary_confidences_path, 'r') as f:
            data = json.load(f)

        fraction_disordered = data.get('fraction_disordered', 'N/A')
        has_clash = data.get('has_clash', 'N/A')
        iptm = data.get('iptm', 'N/A')
        ptm = data.get('ptm', 'N/A')
        ranking_score = data.get('ranking_score', 'N/A')

        html_content += f"<tr><td>Model {i}</td><td>{fraction_disordered}</td><td>{has_clash}</td><td>{iptm}</td><td>{ptm}</td><td>{ranking_score}</td></tr>"

    html_content += "</table>"

    # 各ファイルの種類別に見出しを追加
    html_content += "<h2>3D Models of Protein Structure</h2>"
    html_content += f"<h3>chain.gif_structure.gif</h3>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('chain.gif')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('chain_rotated.gif')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('structure.gif')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"

    html_content += "<h2>Protein Contacts</h2>"
    gif_iptm    = [f for f in os.listdir(folder_name) if f.endswith('iptm.gif')]
    gif_pae     = [f for f in os.listdir(folder_name) if f.endswith('pae.gif')]
    gif_contact = [f for f in os.listdir(folder_name) if f.endswith('contact.gif')]
    html_content += f"<h3>iptm.gif_pae.gif_contact.gif</h3>"
    for file in gif_iptm:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in gif_pae:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in gif_contact:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"

    html_content += "<h2>Per-residue Local Confidence</h2>"
    html_content += f"<h3>plddt.gif</h3>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('plddt.gif')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"

    html_content += "<h2>CSV Files</h2>"
    csv_files = [f for f in os.listdir(folder_name) if f.endswith('.csv')]
    for file in csv_files:
        html_content += f"<a href='{file}' download>Download {file}</a><br>"

    html_content += "<h2>Image model 0</h2>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('rank_0.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"
    png_iptm    = [f for f in os.listdir(folder_name) if f.endswith('iptm_0.png')]
    png_pae     = [f for f in os.listdir(folder_name) if f.endswith('pae_0.png')]
    png_contact = [f for f in os.listdir(folder_name) if f.endswith('probs_0.png')]
    for file in png_iptm:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_pae:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_contact:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('plddt_0.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"

    html_content += "<h2>Image model 1</h2>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('rank_1.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"
    png_iptm    = [f for f in os.listdir(folder_name) if f.endswith('iptm_1.png')]
    png_pae     = [f for f in os.listdir(folder_name) if f.endswith('pae_1.png')]
    png_contact = [f for f in os.listdir(folder_name) if f.endswith('probs_1.png')]
    for file in png_iptm:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_pae:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_contact:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('plddt_1.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"

    html_content += "<h2>Image model 2</h2>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('rank_2.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"
    png_iptm    = [f for f in os.listdir(folder_name) if f.endswith('iptm_2.png')]
    png_pae     = [f for f in os.listdir(folder_name) if f.endswith('pae_2.png')]
    png_contact = [f for f in os.listdir(folder_name) if f.endswith('probs_2.png')]
    for file in png_iptm:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_pae:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_contact:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('plddt_2.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"

    html_content += f"<br><br>"

    html_content += "<h2>Image model 3</h2>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('rank_3.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"
    png_iptm    = [f for f in os.listdir(folder_name) if f.endswith('iptm_3.png')]
    png_pae     = [f for f in os.listdir(folder_name) if f.endswith('pae_3.png')]
    png_contact = [f for f in os.listdir(folder_name) if f.endswith('probs_3.png')]
    for file in png_iptm:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_pae:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_contact:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('plddt_3.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"

    html_content += "<h2>Image model 4</h2>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('rank_4.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"
    png_iptm    = [f for f in os.listdir(folder_name) if f.endswith('iptm_4.png')]
    png_pae     = [f for f in os.listdir(folder_name) if f.endswith('pae_4.png')]
    png_contact = [f for f in os.listdir(folder_name) if f.endswith('probs_4.png')]
    for file in png_iptm:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_pae:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'>"
    for file in png_contact:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:100%;'><br><br>"
    image_files = [f for f in os.listdir(folder_name) if f.endswith('plddt_4.png')]
    for file in image_files:
        html_content += f"<img src='{file}' alt='{file}' style='max-width:50%;'>"
    html_content += f"<br><br>"

    html_content += "</body></html>"
    
    # HTMLコンテンツをファイルに書き出す
    report_path = os.path.join(folder_name, "report.html")
    with open(report_path, "w") as f:
        f.write(html_content)
    print(f"HTML report generated at {report_path}")

def process_files(folder_name, rank_by, align_chains=None):
    confidence_data = load_confidence_data(folder_name)
    sorted_data = rank_and_save_cif(confidence_data, folder_name, rank_by)
    align_models(sorted_data, folder_name, align_chains)
    
    image_files_plddt = []
    image_files_structure = []
    image_files_chain = []
    rimage_files_chain = []
    image_files_iptm = []
    image_files_contact = []
    image_files_pae = []

    for rank, data in enumerate(sorted_data, start=1):
        align_cif_file = os.path.join(folder_name, f'{folder_name}_model_{data["model_index"]}_rank_{rank}_align.cif')
        plot_ca_structure_from_cif(align_cif_file, folder_name, rank-1)
        image_files_structure.append(os.path.join(folder_name, f'{folder_name}_structure_rank_{rank-1}.png'))

        plot_chain_colored_structure(align_cif_file, folder_name, rank-1)
        image_files_chain.append(os.path.join(folder_name, f'{folder_name}_chain_rank_{rank-1}.png'))

        plot_chain_colored_structure(align_cif_file, folder_name, rank-1)
        rimage_files_chain.append(os.path.join(folder_name, f'{folder_name}_rotated_chain_rank_{rank-1}.png'))

        num_chains = len(data.get('chain_pair_iptm', []))
        labels = [chr(65 + i) for i in range(num_chains)]
        
        plot_and_save_chain_pair_matrices(data, folder_name, rank-1, labels)
        image_files_iptm.append(os.path.join(folder_name, f'{folder_name}_chain_pair_iptm_{rank-1}.png'))
        image_files_contact.append(os.path.join(folder_name, f'{folder_name}_contact_probs_{rank-1}.png'))

    for i in range(5):
        file_path = os.path.join(folder_name, f'{folder_name}_full_data_{i}.json')
        with open(file_path, 'r') as f:
            data = json.load(f)

        plot_contact_probs(data.get('contact_probs', []), folder_name, i)
        plot_pae(data.get('pae', []), folder_name, i)
        image_files_pae.append(os.path.join(folder_name, f'{folder_name}_pae_{i}.png'))

        cif_file_path = os.path.join(folder_name, f'{folder_name}_model_{i}.cif')
        plot_plddt_from_cif(cif_file_path, folder_name, i)
        image_files_plddt.append(os.path.join(folder_name, f'{folder_name}_cif_plddt_{i}.png'))

    create_gif(image_files_plddt, os.path.join(folder_name, 'plddt.gif'))
    create_gif(image_files_structure, os.path.join(folder_name, 'structure.gif'))
    create_gif(image_files_chain, os.path.join(folder_name, 'chain.gif'))
    create_gif(rimage_files_chain, os.path.join(folder_name, 'chain_rotated.gif'))
    create_gif(image_files_iptm, os.path.join(folder_name, 'chain_pair_iptm.gif'))
    create_gif(image_files_contact, os.path.join(folder_name, 'contact.gif'))
    create_gif(image_files_pae, os.path.join(folder_name, 'pae.gif'))
    
    generate_html_report(folder_name)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process AlphaFold JSON and CIF data.')
    parser.add_argument('-i', '--input', required=True, help='Input folder name containing the data files.')
    parser.add_argument('-r', '--rank_by', default='ranking_score', help='Criteria for ranking (fraction_disordered, iptm, ptm, ranking_score)')
    parser.add_argument('-a', '--align_chains', help='Chains to use for alignment (e.g., "A", "AB")')
    
    args = parser.parse_args()
    align_chains = list(args.align_chains) if args.align_chains else None
    process_job_request(args.input)
    process_files(args.input, args.rank_by, align_chains)
