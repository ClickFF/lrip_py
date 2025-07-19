def read_list(ligs, miss_list_file):
    """
    Read a list of ligands from a file.
    
    Args:
        miss_list_file (str): The file containing the list of ligands.
        
    Returns:
        dict: A dictionary with ligand IDs as keys and their net charges as values.
    """

    with open(miss_list_file, 'r') as f:
        for line in f:
            if len(line.strip()) > 0:
                line = line.strip()
                if line in ligs['id']:
                    idx = ligs['id'].index(line)
                    ligs['id'].pop(idx)
                    ligs['net_charge'].pop(idx)
    return ligs