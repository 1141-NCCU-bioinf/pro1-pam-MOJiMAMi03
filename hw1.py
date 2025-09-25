import numpy as np
import pandas as pd

def generate_pam(x, input_path, output_path):
    """
    Generate PAMx matrix from mutation probability matrix.
    
    Args:
        x (int): The PAM distance (e.g., 250 for PAM250)
        input_path (str): Path to the input mutation probability matrix file
        output_path (str): Path to save the output PAMx matrix
    """
    
    # Read the mutation probability matrix
    # Assuming the input file has amino acids as headers and index
    df = pd.read_csv(input_path, sep='\t', index_col=0)
    
    # Convert to numpy array for matrix operations
    M = df.values.astype(float)
    
    # Get amino acid labels
    amino_acids = df.columns.tolist()
    
    # Ensure the matrix is square
    assert M.shape[0] == M.shape[1], "Input matrix must be square"
    
    # Calculate M^x (matrix power)
    # This gives us the probability matrix for x evolutionary time units
    Mx = np.linalg.matrix_power(M, x)
    
    # Calculate amino acid frequencies (assuming equal frequencies as starting point)
    # In practice, these might be given or calculated from a database
    # Using equal frequencies for simplicity: 1/20 for each amino acid
    frequencies = np.ones(len(amino_acids)) / len(amino_acids)
    
    # Alternative: If frequencies are provided in the data, extract them
    # For now, we'll use uniform distribution
    
    # Calculate the PAM matrix as log-odds scores
    # PAM(i,j) = log(M^x(i,j) / f_j) where f_j is frequency of amino acid j
    pam_matrix = np.zeros_like(Mx)
    
    for i in range(len(amino_acids)):
        for j in range(len(amino_acids)):
            if frequencies[j] > 0 and Mx[i, j] > 0:
                # Log-odds score: log(observed/expected)
                pam_matrix[i, j] = np.log(Mx[i, j] / frequencies[j])
            else:
                # Handle zero probabilities with a very negative score
                pam_matrix[i, j] = -10
    
    # Convert to log base 10 and scale (common practice is to multiply by 10)
    pam_matrix = pam_matrix / np.log(10) * 10
    
    # Round to integers
    pam_matrix = np.round(pam_matrix).astype(int)
    
    # Create DataFrame for output
    pam_df = pd.DataFrame(pam_matrix, index=amino_acids, columns=amino_acids)
    
    # Save to output file
    pam_df.to_csv(output_path, sep='\t')
    
    print(f"PAM{x} matrix generated and saved to {output_path}")
    return pam_df
