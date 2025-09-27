import numpy as np

def generate_pam(x, input_path, output_path):

    # Order: A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W, Y, V
    frequencies = np.array([
        0.087, 0.041, 0.040, 0.047, 0.033, 0.038, 0.050, 0.089, 0.034, 0.037,
        0.085, 0.081, 0.015, 0.040, 0.051, 0.070, 0.058, 0.010, 0.030, 0.065
    ])
    
    amino_acids = []
    pam1_matrix = []
    
    with open(input_path, 'r') as f:
        lines = f.readlines()
        

        for i, line in enumerate(lines):
            if line.strip().startswith('A') and 'R' in line and 'N' in line:
                amino_acids = line.strip().split()
                break

        for line in lines[i+1:]:
            line = line.strip()
            if line and not line.startswith('#'):
                parts = line.split()
                if len(parts) >= 21: 
                    row_data = [float(x) for x in parts[1:21]]  
                    pam1_matrix.append(row_data)

    pam1_matrix = np.array(pam1_matrix) / 10000.0

    pamx_matrix = np.linalg.matrix_power(pam1_matrix, x)



    
    for i in range(20):
        for j in range(20):
            pamx_matrix[i, j] = np.log10(pamx_matrix[i, j]/frequencies[i])*10    
    
    pamx_matrix = np.round(1 * pamx_matrix).astype(int)
    
    with open(output_path, 'w') as f:
        f.write("   " + " ".join(f"{aa}" for aa in amino_acids) + "\n")
        for i, aa in enumerate(amino_acids):
            f.write(f"{aa} ")
            row_values = " ".join(f"{pamx_matrix[i, j]}" for j in range(20))
            f.write(row_values + "\n")
    
    return pamx_matrix

if __name__ == "__main__":
    path="example/"
    generate_pam(250, path+"mut.txt", path+"pam250.txt")