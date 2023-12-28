def read_fasta(file):
    sequences = {}
    with open(file) as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].rstrip('\n')
                sequences[name] = ''
            else:
                sequences[name] += line.rstrip('\n')
    return sequences

def read_blosum62(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Remove the newline character at the end of each line
    lines = [line.rstrip() for line in lines]

    # The first line contains the amino acids
    amino_acids = lines[0].split()

    # The rest of the lines contain the matrix values
    matrix_values = [line.split()[1:] for line in lines[1:]]

    # Convert the matrix values to integers
    matrix_values = [[int(value) for value in row] for row in matrix_values]

    # Create the BLOSUM62 dictionary
    blosum62 = {amino_acids[i]: {amino_acids[j]: matrix_values[i][j] for j in range(len(amino_acids))} for i in range(len(amino_acids))}
    return blosum62

def global_alignment(seq1, seq2, blosum62, gap_penalty):
    m, n = len(seq1[0]), len(seq2[0])
    similarity_score = 0
    score_matrix = []
    pointer_matrix = []
    for i in range(m+1):
        row = []
        for j in range(n+1):
            row.append(0)
        score_matrix.append(row.copy())
        pointer_matrix.append(row.copy())
    for i in range(m+1):
        score_matrix[i][0] = i * gap_penalty
    for j in range(n+1):
        score_matrix[0][j] = j * gap_penalty
        
    calculate_matrices(m, n, seq1, seq2, score_matrix, pointer_matrix, gap_penalty, blosum62)
    global_alignment = backtracking(seq1, seq2, m, n, pointer_matrix)

    if len(global_alignment) == 2:
        similarity_score = calculate_similarity_score(global_alignment)
    
    return global_alignment, similarity_score

# Computes score and pointer matrices for sequence alignment using given scoring parameters.
def calculate_matrices(m, n, seq1, seq2, score_matrix, pointer_matrix, gap_penalty, blosum62):
    for i in range(1, m+1):
        for j in range(1, n+1): 
            match = score_matrix[i - 1][j - 1]
            delete = score_matrix[i - 1][j]
            insert = score_matrix[i][j - 1]
            for k in range(len(seq1)):
                for l in range(len(seq2)):
                    if (seq1[k][i-1] == '-' or seq2[l][j-1] == '-') and seq1[k][i-1] != seq2[l][j-1]:
                        match += gap_penalty
                    elif seq1[k][i-1] == '-' and seq2[l][j-1] == '-':
                        match += blosum62['*']['*']
                    else:
                        match += blosum62[seq1[k][i-1]][seq2[l][j-1]]
                    if seq1[k][i-1] == '-':
                        delete += blosum62['*']['*']
                    else:
                        delete += gap_penalty
                    if seq2[l][j-1] == '-':
                        insert += blosum62['*']['*']
                    else:
                        insert += gap_penalty
            score_matrix[i][j] = max(match, delete, insert)
            
            if score_matrix[i][j] == match:
                pointer_matrix[i][j] = 'diag'                 
            elif score_matrix[i][j] == delete:
                pointer_matrix[i][j] = 'up'
            else:
                pointer_matrix[i][j] = 'left'

# Backtrack the pointer matrix to compute the alignment
def backtracking(seq1, seq2, m, n, pointer_matrix):
    global_alignment = []
    for _ in range(len(seq1) + len(seq2)):
        global_alignment.append("")

    while m > 0 or n > 0:
        if pointer_matrix[m][n] == 'diag':
            for k in range(len(seq1)):
                global_alignment[k] = seq1[k][m-1] + global_alignment[k]
            for l in range(len(seq2)):
                global_alignment[len(seq1)+l] = seq2[l][n-1] + global_alignment[len(seq1)+l]
            m -= 1
            n -= 1
        elif n == 0 or pointer_matrix[m][n] == 'up':
            for k in range(len(seq1)):
                global_alignment[k] = seq1[k][m-1] + global_alignment[k]
            for l in range(len(seq2)):
                global_alignment[len(seq1)+l] = '-' + global_alignment[len(seq1)+l]
            m -= 1
        elif m == 0 or pointer_matrix[m][n] == 'left':
            for k in range(len(seq1)):
                global_alignment[k] = '-' + global_alignment[k]
            for l in range(len(seq2)):
                global_alignment[len(seq1)+l] = seq2[l][n-1] + global_alignment[len(seq1)+l]
            n -= 1

    return global_alignment

def calculate_similarity_score(Profile):
    matches = 0
    for i in range(len(Profile[0])):
        if Profile[0][i] == Profile[1][i]:
            matches += 1
    return matches / len(Profile[0])

def generate_similarity_matrix(sequences, names):
    similarity_matrix = [[0 for _ in range(len(names))] for _ in range(len(names))]

    for i in range(len(names)):
        for j in range(len(names)):
            seq1 = [sequences[names[i]]]
            seq2 = [sequences[names[j]]]
            similarity_score = global_alignment(seq1, seq2, blosum62, gap_penalty)[1]
            similarity_matrix[i][j] = similarity_score
    print_similarity_matrix(similarity_matrix, names)
    return similarity_matrix

def print_similarity_matrix(similarity_matrix, sequence_names):
    print()
    # Determine the width of each cell
    cell_width = max(len(name) for name in sequence_names)
    cell_width = max(cell_width, max(len(str(score)) for row in similarity_matrix for score in row))

    # Print the header
    print(' ' * (cell_width + 1) + ' '.join(name.ljust(cell_width) for name in sequence_names))

    # Print each row
    for i, name in enumerate(sequence_names):
        print(name.ljust(cell_width + 1) + ' '.join(str(score).ljust(cell_width) for score in similarity_matrix[i]))
        
def upgma(similarity_matrix, sequence_names):
    # Find the max similarity score
    max_score = 0
    for i in range(len(similarity_matrix)):
        for j in range(len(similarity_matrix)):
            if i != j and similarity_matrix[i][j] > max_score:
                max_score = similarity_matrix[i][j]
                max_i = i
                max_j = j     
    
    # Combine the two most similar sequences 
    new_sequence = '(' + sequence_names[max_i] + '-' + sequence_names[max_j] + ')'
    
    # Update the sequence names
    sequence_names.append(new_sequence)
    
    # Update the similarity matrix
    new_similarity_matrix = [[0 for _ in range(len(sequence_names))] for _ in range(len(sequence_names))]
    for i in range(len(sequence_names)):
        for j in range(len(sequence_names)):
            if i == j:
                new_similarity_matrix[i][j] = 1.0
            elif i == len(sequence_names) - 1:
                new_similarity_matrix[i][j] = (similarity_matrix[max_i][j] + similarity_matrix[max_j][j]) / 2
            elif j == len(sequence_names) - 1:
                new_similarity_matrix[i][j] = (similarity_matrix[i][max_i] + similarity_matrix[i][max_j]) / 2
            else:
                new_similarity_matrix[i][j] = similarity_matrix[i][j]
    
    # Remove the two most similar sequences from the similarity matrix
    # Ensure max_i is the larger index
    if max_j > max_i:
        max_i, max_j = max_j, max_i

    # Remove the elements at max_i and max_j
    new_similarity_matrix.pop(max_i)
    new_similarity_matrix.pop(max_j)
    for row in new_similarity_matrix:
        row.pop(max_i)
        row.pop(max_j)
    
    # Remove the names of the two most similar sequences
    sequence_names.pop(max_i)
    sequence_names.pop(max_j)
    
    # Print the new similarity matrix
    print_similarity_matrix(new_similarity_matrix, sequence_names)
    return len(sequence_names) != 2, new_similarity_matrix, 1 - max_score

def perform_alignment(distance_dict, sequences, blosum62, gap_penalty):
    alignment_dict = {}
    for key in distance_dict:
        sequence_list = (key.replace('(', '').replace(')', '')).split("-")
        if len(sequence_list) == 2:
            alignment_dict[key] = global_alignment([sequences[sequence_list[0]]], [sequences[sequence_list[1]]], blosum62, gap_penalty)[0]
            print_group(key, alignment_dict)
        else:
            profiles = []
            removed_sequence = ""
            for key2 in alignment_dict:
                if key2 in key:
                    profiles.append(alignment_dict[key2])
                    removed_sequence = key2
                    keys = key2.replace('(', '').replace(')', '').split("-")
                    for key3 in keys:
                        if key3 in sequence_list:
                            sequence_list.remove(key3)
            alignment_dict.pop(removed_sequence)
            if len(sequence_list) == 1:
                alignment_dict[key] = global_alignment(profiles[0], [sequences[sequence_list[0]]], blosum62, gap_penalty)[0]
                print_group(key, alignment_dict)
            else:
                alignment_dict[key] = global_alignment(profiles[0], profiles[1], blosum62, gap_penalty)[0]
                print_group(key, alignment_dict)
    return alignment_dict

def print_group(key, alignment_dict):
    print("Group   : " + key)
    for value in alignment_dict[key]:
        matched_sequence = [name for name, sequence in sequences.items() if sequence == value.replace('-', '')]
        if matched_sequence:
            print(f"{matched_sequence[0].ljust(8)}: {value}")
    print()
    
sequences = read_fasta('Input.txt')
blosum62 = read_blosum62('Blosum62.txt')
gap_penalty = -int(input("Please enter the gap penalty: "))
names = list(sequences.keys())
similarity_matrix = generate_similarity_matrix(sequences, names)
distance_dict = {}
flag = True
while flag:
    flag, similarity_matrix, distance = upgma(similarity_matrix, names)
    distance_dict[names[-1]] = distance
distance_dict['(' + names[0] + '-' + names[1] + ')'] = similarity_matrix[0][1]
print("\nLast Group: " + '(' + names[0] + '-' + names[1] + ')\n')

alignment_dict = perform_alignment(distance_dict, sequences, blosum62, gap_penalty)


