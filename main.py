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
    #print_blosum62(blosum62)
    return blosum62

def print_blosum62(blosum62):
    # Get the amino acids
    amino_acids = list(blosum62.keys())

    # Print the header
    print(''.join(amino_acids))

    # Print each row
    for amino_acid1 in amino_acids:
        row = [blosum62[amino_acid1][amino_acid2] for amino_acid2 in amino_acids]
        print(amino_acid1 + ''.join(map(str, row)))

def global_alignment(seq1, seq2, blosum62, gap_penalty):
    m, n = len(seq1), len(seq2)
    score = [[0 for _ in range(n+1)] for _ in range(m+1)]
    for i in range(m+1):
        score[i][0] = -i * gap_penalty
    for j in range(n+1):
        score[0][j] = -j * gap_penalty
    for i in range(1, m+1):
        for j in range(1, n+1):
            match = score[i-1][j-1] + blosum62[seq1[i-1]][seq2[j-1]]
            delete = score[i-1][j] - gap_penalty
            insert = score[i][j-1] - gap_penalty
            score[i][j] = max(match, delete, insert)
    return score[m][n]

def similarity_score(seq1, seq2, alignment_score):
    matches = sum([1 for a, b in zip(seq1, seq2) if a == b])
    return matches / len(seq1) if len(seq1) > 0 else 0

def print_similarity_matrix(similarity_matrix, sequence_names):
    # Determine the width of each cell
    cell_width = max(len(name) for name in sequence_names)
    cell_width = max(cell_width, max(len(str(score)) for row in similarity_matrix for score in row))

    # Print the header
    print(' ' * (cell_width + 1) + ' '.join(name.ljust(cell_width) for name in sequence_names))

    # Print each row
    for i, name in enumerate(sequence_names):
        print(name.ljust(cell_width + 1) + ' '.join(str(score).ljust(cell_width) for score in similarity_matrix[i]))
     
def print_alignment_score(name1, name2, score):
    # Determine the width of each column
    column_widths = [
        len(name1),  # Sequence 1
        len(name2),  # Sequence 2
        len(str(score))  # Score
    ]

    # Print the header
    print('Sequence 1'.ljust(column_widths[0] + 1) + 'Sequence 2'.ljust(column_widths[1] + 1) + 'Score'.ljust(column_widths[2]))

    # Print the row
    print(name1.ljust(column_widths[0] + 1) + name2.ljust(column_widths[1] + 1) + str(score).ljust(column_widths[2]))

def calculate_q_matrix(distance_matrix):
    # Calculate the Q matrix from the distance matrix
    q_matrix = [[0 for _ in range(len(distance_matrix))] for _ in range(len(distance_matrix))]
    for i in range(len(distance_matrix)):
        for j in range(len(distance_matrix)):
            if i != j:
                q_matrix[i][j] = (len(distance_matrix) - 2) * distance_matrix[i][j] - sum(distance_matrix[i]) - sum(distance_matrix[j])
    return q_matrix

def find_min_q(q_matrix):
    # Find the pair of sequences with the minimum Q value
    min_i, min_j = 0, 1
    min_q = q_matrix[0][1]
    for i in range(len(q_matrix)):
        for j in range(i + 1, len(q_matrix)):
            if q_matrix[i][j] < min_q:
                min_q = q_matrix[i][j]
                min_i, min_j = i, j
    return min_i, min_j

def neighbor_joining(distance_matrix, names, tree=''):
    # Base case: if there are only two sequences left, return the tree
    if len(names) == 2:
        return tree + f'({names[0]}:{distance_matrix[0][1] / 2}, {names[1]}:{distance_matrix[0][1] / 2})'

    # Calculate the Q matrix
    q_matrix = calculate_q_matrix(distance_matrix)

    # Find the pair of sequences with the minimum Q value
    i, j = find_min_q(q_matrix)

    # Calculate the distances from the new node to the sequences i and j
    distance_i = 0.5 * distance_matrix[i][j] + 0.5 * (sum(distance_matrix[i]) - sum(distance_matrix[j])) / (len(distance_matrix) - 2)
    distance_j = distance_matrix[i][j] - distance_i

    # Create the new node
    new_node = {names[i]: distance_i, names[j]: distance_j}

    # Update the distance matrix and the names list
    new_distances = []
    for k in range(len(distance_matrix)):
        if k != i and k != j:
            distance = (distance_matrix[k][i] + distance_matrix[k][j] - distance_matrix[i][j]) / 2
            new_distances.append(distance)
    new_distances.append(0)
    new_distance_matrix = []
    for k in range(len(distance_matrix)):
        if k != i and k != j:
            row = []
            for l in range(len(distance_matrix)):
                if l != i and l != j:
                    row.append(distance_matrix[k][l])
            new_distance_matrix.append(row)
    new_distance_matrix.append(new_distances)
    for k in range(len(new_distance_matrix) - 1):
        new_distance_matrix[k].append(new_distances[k])
    new_names = [name for k, name in enumerate(names) if k != i and k != j] + [names[i] + '-' + names[j]]

    # Print the updated distance matrix and names list
    print_similarity_matrix(new_distance_matrix, new_names)
    
    # Update the tree
    tree += f'({names[i]}:{distance_i}, {names[j]}:{distance_j})'
    print(tree)
    # Recurse on the updated distance matrix and names list
    return neighbor_joining(new_distance_matrix, new_names, tree)

sequences = read_fasta('Input.txt')
blosum62 = read_blosum62('Blosum62.txt')

gap_penalty = int(input("Please enter the gap penalty: "))

names = list(sequences.keys())
similarity_matrix = [[0 for _ in range(len(names))] for _ in range(len(names))]

for i in range(len(names)):
    for j in range(i+1, len(names)):
        alignment_score = global_alignment(sequences[names[i]], sequences[names[j]], blosum62, gap_penalty)
        print_alignment_score(names[i], names[j], alignment_score)
        similarity = similarity_score(sequences[names[i]], sequences[names[j]], alignment_score)
        similarity_matrix[i][j] = similarity
        similarity_matrix[j][i] = similarity


print_similarity_matrix(similarity_matrix, names)


# Build the guide tree
guide_tree = neighbor_joining(similarity_matrix, names, '')

# Print the guide tree
print(guide_tree)