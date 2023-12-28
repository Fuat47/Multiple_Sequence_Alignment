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

def similarity_score(seq1, seq2):
    matches = 0
    min_length = min(len(seq1[0]), len(seq2[0]))
    for i in range(min_length):
        if seq1[0][i] == seq2[0][i]:
            matches += 1
    return matches / min_length if min_length > 0 else 0

def print_similarity_matrix(similarity_matrix, sequence_names):
    # Determine the width of each cell
    cell_width = max(len(name) for name in sequence_names)
    cell_width = max(cell_width, max(len(str(score)) for row in similarity_matrix for score in row))

    # Print the header
    print(' ' * (cell_width + 1) + ' '.join(name.ljust(cell_width) for name in sequence_names))

    # Print each row
    for i, name in enumerate(sequence_names):
        print(name.ljust(cell_width + 1) + ' '.join(str(score).ljust(cell_width) for score in similarity_matrix[i]))

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

def neighbor_joining(distance_matrix, names, tree='', order=None):
    # Initialize the order array if it's None
    if order is None:
        order = []

    # Base case: if there are only two sequences left, return the tree
    if len(names) == 2:
        return tree + f'({names[0]}:{distance_matrix[0][1] / 2}, {names[1]}:{distance_matrix[0][1] / 2})', order

    # Calculate the Q matrix
    q_matrix = calculate_q_matrix(distance_matrix)

    # Find the pair of sequences with the minimum Q value
    i, j = find_min_q(q_matrix)

    # Add the names to the order array
    order.append((names[i], names[j]))

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
    new_names = [name for k, name in enumerate(names) if k != i and k != j] + ['(' + names[i] + '-' + names[j] + ')']

    # Print the updated distance matrix and names list
    print_similarity_matrix(new_distance_matrix, new_names)
    
    # Update the tree
    tree += f'({names[i]}:{distance_i}, {names[j]}:{distance_j})'
    print(tree)
    # Recurse on the updated distance matrix and names list
    return neighbor_joining(new_distance_matrix, new_names, tree, order)

def align_sequences(sequences, order):
    for i in range(len(order)):
        # Get the current pair of sequences
        first_sequence, second_sequence = order[i]
        multiple_alignment_first = sequences[first_sequence]
        multiple_alignment_second = sequences[second_sequence]
        tupple = ()
        tupple2 = ()
        alignment = global_alignment(tupple + (multiple_alignment_first,),tupple2 + (multiple_alignment_second,), blosum62, gap_penalty)

        # Compute the alignment between the current multiple alignment and the next sequence
        
        # alignment = compute_alignment(multiple_alignment_first, multiple_alignment_second, blosum62, gap_penalty)
        
        # sequences dictionary add new alignments with name order
        # sequences['(' + '-'.join(order[i]) + ')'] = alignment
        
        
        
        
#         # Update the multiple alignment
#         # multiple_alignment = merge_alignments(multiple_alignment, alignment)
#         for element in alignment:
#             print(element)
#     return alignment

# def compute_alignment(seq1, seq2, blosum62, gap_penalty):
#     # Initialize the score matrix and the pointer matrix
#     score_matrix = []
#     for _ in range(len(seq1) + 1):
#         row = []
#         for _ in range(len(seq2) + 1):
#             row.append(0)
#         score_matrix.append(row)

#     pointer_matrix = []
#     for _ in range(len(seq1) + 1):
#         row = []
#         for _ in range(len(seq2) + 1):
#             row.append(0)
#         pointer_matrix.append(row)

#     # Fill the score matrix and the pointer matrix
#     for i in range(1, len(seq1) + 1):
#         for j in range(1, len(seq2) + 1):
#             match = score_matrix[i - 1][j - 1] + blosum62[seq1[i - 1]][seq2[j - 1]]
#             delete = score_matrix[i - 1][j] + gap_penalty
#             insert = score_matrix[i][j - 1] + gap_penalty

#             score_matrix[i][j] = max(match, delete, insert)

#             if score_matrix[i][j] == match:
#                 pointer_matrix[i][j] = 'diag'
#             elif score_matrix[i][j] == delete:
#                 pointer_matrix[i][j] = 'up'
#             else:
#                 pointer_matrix[i][j] = 'left'

#     # Backtrack the pointer matrix to compute the alignment
#     alignment1, alignment2 = '', ''
#     i, j = len(seq1), len(seq2)

#     while i > 0 and j > 0:
#         if pointer_matrix[i][j] == 'diag':
#             alignment1 = seq1[i - 1] + alignment1
#             alignment2 = seq2[j - 1] + alignment2
#             i -= 1
#             j -= 1
#         elif pointer_matrix[i][j] == 'up':
#             alignment1 = seq1[i - 1] + alignment1
#             alignment2 = '-' + alignment2
#             i -= 1
#         else:
#             alignment1 = '-' + alignment1
#             alignment2 = seq2[j - 1] + alignment2
#             j -= 1

#     while i > 0:
#         alignment1 = seq1[i - 1] + alignment1
#         alignment2 = '-' + alignment2
#         i -= 1

#     while j > 0:
#         alignment1 = '-' + alignment1
#         alignment2 = seq2[j - 1] + alignment2
#         j -= 1

#     return alignment1, alignment2


sequences = read_fasta('Input.txt')
blosum62 = read_blosum62('Blosum62.txt')
tupple = ()
for key, value in sequences.items():
    sequences[key] = tupple + (value,)
    tupple = ()

gap_penalty = 2
#int(input("Please enter the gap penalty: "))

names = list(sequences.keys())
similarity_matrix = [[0 for _ in range(len(names))] for _ in range(len(names))]

for i in range(len(names)):
    for j in range(i+1, len(names)):
        alignment_score = global_alignment(sequences[names[i]], sequences[names[j]], blosum62, gap_penalty)[0]
        similarity = similarity_score(sequences[names[i]], sequences[names[j]])
        similarity_matrix[i][j] = similarity
        similarity_matrix[j][i] = similarity


print_similarity_matrix(similarity_matrix, names)


# Build the guide tree
guide_tree, order= neighbor_joining(similarity_matrix, names, '')

# Print the guide tree
print(guide_tree)


