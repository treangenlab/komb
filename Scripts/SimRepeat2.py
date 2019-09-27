import numpy as np 
import numpy.random as rnd 

# Base code
def generate_random_vec(length, 
                        alphabet_size = 4):
    if length > 0:
        vec = rnd.randint(0, alphabet_size, size = length)
    else:
        vec = []
    return vec

def translate_vec_to_seq(vec, 
                         alphabet = ["A", "G", "T", "C"]):
    seq = []
    for x in vec:
        seq.append(alphabet[x])

    return seq

def generate_with_repeat(length, repeat, repeat_count, 
                         alphabet_size = 4, min_repeat_spacing = 100):
    repeat_len = len(repeat)
    
    repeat_positions = rnd.choice(length // (repeat_len + min_repeat_spacing),
                                  repeat_count, replace = False)
    repeat_positions.sort(kind = "heapsort")
    print(repeat_positions)
    for i in range(len(repeat_positions)):
        repeat_positions[i] = i * (repeat_len + min_repeat_spacing) + repeat_positions[i]

    current_pos = 0
    vec = np.asarray([])
    for pos in repeat_positions:
        vec = np.append(vec, generate_random_vec(pos - current_pos))
        vec = np.append(vec, repeat)
        current_pos += pos + repeat_len

    print(vec)
    vec = np.append(vec, generate_random_vec(length - current_pos))
    return vec.astype(int)

def check_repeat_count(vec, repeats, repeat_counts):
    return

def generate_with_repeats(length, repeats, repeat_counts, 
                          alphabet_size = 4, min_repeat_spacing = 100):
    vec = np.asarray([])
    total_repeats = np.sum(repeat_counts)
    total_repeat_length = np.sum([len(repeats[i]) * repeat_counts[i] 
                                  for i in range(len(repeats))])
    padding_available = length - total_repeat_length
    
    # Generate a vector contatining all repeats needed
    repeat_vec = []
    for i in range(len(repeats)):
        for j in range(repeat_counts[i]):
            repeat_vec.append(repeats[i])
    
    # Randomly permute the repeats
    repeat_vec = rnd.permutation(repeat_vec)
    
    # Generate padding needed
    pad_vec = []
    pad_lengths = rnd.multinomial(padding_available, 
                                  np.ones(total_repeats + 1)/(total_repeats + 1),
                                  size = 1)[0]
    for length in pad_lengths:
        pad_vec.append(generate_random_vec(length))
    
    # Insert random gaps between the repeats to form final genome
    vec = np.append(vec, pad_vec[0])
    for i in range(total_repeats ):
        vec = np.append(vec, repeat_vec[i])
        vec = np.append(vec, pad_vec[i + 1])
    
    # Optional: prune down extra occurences of repeat patterns
    check_repeat_count(vec, repeats, repeat_counts)
    return vec.astype(int)

# Sample tests
def test_seq_generator():
    vec = generate_random_vec(16, 4)
    seq = translate_vec_to_seq(vec, ["A", "G", "T", "C"])
    print(seq)

    vec = generate_random_vec(16, 4)
    seq = translate_vec_to_seq(vec, ["A", "G", "T", "C"])
    print(seq)

    vec = generate_random_vec(24, 4)
    seq = translate_vec_to_seq(vec, ["A", "G", "T", "C"])
    print(seq)

    vec = generate_random_vec(24, 4)
    seq = translate_vec_to_seq(vec, ["A", "G", "T", "C"])
    print(seq)

def test_rep_generator():
    vec = generate_with_repeat(24, [0, 0, 0, 0], 2, min_repeat_spacing = 3)
    seq = translate_vec_to_seq(vec)
    print(seq)

    vec = generate_with_repeat(24, [0, 1, 1, 1], 2, min_repeat_spacing = 3)
    seq = translate_vec_to_seq(vec)
    print(seq)

    vec = generate_with_repeat(24, [0, 0, 0], 3, min_repeat_spacing = 3)
    seq = translate_vec_to_seq(vec)
    print(seq)

    vec = generate_with_repeat(36, [0, 1, 1, 0], 4, min_repeat_spacing = 3)
    seq = translate_vec_to_seq(vec)
    print(seq)

def test_mult_rep_generator():
    vec = generate_with_repeats(50, [[0,0,0], [1,1,1]], [2, 1], min_repeat_spacing = 2)
    seq = translate_vec_to_seq(vec)
    print(seq)

def gen3synthetic():
    inter_gen_repeat = generate_random_vec(500)
    intra_gen_repeat_f1_1 = generate_random_vec(400)
    intra_gen_repeat_f2_1 = generate_random_vec(200)
    intra_gen_repeat_f1_2 = generate_random_vec(400)
    intra_gen_repeat_f2_2 = generate_random_vec(200)
    intra_gen_repeat_f1_3 = generate_random_vec(400)
    intra_gen_repeat_f2_3 = generate_random_vec(200)

    gen_1 = generate_with_repeats(5000000, 
            [inter_gen_repeat, intra_gen_repeat_f1_1, intra_gen_repeat_f2_1],
            [500, 200, 400])
    gen_2 = generate_with_repeats(5000000, 
            [inter_gen_repeat, intra_gen_repeat_f1_2, intra_gen_repeat_f2_2],
            [500, 200, 400])
    gen_3 = generate_with_repeats(5000000, 
            [inter_gen_repeat, intra_gen_repeat_f1_3, intra_gen_repeat_f2_3],
            [500, 200, 400])

    seq_1 = translate_vec_to_seq(gen_1)
    seq_2 = translate_vec_to_seq(gen_2)
    seq_3 = translate_vec_to_seq(gen_3)
    
    with open("repeats.txt", 'w') as f:
        f.write("Seed: " + str(seed) + "\n\n")
        f.write("".join(translate_vec_to_seq(inter_gen_repeat)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(intra_gen_repeat_f1_1)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(intra_gen_repeat_f2_1)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(intra_gen_repeat_f1_2)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(intra_gen_repeat_f2_2)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(intra_gen_repeat_f1_3)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(intra_gen_repeat_f2_3)))
        f.write("\n\n")

    with open("seq1.fa", 'w') as f:
        f.write("".join(seq_1))

    with open("seq2.fa", 'w') as f:
        f.write("".join(seq_2))

    with open("seq3.fa", 'w') as f:
        f.write("".join(seq_3))

def gen1gen1rep(seed):
    repeat1 = generate_random_vec(200)
    repeat2 = generate_random_vec(500)
    repeat3 = generate_random_vec(400)
    repeat4 = generate_random_vec(200)
    repeat5 = generate_random_vec(500)

    gen1 = generate_with_repeats(5000000, [repeat1], [500])
    gen2 = generate_with_repeats(5000000, [repeat2], [500])
    gen3 = generate_with_repeats(5000000, [repeat3], [500])
    gen4 = generate_with_repeats(5000000, [repeat4], [1000])
    gen5 = generate_with_repeats(5000000, [repeat5], [1000])

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)
    seq3 = translate_vec_to_seq(gen3)
    seq4 = translate_vec_to_seq(gen4)
    seq5 = translate_vec_to_seq(gen5)

    with open("seq1a_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.A\n")
        f.write("".join(seq1))

    with open("seq1b_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.B\n")
        f.write("".join(seq2))

    with open("seq1c_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.C\n")
        f.write("".join(seq3))

    with open("seq1d_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.D\n")
        f.write("".join(seq4))

    with open("seq1e_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.E\n")
        f.write("".join(seq5))

    with open("repeats_{}.txt".format(seed), 'w') as f:
        f.write("Seed: " + str(seed) + "\n\n")
        f.write("".join(translate_vec_to_seq(repeat1)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat2)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat3)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat4)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat5)))
        f.write("\n\n")

# Main sub-routine
def main():
    seed = 35432
    rnd.seed(seed)

    gen1gen1rep(seed)

if __name__ == "__main__":
    main()