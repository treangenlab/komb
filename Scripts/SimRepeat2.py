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

def fasta_to_vec(filename,
                 alphabet = {"A": 0, "G": 1, "T": 2, "C": 3}):
    vec = np.asarray([])
    i = 0
    seq = ""
    with open(filename, "r") as f:
        for line in f:
            if i > 0:
                seq = seq + line.strip("\n")
            i += 1
    vec = np.append(vec, [alphabet[c] for c in seq])
    return vec


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

def generate_with_repeats_real_genome(genome, repeats, repeat_counts, 
                                      alphabet_size = 4):
    vec = np.asarray([])
    total_repeats = np.sum(repeat_counts)
    total_repeat_length = np.sum([len(repeats[i]) * repeat_counts[i] 
                                  for i in range(len(repeats))])
    padding_available = len(genome) - total_repeat_length

    # Generate a vector contatining all repeats needed
    repeat_vec = []
    for i in range(len(repeats)):
        for j in range(repeat_counts[i]):
            repeat_vec.append(repeats[i])
    
    # Randomly permute the repeats
    repeat_vec = rnd.permutation(repeat_vec)

    # Generate padding lengths
    pad_lengths = rnd.multinomial(padding_available, 
                                  np.ones(total_repeats + 1)/(total_repeats + 1),
                                  size = 1)[0]
    
    # Generate padding out of provided genome and form the final genome
    cur_pos = 0
    vec = np.append(vec, genome[cur_pos:pad_lengths[0]])
    for i in range(total_repeats):
        vec = np.append(vec, repeat_vec[i])
        cur_pos = len(vec)
        vec = np.append(vec, genome[cur_pos:cur_pos+pad_lengths[i + 1]])
    
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

def genecoli(seed):
    repeat1 = generate_random_vec(200)
    repeat2 = generate_random_vec(400)
    repeat3 = generate_random_vec(200)
    repeat4 = generate_random_vec(400)
    repeat5 = generate_random_vec(200)
    repeat6 = generate_random_vec(400)
    repeat7 = generate_random_vec(200)
    repeat8 = generate_random_vec(400)
    repeat9 = generate_random_vec(200)
    repeat10 = generate_random_vec(400)
    repeat11 = generate_random_vec(200)
    inter   = generate_random_vec(500)

    ecoli1 = fasta_to_vec("genome_EscherichiacoliO157strainAR-0427chromosomecompletegenome.fa")
    ecoli2 = fasta_to_vec("genome_EscherichiacoliO157strainAR-0428chromosomecompletegenome.fa")
    ecoli3 = fasta_to_vec("genome_EscherichiacoliO157strainAR-0429chromosomecompletegenome.fa")
    ecoli4 = fasta_to_vec("genome_EscherichiacoliO157strainAR-0430chromosomecompletegenome.fa")
    ecoli5 = fasta_to_vec("genome_EscherichiacoliO157strainFDAARGOS_293chromosomecompletegenome.fa")

    # One repeat family one genome
    gen1 = generate_with_repeats_real_genome(ecoli1, [repeat1], [500])

    # Two intra, one inter, two genomes
    gen2 = generate_with_repeats_real_genome(ecoli2, [repeat2, repeat3, inter], [200, 400, 700])
    gen3 = generate_with_repeats_real_genome(ecoli3, [repeat4, repeat5, inter], [200, 400, 700])
    
    # Two intra, one inter, five genomes
    gen4 = generate_with_repeats_real_genome(ecoli1, [repeat6, repeat7, inter], [200, 400, 700])
    gen5 = generate_with_repeats_real_genome(ecoli4, [repeat8, repeat9, inter], [200, 400, 700])
    gen6 = generate_with_repeats_real_genome(ecoli5, [repeat10, repeat11, inter], [200, 400, 700])

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)
    seq3 = translate_vec_to_seq(gen3)
    seq4 = translate_vec_to_seq(gen4)
    seq5 = translate_vec_to_seq(gen5)
    seq6 = translate_vec_to_seq(gen6)

    with open("seq1_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1\n")
        f.write("".join(seq1))

    with open("seq2_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2\n")
        f.write("".join(seq2))

    with open("seq3_{}.fa".format(seed), "w") as f:
        f.write(">SEQ3\n")
        f.write("".join(seq3))

    with open("seq4_{}.fa".format(seed), "w") as f:
        f.write(">SEQ4\n")
        f.write("".join(seq4))

    with open("seq5_{}.fa".format(seed), "w") as f:
        f.write(">SEQ5\n")
        f.write("".join(seq5))

    with open("seq6_{}.fa".format(seed), "w") as f:
        f.write(">SEQ6\n")
        f.write("".join(seq6))

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
        f.write("".join(translate_vec_to_seq(repeat6)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat7)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat8)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat9)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat10)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat11)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(inter)))
        f.write("\n\n")


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

def gen1gen1rep_long(seed):
    repeat1 = generate_random_vec(1000)
    repeat2 = generate_random_vec(1500)
    
    gen1 = generate_with_repeats(5000000, [repeat1], [100])
    gen2 = generate_with_repeats(5000000, [repeat2], [50])

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)
    
    with open("seq1f_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.A\n")
        f.write("".join(seq1))

    with open("seq1g_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1.B\n")
        f.write("".join(seq2))

    with open("repeats_{}.txt".format(seed), 'w') as f:
        f.write("Seed: " + str(seed) + "\n\n")
        f.write("".join(translate_vec_to_seq(repeat1)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat2)))
        f.write("\n\n")

def gen1gen2rep(seed):
    repeat1 = generate_random_vec(200)
    repeat2 = generate_random_vec(400)

    gen1 = generate_with_repeats(5000000, [repeat1, repeat2], [400, 200])
    
    gen2 = generate_with_repeats(1000000, [repeat1], [400]) 
    gen2 = np.append(gen2, generate_random_vec(3000000))
    gen2 = np.append(gen2, generate_with_repeats(1000000, [repeat2], [200]))

    gen3 = generate_with_repeats(1000000, [repeat2], [200]) 
    gen3 = np.append(gen3, generate_random_vec(3000000))
    gen3 = np.append(gen3, generate_with_repeats(1000000, [repeat1], [400]))

    gen4 = generate_with_repeats(1000000, [repeat1, repeat2], [400, 200]) 
    gen4 = np.append(gen4, generate_random_vec(4000000))

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)
    seq3 = translate_vec_to_seq(gen3)
    seq4 = translate_vec_to_seq(gen4)

    with open("seq2a_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.A\n")
        f.write("".join(seq1))

    with open("seq2b_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.B\n")
        f.write("".join(seq2))

    with open("seq2c_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.C\n")
        f.write("".join(seq3))

    with open("seq2d_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.D\n")
        f.write("".join(seq4))    

    with open("repeats_{}.txt".format(seed), 'w') as f:
        f.write("Seed: " + str(seed) + "\n\n")
        f.write("".join(translate_vec_to_seq(repeat1)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat2)))
        f.write("\n\n")

def gen1gen2rep_long(seed):
    repeat1 = generate_random_vec(1000)
    repeat2 = generate_random_vec(1500)

    gen1 = generate_with_repeats(5000000, [repeat1, repeat2], [200, 100])
    
    gen2 = generate_with_repeats(1000000, [repeat1], [200]) 
    gen2 = np.append(gen2, generate_random_vec(3000000))
    gen2 = np.append(gen2, generate_with_repeats(1000000, [repeat2], [100]))

    gen3 = generate_with_repeats(1000000, [repeat2], [100]) 
    gen3 = np.append(gen3, generate_random_vec(3000000))
    gen3 = np.append(gen3, generate_with_repeats(1000000, [repeat1], [200]))

    gen4 = generate_with_repeats(1000000, [repeat1, repeat2], [200, 100]) 
    gen4 = np.append(gen4, generate_random_vec(4000000))

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)
    seq3 = translate_vec_to_seq(gen3)
    seq4 = translate_vec_to_seq(gen4)

    with open("seq2e_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.E\n")
        f.write("".join(seq1))

    with open("seq2f_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.F\n")
        f.write("".join(seq2))

    with open("seq2g_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.G\n")
        f.write("".join(seq3))

    with open("seq2h_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2.H\n")
        f.write("".join(seq4))    

    with open("repeats_{}.txt".format(seed), 'w') as f:
        f.write("Seed: " + str(seed) + "\n\n")
        f.write("".join(translate_vec_to_seq(repeat1)))
        f.write("\n\n")
        f.write("".join(translate_vec_to_seq(repeat2)))
        f.write("\n\n")


def gen1gen0rep(seed):
    gen1 = generate_random_vec(5000000)
    gen2 = generate_random_vec(5000000)
    gen3 = generate_random_vec(5000000)

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)
    seq3 = translate_vec_to_seq(gen3)
    
    with open("seq0a_{}.fa".format(seed), "w") as f:
        f.write(">SEQ0.A\n")
        f.write("".join(seq1))

    with open("seq0b_{}.fa".format(seed), "w") as f:
        f.write(">SEQ0.B\n")
        f.write("".join(seq2))

    with open("seq0c_{}.fa".format(seed), "w") as f:
        f.write(">SEQ0.C\n")
        f.write("".join(seq3))

def gensmall(seed):
    repeat1 = generate_random_vec(500)

    gen1 = generate_random_vec(2000000)
    gen2 = generate_with_repeats_real_genome(gen1, [repeat1], [10])

    seq1 = translate_vec_to_seq(gen1)
    seq2 = translate_vec_to_seq(gen2)

    with open("seq1_{}.fa".format(seed), "w") as f:
        f.write(">SEQ1\n")
        f.write("".join(seq1))

    with open("seq2_{}.fa".format(seed), "w") as f:
        f.write(">SEQ2\n")
        f.write("".join(seq2))

    with open("repeats_{}.txt".format(seed), 'w') as f:
        f.write("Seed: " + str(seed) + "\n\n")
        f.write("".join(translate_vec_to_seq(repeat1)))
        f.write("\n\n")

# Main sub-routine
def main():
    seed = 0
    rnd.seed(seed)

    # genecoli(seed)
    # gensmall(seed)

    gen3 = generate_random_vec(100000)

    seq3 = translate_vec_to_seq(gen3)

    with open("seq5_{}.fa".format(seed), "w") as f:
        f.write(">SEQ5\n")
        f.write("".join(seq3))

    # gen1gen2rep(seed)
    # gen1gen0rep(seed)
    # gen1gen1rep_long(seed)
    # gen1gen2rep_long(seed)

    # repeat1 = generate_random_vec(500)
    # gen1 = generate_with_repeats(5000000, [repeat1], [500])
    # seq1 = translate_vec_to_seq(gen1)

    # with open("seq1a_{}.fa".format(seed), "w") as f:
    #     f.write(">SEQ1.A\n")
    #     f.write("".join(seq1))

    # with open("repeats_{}.txt".format(seed), 'w') as f:
    #     f.write("Seed: " + str(seed) + "\n\n")
    #     f.write("".join(translate_vec_to_seq(repeat1)))
    #     f.write("\n\n")
    
if __name__ == "__main__":
    main()