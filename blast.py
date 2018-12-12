# author: Anthony Ho & Simon Hsu
# run : python blast.py 1 chr22

import random
import sys
import pickle
import time
import os
import matplotlib.pyplot as plt

### test data collections
bps = []
b_time = []
accuracys = []
mutations = []
m_time = []
m_accuracys = []


# constant
NUCLEOTIDE = "ACTG"
MATCH = 1.0
MISSMATCH = -2.0
GAP_PENALTY = -1.0
UNGAPPED_TRESHOLD = 8.0
HSP_TRESHOLD = 100
GAPPED_MULTIPLE = 4
class Seed:
    def __init__(self, ind, prob):
        self.ind = ind
        self.prob = prob

class HSPair:
    def __init__(self, qstart, qend, dstart, dend, score):
        self.qstart = qstart
        self.qend = qend
        self.dstart = dstart
        self.dend = dend
        self.score = score
    
    def __eq__(self, other):
        return self.qstart == other.qstart and self.qend == other.qend and self.dstart == other.dstart and self.dend == other.dend

class Alignment:
    def __init__(self):
        self.query_alignment = None
        self.db_alignment = None
        self.query_start = None
        self.query_end = None
        self.db_start = None
        self.db_end = None
        self.score = None
    
    def __str__(self):
        return "score: " + str(self.score) + " query_pos: [" + str(self.query_start) + ":" + str(self.query_end) + "] db_pos: [" + str(self.db_start) + ":" + str(self.db_end) + "] " #+ self.query_alignment + "\n" + self.db_alignment

class BlastPreprocessor:
    def __init__(self, seq, probs):
        self.sequence = seq
        self.probabilites = probs
        self.seeds = {}

    def preprocess(self):
        if os.path.isfile("seeds.pkl"):
            f = open('seeds.pkl', 'rb')
            self.seeds = pickle.load(f)
            return
        self.seeds = self.find_seeds()
        self.save_seeds()
        return
    
    def find_seeds(self):
        len_seq = len(self.sequence)
        for i in range(len_seq):
            if i + 10 > len_seq:
                break
            wd = "".join(self.sequence[i:i+11])
            probs = sum(self.probabilites[i:i+11])
            if wd not in self.seeds:
                self.seeds[wd] = []
            self.seeds[wd].append(Seed(i, probs))
        return self.seeds
    
    def get_seeds(self):
        return self.seeds

    def save_seeds(self):
        f = open('seeds.pkl', 'wb')
        pickle.dump(self.seeds, f, pickle.HIGHEST_PROTOCOL)
        return

class Blast:
    def __init__(self, seqs, probs, seeds):
        self.database_seq = seqs
        self.database_probs = probs
        self.seeds = seeds
    
    def start(self, seq):
        hsps = self.ungapp(seq)
        # print hsps
        max_score = 0
        max_align = None
        for hsp in hsps:
            aligment = self.gapped(hsp, seq)
            if aligment.score > max_score:
                max_score = aligment.score
                max_align = aligment
        return max_align

    def gapped(self, hsp, query):
        l_score, l_query_seq, l_query_start, l_db_seq, l_db_start = self.align(hsp, query, -1)
        r_score, r_query_seq, r_query_start, r_db_seq, r_db_start = self.align(hsp, query, 1)
        alignment = Alignment()
        alignment.query_alignment = l_query_seq + query[hsp.qstart:hsp.qend+1] + r_query_seq
        alignment.db_alignment = l_db_seq + self.database_seq[hsp.dstart:hsp.dend+1] + r_db_seq
        alignment.score = hsp.score + l_score + r_score
        alignment.query_start = l_query_start
        alignment.query_end = r_query_start
        alignment.db_start = l_db_start
        alignment.db_end = r_db_start

        return alignment

    def align(self, hsp, query, direction):
        best_score = 0
        best_i = 0
        best_j = 0

        if direction == -1:
            query_start = hsp.qstart
            db_start = hsp.dstart
            query_range = query_start + 1
            db_range =  GAPPED_MULTIPLE * (query_start + 1) if db_start + 1 > GAPPED_MULTIPLE*(query_start + 1) else db_start + 1
        else:
            query_start = hsp.qend
            db_start = hsp.dend
            query_range = len(query) - query_start
            db_range = GAPPED_MULTIPLE * (query_range) if (len(self.database_seq) - db_start) > GAPPED_MULTIPLE * query_range else len(self.database_seq) - db_start 


        matrix = [[0.0]*db_range for i in range(query_range)]
        pointer = [[0]*db_range for i in range(query_range)]

        for i in range(query_range):
            for j in range(db_range):
                # initialize gap penality
                if i == 0 and j == 0:
                    matrix[i][j] = 0.0
                elif i == 0:
                    matrix[i][j] = matrix[i][j-1] + GAP_PENALTY
                    pointer[i][j] = 1
                elif j == 0:
                    matrix[i][j] = matrix[i-1][j] + GAP_PENALTY
                    pointer[i][j] = 2
                # build dp matrix
                else:
                    match_score = matrix[i-1][j-1] + self.score_gapped(db_start + direction * j, query[query_start + direction * i])
                    i_gap_score = matrix[i-1][j] + GAP_PENALTY
                    j_gap_score = matrix[i][j-1] + GAP_PENALTY

                    if match_score >= i_gap_score and match_score >= j_gap_score:
                        matrix[i][j] = match_score
                        pointer[i][j] = 0
                    elif j_gap_score > match_score and j_gap_score > i_gap_score:
                        matrix[i][j] = j_gap_score
                        pointer[i][j] = 1
                    else:
                        matrix[i][j] = i_gap_score
                        pointer[i][j] = 2

                    if matrix[i][j] > best_score:
                        best_score, best_i, best_j = matrix[i][j], i, j
        
        # backtracking
        i, j = best_i, best_j
        query_seq = ""
        db_seq = ""
        while i > 0 or j > 0:
            if pointer[i][j] == 0:
                query_seq += query[query_start + direction * i]
                db_seq += self.database_seq[db_start + direction * j]
                i -= 1
                j -= 1
            elif pointer[i][j] == 1:
                query_seq += "-"
                db_seq += self.database_seq[db_start + direction * j]
                j -= 1
            else:
                query_seq += query[query_start + direction * i]
                db_seq += "-"
                i -= 1
        
        return best_score, query_seq, query_start + direction * best_i, db_seq, db_start + direction * best_j

    def ungapp(self, seq):
        seq_len = len(seq)
        hsps = []
        for i in range(seq_len):
            if i + 10 >= seq_len:
                break
            word = seq[i:i+11]
            if word not in self.seeds:
                continue
            seeds = self.seeds[word]

            for seed in seeds:
                # left extension
                score = seed.prob
                max_score = seed.prob
                lctr = 1
                max_lctr = 0
                db_p = seed.ind
                seq_p = i
                while db_p - lctr >= 0 and seq_p - lctr >= 0:
                    score += self.score_ungapped(db_p-lctr, seq[seq_p-lctr])
                    if score > max_score:
                        max_score = score
                        max_lctr = lctr
                    
                    if max_score - score > UNGAPPED_TRESHOLD:
                        break
                    lctr += 1

                # right extension
                score = max_score
                rctr = 1
                max_rctr = 0
                db_p = seed.ind + 10
                seq_p = i + 10
                while db_p + rctr < len(self.database_seq) and seq_p + rctr < seq_len:
                    score += self.score_ungapped(db_p+rctr, seq[seq_p+rctr])
                    if score > max_score:
                        max_score = score
                        max_rctr = rctr
                    
                    if max_score - score > UNGAPPED_TRESHOLD:
                        break
                    rctr += 1
                hsp = HSPair(i-max_lctr, i+10+max_rctr, seed.ind-max_lctr, seed.ind+10+max_rctr, max_score)
                if hsp not in hsps and hsp.score > HSP_TRESHOLD:
                    hsps.append(HSPair(i-max_lctr, i+10+max_rctr, seed.ind-max_lctr, seed.ind+10+max_rctr, max_score))
        return hsps
    
    def score_ungapped(self, ind, nucleotide):
        if nucleotide == self.database_seq[ind]:
            return MATCH * self.database_probs[ind] + MISSMATCH * (1.0 - self.database_probs[ind])
        else:
            return MISSMATCH * (self.database_probs[ind]) + MATCH * ((1.0 - self.database_probs[ind]) / 3.0)
    def score_gapped(self, ind, nucleotide):
        if nucleotide == self.database_seq[ind]:
            return MATCH * self.database_probs[ind] + MISSMATCH * (1.0 - self.database_probs[ind])
        else:
            return MISSMATCH * (self.database_probs[ind]) + MATCH * ((1.0 - self.database_probs[ind]) / 3.0)
    


class Test:
    def __init__(self, iteration, filename):
        self.database_seq = None
        self.database_prob = None
        self.filename = filename
        self.iteration = iteration
        
        #import sequence database
        self.import_sequence(self.filename)
        self.database_len = len(self.database_seq)

        # preprocess database sequence
        self.preprocessor = BlastPreprocessor(self.database_seq, self.database_prob)
        self.preprocessor.preprocess()
        # initialize blast instance
        self.blast = Blast(self.database_seq, self.database_prob, self.preprocessor.get_seeds())

    
    def import_sequence(self, file):
        # import sequence
        f = open(self.filename + ".fa", 'r')
        sequence = ""
        for line in f:
            sequence += line
        self.database_seq = sequence

        # import probs
        f = open(self.filename + '.conf', 'r')
        probs = []
        for line in f:
            ps = map(float, line.split())
            probs += ps
        self.database_prob = probs

    def generate_sequence(self, seq_len):
        # generate sequence based on the probability provided
        start_pos = random.randint(0, self.database_len - seq_len)
        sequence = ""
        for i in range(seq_len):
            p = random.random()
            if p <= self.database_prob[start_pos+i]:
                sequence += self.database_seq[start_pos+i]
            else:
                possible_nucleotide = NUCLEOTIDE.replace(self.database_seq[start_pos+i], "")
                sequence += possible_nucleotide[random.randint(0,2)]
        return start_pos, sequence
    
    def permute_sequence(self, sequence, error_rate):
        # permute it base on error rate
        result = ""
        for s in sequence:
            if random.random() <= error_rate:
                if random.randint(0,1):
                    # delete
                    continue
                else:
                    # insert
                    result += NUCLEOTIDE[random.randint(0,3)]
            else:
                result += s

        return result
    
    def run(self):
        bps_len = [50, 100, 200, 400, 800, 1600, 3200, 5000]
        mutation_rate = [0.00,0.03, 0.05, 0.07, 0.1]

        for i in range(len(bps_len)):
            start_pos, seq = self.generate_sequence(bps_len[i])
            seq_len = len(seq)
            end_pos = start_pos + seq_len
            seq = self.permute_sequence(seq, 0.00)

            start_time = time.time()

            alignment = self.blast.start(seq)
            if alignment == None:
                accuracy = 0.0
            else:
                accuracy = 1.0 * (seq_len - abs(alignment.db_start - start_pos) - abs(alignment.db_end - end_pos))/ seq_len
            diff_time = time.time() - start_time
            bps.append(bps_len[i])
            b_time.append(diff_time)
            accuracys.append(accuracy)

        for i in range(len(mutation_rate)):
            start_pos, seq = self.generate_sequence(4000)
            seq_len = len(seq)
            end_pos = start_pos + seq_len
            seq = self.permute_sequence(seq, mutation_rate[i])

            start_time = time.time()

            alignment = self.blast.start(seq)
            if alignment == None:
                accuracy = 0.0
            else:
                accuracy = 1.0 * (seq_len - abs(alignment.db_start - start_pos) - abs(alignment.db_end - end_pos))/ seq_len
            diff_time = time.time() - start_time
            mutations.append(mutation_rate[i])
            m_time.append(diff_time)
            m_accuracys.append(accuracy)

if __name__ == "__main__":
    # argv 1 iteration argv 2 filename
    testsuit = Test(int(sys.argv[1]), sys.argv[2])
    testsuit.run()

    plt.suptitle('Graph of bps length vs accuracy')
    plt.xlabel('bps length')
    plt.ylabel('accuracy')
    # plt.axis([50, 5000, 0.0, 1.0])
    plt.plot(bps, accuracys)
    plt.show()

    plt.suptitle('Graph of mutation rate vs accuracy')
    plt.xlabel('mutation rate')
    plt.ylabel('accuracy')
    plt.plot(mutations, m_accuracys)
    # plt.axis([0.0, 0.1, 0.0, 1.0])
    plt.show()

    plt.suptitle('Graph of bps length vs time')
    plt.xlabel('bps length')
    plt.ylabel('time(seconds)')
    plt.plot(bps, b_time)
    plt.show()

    plt.suptitle('Graph of mutation rate vs time')
    plt.xlabel('mutation rate')
    plt.ylabel('time (seconds)')
    plt.plot(mutations, m_time)
    plt.show()

