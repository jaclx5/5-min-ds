from collections import Counter

from tqdm import tqdm

from wordle import *

def compute_response(solution, guess):
    resp = ""

    for cw, cg in zip(solution, guess):
        if cw == cg:
            resp += "g"
        elif cg in solution:
            resp += "y"
        else:
            resp += "x"

    return resp



words = get_words("../../../data/wordle/words.txt")
counts = []

for solution in tqdm(words):
    solver = WordleSolver(words, verbose=False)
    response = None

    count = 0
    while response != "ggggg":
        word = solver.next_try(response)
        response = compute_response(solution, word)
        count += 1

    counts.append(count)

x = 0
c = 0
for k, v in Counter(counts).items():
    x += k * v
    c += v
print(x/c)