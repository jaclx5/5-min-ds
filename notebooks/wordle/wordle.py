from collections import Counter
from functools import reduce
import copy
import re

import random


class WordleSolver:
    def __init__(self, word_list, verbose=True):
        self.word_list = word_list
        self.verbose = verbose

        self.solution = [""] * 5
        self.contains = [set() for _ in range(5)]
        self.discards = set()

        self.last_word = ""

    def next_try(self, response=None):
        if response == "ggggg":
            return "IT'S DONE !"

        if response:
            self._parse_response(response)
            self._crop_words()

            repeat = True
        else:
            repeat = False

        self.last_word = self._best_candidate_1(repeat=repeat)

        return self.last_word

    def _parse_response(self, resp):
        """
        Update the state variables:
            - solution
            - contains
            - discards
        """
        for i, (c, r) in enumerate(zip(self.last_word, resp)):
            if r == "g":
                # test if the letter in compatible with previous responses
                assert (self.solution[i] in [c, ""]), f"Letter in position {i} must be {self.solution[i]}"
                assert (c not in self.contains[i]) and (c not in self.discards), f"Letter in position {i} can't be {c}"

                self.solution[i] = c

            elif r == "y":
                # test if the letter in compatible with previous response
                assert c not in self.discards, f"Letter in position {i} can't be {c}"

                self.contains[i].add(c)

            else:
                self.discards.add(c)

    def _crop_words(self):
        """
        Update the word list so far
        """

        # compute the pattern expression taking into account the information so far:
        pat = ""
        for (s, c) in zip(self.solution, self.contains):
            if s != "":
                pat += s
            else:
                pat += f"[^{''.join(self.discards) + ''.join(c)}]"

        # keeps the words that match the pattern
        self.word_list = list(filter(lambda w: re.match(pat, w), self.word_list))

        # the misplaced letters must exist in the word
        missplaced = reduce(lambda x, y: x | y, self.contains)

        # keeps the words that contain missplaced letters
        self.word_list = list(filter(lambda w: not (missplaced - set(w)), self.word_list))

    @staticmethod
    def most_common(words, pos):
        max_c, max_i, max_v = "", 0, 0
        
        for i in pos:
            c, v = Counter([w[i] for w in words]).most_common()[0]
            
            if v > max_v:
                max_c, max_i, max_v = c, i, v
                
        return max_c, max_i, max_v

    @staticmethod
    def filter_words(words, c, i, repeat):
        words = list(filter(lambda w: w[i] == c, words))

        if not repeat:
            words = list(filter(lambda w: len(set(w)) == 5, words))

        return words


    def _best_candidate_1(self, repeat):
        """
        Returns the word that contains the most common characters in each position without repeated
        letters
        """

        if self.verbose:
            print("Searching for the best word...")

        # we don't want to mess with the original list
        words = copy.copy(self.word_list)
        
        # list of positions to count    
        pos = list(range(5))

        while pos and len(words) > 1:
            # get the most common positional letter and the position in which it was detected
            (c, i, v) = WordleSolver.most_common(words, pos)
            
            if self.verbose:
                print(f"{c} is the most common letter in position {i+1}. It appears {v:<5d} times in {len(words):<5d} words")
            
            # gets only the words with the MCP letter
            words = WordleSolver.filter_words(words, c, i, repeat)
            
            # remove the detected position        
            del pos[pos.index(i)]

        assert words, "Can't find the word, sorry!"

        return words[0]

    def _best_candidate_2(self, repeat):
        return random.sample(self.word_list, 1)[0]

def get_words(words_file):
    # get all the words in the Ubuntu american-english dictionary
    words = open(words_file).read().split("\n")

    # keep only the 5 letter words
    words = list(filter(lambda w: re.match(r"^[a-z]{5}$", w), words))

    # upper cases
    return list(map(lambda w: w.upper(), words))


#
# Main
#

if __name__ == "__main__":
    print("""
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    When typing WORDLE response use:
        - g - For right letter in the right place (green letters)
        - y - For right letter in the wrong place (yellow letters)
        - x - Wrong letters (gray letters)

    For example:
        If the solution is:   DRINK
        and we try:           FROND
        The response will be: xgxgy

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    """)

    solver = WordleSolver(get_words("../../../data/wordle/words.txt"))

    response = None

    while True:
        word = solver.next_try(response)

        print(f"\nTry the word: {word}\n")

        response = input("Type the WORDLE response: ")

    print("Boss!")