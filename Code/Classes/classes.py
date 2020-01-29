import random
from pprint import pprint

class Protein(object):
    def __init__(self, protein_str):
        self.name = protein_str
        self.name_list = list(self.name)
        self.length = len(protein_str) 
        self.route = []
        self.grid = [["0"] * 100 for _ in range(100)]
        self.start_x = 0
        self.start_y = 0
        self.wrong_protein = False
        self.score = 0
        self.bonds = []
        self.neighbours = []