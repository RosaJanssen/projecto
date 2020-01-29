from queue import Queue
import copy
from ..Classes.classes import Protein
from math import floor
import matplotlib.pyplot as plt
import numpy as np
import random
import csv

def check_neighbours_h(protein, pos_x, pos_y):
    """
    Given a specific position in a given protein, this function checks whether the neighbouring aminoacids are hydrophobic. 
    It returns a list that includes the corresponding directions to arrive at the hydrophobic neighbours. 
    """
    
    # Initialize an empty list 
    neighbours = []

    # If upper neighbour is hydrophobic, add corresponding direction to the neighbours list 
    if protein.grid[pos_y - 1][pos_x] == "H":
        neighbours.append(2)

    # If lower neighbour is hydrophobic, add corresponding direction to the neighbours list                     
    if protein.grid[pos_y + 1][pos_x] == "H":
        neighbours.append(-2)

    # If left neighbour is hydrophobic, add corresponding direction to the neighbours list             
    if protein.grid[pos_y][pos_x - 1] == "H":
        neighbours.append(-1)

    # If right neighbour is hydrophobic, add corresponding direction to the neighbours list             
    if protein.grid[pos_y][pos_x + 1] == "H":
        neighbours.append(1)

    return neighbours

def check_neighbours_c(protein, pos_x, pos_y):
    """
    Given a specific position in a given protein, this function checks whether the neighbouring aminoacids are cysteine. 
    It returns a list that includes the corresponding directions to arrive at the cysteine neighbours. 
    """
    
    # Initialize an empty list 
    neighbours = []

    # If upper neighbour is cysteine, add corresponding direction to the neighbours list       
    if protein.grid[pos_y - 1][pos_x] == "C":
        neighbours.append(2)

    # If lower neighbour is cysteine, add corresponding direction to the neighbours list       
    if protein.grid[pos_y + 1][pos_x] == "C":
        neighbours.append(-2)

    # If left neighbour is cysteine, add corresponding direction to the neighbours list
    if protein.grid[pos_y][pos_x - 1] == "C":
        neighbours.append(-1)

    # If right neighbour is cysteine, add corresponding direction to the neighbours list        
    if protein.grid[pos_y][pos_x + 1] == "C":
        neighbours.append(1)

    return neighbours

def calculate_score(protein):
    """
    This function calculates the score of a protein. It loops over the protein length to first find HH- and HC-bonds.
    Subsequently, it looks for CC-bonds. The corresponding points are added to the score. For HH- and CC-bonds the 
    addition to the score is halved, in order to prevent duplicate bonds.
    """

    # Define x and y coordinates as starting positions (0,0; defined in Protein class)
    pos_x = protein.start_x
    pos_y = protein.start_y

    # Loop over the protein length
    for x in range(protein.length):

        # Check whether the aminoacid is hydrophobic              
        if protein.name_list[x] == "H":

            # If the aminoacid is hydrophobic, check the neighbouring hydrophobic aminoacids
            neighbours_h = check_neighbours_h(protein = protein, pos_x = pos_x, pos_y = pos_y)
            
            # Loop over the list of (directions to) hydrophobic neighbours 
            for i in neighbours_h:
                
                # If the first aminoacid is not connected to the hydrophobic neighbour, add 0.5 to the score 
                if x == 0:
                    if protein.route[x] != i:
                        protein.score += 0.5

                # If the final aminoacid is not connected to the hydrophobic neighbour, add 0.5 to the score
                elif x >= len(protein.route) - 1:
                    if protein.route[x-1] != - i :
                        protein.score += 0.5                  
                
                # If any aminoacid (except the first and final) is not connected to the hydrophobic neighbour, add 0.5 to the score 
                elif protein.route[x] != i and protein.route[x-1] != -i:
                    protein.score += 0.5
                
            # If the aminoacid is hydrophobic, check the neighbouring cysteine aminoacids
            neighbours_c = check_neighbours_c(protein = protein, pos_x = pos_x, pos_y = pos_y)
            
            # Loop over the list of (directions to) cysteine neighbours 
            for i in neighbours_c:
                
                # If the first aminoacid is not connected to the cysteine neighbour, add 1.0 to the score 
                if x == 0:
                    if protein.route[x] != i:
                        protein.score += 1
                
                # If the final aminoacid is not connected to the cysteine neighbour, add 1.0 to the score
                elif x >= len(protein.route) - 1:
                    if protein.route[x-1] != - i :
                        protein.score += 1                  
                
                # If any aminoacid (but the first and final) is not connected to the cysteine neighbour, add 1.0 to the score 
                elif protein.route[x] != i and protein.route[x-1] != -i:
                    protein.score += 1

        # Check whether the aminoacid is a cysteine aminoacid
        if protein.name_list[x] == "C":
            
            # If the aminoacid is cysteine, check the neighbouring cysteine aminoacids 
            # Neighbouring hydrophobic aminoacids are not checked as HC/CH-bonds are already found in the previous part (line 173 - 191)
            neighbours_c = check_neighbours_c(protein = protein, pos_x = pos_x, pos_y = pos_y)
            
            # Loop over the list of (directions to) cysteine neighbours 
            for i in neighbours_c:

                # If the first aminoacid is not connected to the cysteine neighbour, add 1.0 to the score 
                if x == 0:
                    if protein.route[x] != i:
                        protein.score += 2.5
                      
                # If the final aminoacid is not connected to the cysteine neighbour, add 1.0 to the score  
                elif x >= len(protein.route) - 1:
                    if protein.route[x-1] != - i :
                        protein.score += 2.5                 

                # If any aminoacid (but the first and final) is not connected to the cysteine neighbour, add 1.0 to the score 
                elif protein.route[x] != i and protein.route[x-1] != -i:
                    protein.score += 2.5

        # As long as x is within the range of the route length, 
        if x < len(protein.route):
            
            # If the move is horizontal, adapt the x-coordinate 
            if abs(protein.route[x]) <= 1:
                pos_x += protein.route[x]

            # If the move is vertical, adapt the y-coordinate
            else:
                pos_y -= int(protein.route[x]/2)


def create_queue(constant_string, protein_string):
    """
    This function creates a queue varying all possible directions the protein could fold into, given the protein name. 
    Returns the queue with the complete routes. 
    """

    # Define the desired length of the route 
    depth = len(protein_string) - 1 

    # Create queue
    queue = Queue()

    # Initialize empty list for complete routes
    queue2 = []

    # Start queue with empty string
    queue.put("")

    # Checks if there are elements left in queue
    while not queue.empty():
        
        # Get first element from queue
        state = queue.get()
        
        # Append direction to route as long as the depth is not reached
        if len(state) < depth:
            
            # Iterate over possible directions
            for i in ['A', 'B', 'C', 'D']:
                child = copy.deepcopy(state)
                child += i
                queue.put(child)
                
                # Append to list if route is of desired length
                if len(state) == depth - 1:
                    queue2.append(child)

    return queue2

def make_proteins(path, substring, queue2):
    """
    This function loops through the queue with the proteins found and returns the protein with the highest stability.
    """   

    # Initialize integer to save highest score
    highest_score = 0

    # Initialize empty list to append best proteins to
    all_proteins = []

    # Initialize empty list to append other proteins to
    best_proteins = []

    # Initialize best grid to none
    best_grid = None
    
    # Loop through all proteins
    for i in queue2:

        # Initialize empty list to save coordinates
        coordinates = []    

        # Make protein
        protein = Protein(path + substring)

        # Initialize starting position
        pos_x = protein.start_x
        pos_y = protein.start_y

        # Append first coordinate
        coordinates.append((pos_x, pos_y))

        # Make list variable to iterate over route
        lijst = list(i)

        # Change starting position to first amino acid name
        protein.grid[pos_y][pos_x] = protein.name_list[0]
        
        # Iterate over route
        for j in range(protein.length - 1):
            
            # Find movement
            move = lijst[j]

            # If move is upwards, change y coordinate and append to route
            if move == "A":
                pos_y -= 1
                protein.route.append(2)

            # If move is downwards, change y coordinate and append to route
            if move == "C":
                pos_y += 1
                protein.route.append(-2)

            # If move is to the left, change x coordinate and append to route
            if move == "D":
                pos_x -= 1
                protein.route.append(-1)

            # If first move is to the right, change x coordinate and append to route
            if move == "B":
                pos_x += 1
                protein.route.append(1)

            # Change next position in list to next amino acid element 
            protein.grid[pos_y][pos_x] = protein.name_list[j+1]
            
            # If the coordinate is not taken yet, append to list with coordinates
            if (pos_x, pos_y) not in coordinates:
                coordinates.append((pos_x, pos_y))

        # Check if there are no coordinates double in use to prevent invalid proteins
        if len(coordinates) == protein.length:
            
            # Calculate stability of protein
            calculate_score(protein = protein)

            # Check if stability is higher then the ones found yet
            if protein.score > highest_score:
                
                # If so, save protein
                best_proteins = []
                highest_score = protein.score
                best_proteins.append(protein)

            # Check if stability is as high as the ones found, if so: save to list
            elif protein.score == highest_score  :
                best_proteins.append(protein)
            
            # If not, save protein to list of other proteins
            else:
                all_proteins.append(protein)

    # If there are proteins with a higher score than 0, choose one randomly 
    if best_proteins != []:
        best_protein = random.choice(best_proteins)

    # If not, choose one random other protein
    else:
        if all_proteins == []:
            return 0 
        best_protein = random.choice(all_proteins)

    return(best_protein)

def generate_data(times, protein_string):
    """
    This function creates all possible routes for a protein x times (indicated by the user) and 
    writes the results (protein route and score) in the result file. This file can be used to create histograms.
    """ 
    best_protein = None
    # Loop over number of times to run the algoritm (indicated by user)
    for a in range(times):

        # Set look ahead variable to four aminoacids
        lookahead = 4

        # Initialize path 
        path = ""

        # Initialize constant string: so far fixed part of the route
        constant_string = ""

        # Loop over the protein characters
        for i in range(len(protein_string)-1):

            # Define substring as final four characters
            substring = protein_string[i:i+lookahead]

            # Create all possible routes
            queue = create_queue(constant_string, substring)
            
            # Add the previous steps to the route
            for j in range(len(queue)):
                queue[j] = constant_string + queue[j]

            # Find protein with highest stability 
            best_protein = make_proteins(path, substring, queue)
            
            
            # If protein is invalid, skip to next iteration
            if best_protein == 0:
                break

            # Define next step in route
            else:
                if best_protein.route[i] == 1:
                    add = "B"
                if best_protein.route[i] == -1:
                    add = "D"
                if best_protein.route[i] == 2:
                    add = "A"
                if best_protein.route[i] == -2:
                    add = "C"

                # Add step to definite route
                constant_string += add  
                path += str(best_protein.name[i])

        if best_protein != 0 and best_protein != None:
            # Print the best score
            print("hoogste score", best_protein.score)
            
            # Write best routes and scores in the result file 
            f = open("resultfile", "a")
            with f:
                writer = csv.writer(f)
                writer.writerow(["Lookahead", best_protein.route, best_protein.score])


# if __name__ == "__main__":

#     # Open result file and initialize rows for the algoritm, route and stability 
#     f = open("resultfile", "w")
#     with f:
#         writer = csv.writer(f)
#         writer.writerow(["Algoritme", "Route", "Stability"])
    
#     # Ask user for the protein structure
#     protein_string = input("Please give the protein structure: ")

#     # Ask user for the number of times to run algorithm 
#     times = input("Please indicate the number of times to repeat the algoritm: ")

#     # Convert number of times to an integer
#     times_number = int(times)

#     # Call the generate data function and pass along the number of times to run the algoritm and the protein structure
#     generate_data(times = times_number, protein_string = protein_string)