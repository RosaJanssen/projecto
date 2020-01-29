from ..Classes.classes import Protein
import random
from pprint import pprint
import matplotlib.pyplot as plt
import csv
import numpy as np

def create_random(protein):
    """
    This function creates a random protein route, given the protein structure (protein.name_list) and the start coordinates. 
    It first places the first aminoacid in the middle of a grid. Subsequently, the algorithm checks 
    whether the neighbours of the aminoacid (up, down, left and right) in the grid are free. 
    Finally, the algorithm randomly chooses one of the free neighbours to place the succeeding
    aminoacid. The latter process is repeated until all aminoacids from the protein structure have been placed.
    """

    # Initialize the grid by placing the first aminoacid in the middle of the grid 
    protein.grid[protein.start_y][protein.start_x] = protein.name_list[0]
    
    # Define x and y coordinates as starting positions (0,0; defined in Protein class)
    pos_x = protein.start_x
    pos_y = protein.start_y

    # Loop over the length of the protein name list 
    for x in range(len(protein.name_list) - 1):
        
        # Initialize empty list for available neighbouring places
        neighbours = []

        # If upper neighbour is available (value = 0) add "above" the the neighbour list 
        if protein.grid[pos_y - 1][pos_x] == "0":
            neighbours.append("above")
        
        # If lower neighbour is available (value = 0) add "down" the the neighbour list 
        if protein.grid[pos_y + 1][pos_x] == "0":
            neighbours.append("down")
        
        # If left neighbour is available (value = 0) add "left" the the neighbour list 
        if protein.grid[pos_y][pos_x - 1] == "0":
            neighbours.append("left")

        # If right neighbour is available (value = 0) add "right" the the neighbour list         
        if protein.grid[pos_y][pos_x + 1] == "0":
            neighbours.append("right")

        # If there is no neighbouring spot available at all, define the protein object to be a wrong protein 
        if neighbours == []:
            protein.wrong_protein = True
            return

        # Randomly choose one available neighbour from all available neighbours       
        move = random.choice(neighbours)

        # If upper neighbour is chosen, update y-coordinate and add corresponding direction to the protein route 
        if move == "above":
            pos_y -= 1
            protein.route.append(2)

        # If lower neighbour is chosen, update y-coordinate and add corresponding direction to the protein route 
        if move == "down":
            pos_y += 1
            protein.route.append(-2)

        # If left neighbour is chosen, update x-coordinate and add corresponding direction to the protein route 
        if move == "left":
            pos_x -= 1
            protein.route.append(-1)

        # If right neighbour is chosen, update x-coordinate and add corresponding direction to the protein route 
        if move == "right":
            pos_x += 1
            protein.route.append(1)

        # Place succeeding aminoacid at the updated x- and y-coordinate
        protein.grid[pos_y][pos_x] = protein.name_list[x+1]



def check_neighbours_h(protein, pos_x, pos_y):
    """
    Given a specific position in a given protein, this function checks whether one or more of the neighbouring aminoacids are hydrophobic. 
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
    Given a specific position in a given protein, this function checks whether one or more of the neighbouring aminoacids 
    are cysteine aminoacids. It returns a list that includes the corresponding directions to arrive at the cysteine neighbours. 
    """

    # Initialize an empty list 
    neighbours = []

    # If upper neighbour is a cysteine aminoacid, add corresponding direction to the neighbours list 
    if protein.grid[pos_y - 1][pos_x] == "C":
        neighbours.append(2)
    
    # If lower neighbour is a cysteine aminoacid, add corresponding direction to the neighbours list 
    if protein.grid[pos_y + 1][pos_x] == "C":
        neighbours.append(-2)
    
    # If left neighbour is a cysteine aminoacid, add corresponding direction to the neighbours list 
    if protein.grid[pos_y][pos_x - 1] == "C":
        neighbours.append(-1)

    # If right neighbour is a cysteine aminoacid, add corresponding direction to the neighbours list       
    if protein.grid[pos_y][pos_x + 1] == "C":
        neighbours.append(1)

    return neighbours


def calculate_score(protein):
    """
    This function calculates the score of a protein route. It loops over the protein length to first find HH- and HC-bonds.
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
                        protein.bonds.append(x)
                
                # If the final aminoacid is not connected to the hydrophobic neighbour, add 0.5 to the score
                elif x >= len(protein.route) - 1:
                    if protein.route[x-1] != - i:
                        protein.score += 0.5      
                        protein.bonds.append(x)            

                # If any aminoacid (but the first and final) is not connected to the hydrophobic neighbour, add 0.5 to the score 
                elif protein.route[x] != i and protein.route[x-1] != -i:
                    protein.score += 0.5
                    protein.bonds.append(x)

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
            # Nneighbouring hydrophobic aminoacids are not checked as HC/CH-bonds are already found in the previous part (line 173 - 191)
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
            
            # If the move is horizontal, update the x-coordinate 
            if abs(protein.route[x]) <= 1:
                pos_x += protein.route[x]
            # If the move is vertical, update the y-coordinate
            else:
                pos_y -= int(protein.route[x]/2)



def find_optimum(times, protein_structure):
    """
    This function finds the optimal protein route, number of bonds and score, given a number of times to run the algorithm 
    and a protein structure as specified by the user. It first opens a result file meant to write the results in. 
    Subsequently, the function creates random protein routes and their scores x times (indicated by the user), and writes
    these outcomes in the result file. The function remembers the highest score and prints the corresponding route, score and number of bonds.
    The route that corresponds to the highest score is visualised in a scatterplot.  
    """

    # Open the resultfile
    f = open("resultfile", "w")

    # Initialize rows to order Algoritm, Route and Stability values.
    with f:
        writer = csv.writer(f)
        writer.writerow(["Algoritme", "Route", "Stability"])

    # Set the highest score to zero and the best protein and grid to none
    highest_score = 0
    best_proteins = []
    all_proteins = []
    best_grid = None

    # Repeat the following process as many times as indicated by the user
    for i in range(times):
        
        # Call the protein structure
        protein = Protein(protein_structure)

        # Call the function to create a random protein 
        create_random(protein = protein)

        # If the latter function makes a valid route, calculate its score and write the route and its score in the result file 
        if protein.wrong_protein == False:
            calculate_score(protein = protein)
            f = open("resultfile", "a")
            with f:
                writer = csv.writer(f)
                writer.writerow(["Random", protein.route, protein.score])

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

    print("beste protein:", best_protein.route)
    print("beste score:", highest_score)
    print("bonds:", best_protein.bonds)
    print(best_protein.bonds)

    # Initialize axis
    plt.axis([-10,10,-10,10])

    # Initialize staring position
    pos_x = best_protein.start_x
    pos_y = best_protein.start_y

    # Check what type of amino acid the first amino acid of the protein is and plot this one
    if best_protein.name_list[0] == "H":
        plt.scatter(pos_x, pos_y, c= "blue")
        plt.annotate("0", (pos_x, pos_y))

    elif best_protein.name_list[0] == "C":
        plt.scatter(pos_x, pos_y, c= "green")
        plt.annotate("0", (pos_x, pos_y))

    else:
        plt.scatter(pos_x, pos_y, c= "red")
        plt.annotate("0", (pos_x, pos_y))

    # Iterates over the whole protein
    for x in range(len(best_protein.route)):

        # Checks if the direction is upwards
        if best_protein.route[x] == 2:

            # If so, changes y coordinate
            pos_y += 1

            # Check what type of amino acid the first amino acid of the protein is and plot this one
            if best_protein.name_list[x+1] == "H":
                plt.scatter(pos_x, pos_y, c= "blue")
                plt.annotate(x+1, (pos_x, pos_y))

            elif best_protein.name_list[x+1] == "C":
                plt.scatter(pos_x, pos_y, c= "green")
                plt.annotate(x+1, (pos_x, pos_y))

            else:
                plt.scatter(pos_x, pos_y, c= "red")
                plt.annotate(x+1, (pos_x, pos_y))

            # Plot the connection between the elements
            plt.plot([pos_x, pos_x],[pos_y -1, pos_y], c="grey")
            
        # Checks if the direction is upwards
        if best_protein.route[x] == -2:

            # If so, change y-coordinate
            pos_y -= 1

            # Check what type of amino acid the first amino acid of the protein is and plot this one
            if best_protein.name_list[x+1] == "H":
                plt.scatter(pos_x, pos_y, c= "blue")
                plt.annotate(x+1, (pos_x, pos_y))

            elif best_protein.name_list[x+1] == "C":
                plt.scatter(pos_x, pos_y, c= "green")
                plt.annotate(x+1, (pos_x, pos_y))
                
            else:
                plt.scatter(pos_x, pos_y, c= "red")
                plt.annotate(x+1, (pos_x, pos_y))

            # Plot the connection between the elements
            plt.plot([pos_x, pos_x],[pos_y +1, pos_y], color="grey")

        # Checks if the direction is to the right
        if best_protein.route[x] == 1:

            # If so, change x-coordinate
            pos_x += 1

            # Check what type of amino acid the first amino acid of the protein is and plot this one
            if best_protein.name_list[x+1] == "H":
                plt.scatter(pos_x, pos_y, c= "blue")
                plt.annotate(x+1, (pos_x, pos_y))

            elif best_protein.name_list[x+1] == "C":
                plt.scatter(pos_x, pos_y, c= "green")
                plt.annotate(x+1, (pos_x, pos_y))

            else:
                plt.scatter(pos_x, pos_y, c= "red")
                plt.annotate(x+1, (pos_x, pos_y))

            # Plot the connection between the elements
            plt.plot([pos_x - 1, pos_x],[pos_y , pos_y], c="grey")

        # Checks if the direction is to the left
        if best_protein.route[x] == -1:

            # If so, change x-coordinate
            pos_x -= 1

            # Check what type of amino acid the first amino acid of the protein is and plot this one
            if best_protein.name_list[x+1] == "H":
                plt.scatter(pos_x, pos_y, c= "blue")
                plt.annotate(x+1, (pos_x, pos_y))

            elif best_protein.name_list[x+1] == "C":
                plt.scatter(pos_x, pos_y, c= "green")
                plt.annotate(x+1, (pos_x, pos_y))

            else:
                plt.scatter(pos_x, pos_y, c= "red")
                plt.annotate(x+1, (pos_x, pos_y))

            # Plot connection between the elements
            plt.plot([pos_x + 1, pos_x],[pos_y , pos_y], color="grey")
    
    # Set interval between axis points to one
    plt.xticks(np.arange(-10, 10, 1.0))
    plt.yticks(np.arange(-10, 10, 1.0))

    # Plot a grid
    plt.grid()

    # Show plot
    plt.show()

    


if __name__ == "__main__":
    """
    This function first asks the user to give a number of times to repeat the algoritm and the protein structure.
    Subsequently, it calls the find_optimum function giving along the times and protein structure as provided by the user. 
    """
  
    # Asks user for the number of times to run algorithm 
    times = input("Please indicate the number of times to repeat the algoritm: ")

    # Converts number of times to an integer
    times_number = int(times)

    # Asks user for the protein structure
    protein_structure = input("Please give the protein structure: ")

    # Calls the find_optimum function giving along the times and protein structure as provided by the user
    find_optimum(times = times_number, protein_structure = protein_structure)

