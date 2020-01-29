from ..Classes.classes import Protein
import random
from math import floor
import matplotlib.pyplot as plt
import numpy as np


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


def create_random(protein):
    """
    This function creates a random protein route, given the protein structure (protein.name_list) and the start coordinates. 
    It first places the first aminoacid in the middle of a grid. Subsequently, the algorithm checks 
    whether the neighbours of the aminoacid (up, down, left and right) in the grid are free. In addition, it checks whether 
    these free neighbours have H- or C-neighbours. If the latter is true, the algorithm randomly chooses from one or more of
    the neighbours that have H- Or C-neighbours. If this is not the case, the algorithm randomly chooses from the available
    neighbours. The succeeding aminoacid is placed at the chosen spot. The whole latter process is repeated until all
    aminoacids from the protein structure have been placed.
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

        # Initialize empty preference list for preferable neighbouring places to go to
        pref = []

        # If upper neighbour is available (value = 0) add direction to the neighbour list 
        if protein.grid[pos_y - 1][pos_x] == "0":
            neighbours.append(2)
            # If upper neighbour has (one or more) hydrophobic or cysteine neighbour(s), add direction to the preference list 
            if len(check_neighbours_h(protein = protein, pos_x = pos_x, pos_y = pos_y -1)) > 0 or len(check_neighbours_c(protein = protein, pos_x = pos_x, pos_y = pos_y -1)) > 0:
                pref.append(2)
        
        # If lower neighbour is available (value = 0) add direction to the neighbour list 
        if protein.grid[pos_y + 1][pos_x] == "0":
            neighbours.append(-2)
            # If upper neighbour has (one or more) hydrophobic or cysteine neighbour(s), add direction to the preference list 
            if len(check_neighbours_h(protein = protein, pos_x = pos_x, pos_y = pos_y +1)) > 0 or len(check_neighbours_c(protein = protein, pos_x = pos_x, pos_y = pos_y +1)) > 0:
                pref.append(-2)
        
        # If left neighbour is available (value = 0) add direction to the neighbour list 
        if protein.grid[pos_y][pos_x - 1] == "0":
            neighbours.append(-1)
            # If left neighbour has (one or more) hydrophobic or cysteine neighbour(s), add direction to the preference list 
            if len(check_neighbours_h(protein = protein, pos_x = pos_x - 1, pos_y = pos_y)) > 0 or len(check_neighbours_c(protein = protein, pos_x = pos_x -1, pos_y = pos_y)) > 0:
                pref.append(-1)
        
        # If right neighbour is available (value = 0) add direction to the neighbour list 
        if protein.grid[pos_y][pos_x + 1] == "0":
            neighbours.append(1)
            # If right neighbour has (one or more) hydrophobic or cysteine neighbour(s), add direction to the preference list 
            if len(check_neighbours_h(protein = protein, pos_x = pos_x + 1, pos_y = pos_y)) > 0 or len(check_neighbours_c(protein = protein, pos_x = pos_x +1, pos_y = pos_y)) > 0:
                pref.append(1)

        # If there is no neighbouring spot available at all, define the protein object to be a wrong protein 
        if neighbours == []:
            protein.wrong_protein = True
            return

        # If the preference list is empty, choose randomly from the neighbour list
        if pref == []:
            move = random.choice(neighbours)     
        # If there are preferenced neighbours available, choose randomly from those
        else:
            move = random.choice(pref)

        # If the move is horizontal, adapt the x-coordinate 
        if abs(move) <= 1:
            pos_x += move

        # If the move is vertical, adapt the y-coordinate
        else:
            pos_y -= int(move/2)

        # Append move to the protein route
        protein.route.append(move)

        # Place succeeding aminoacid at the updated x- and y-coordinate
        protein.grid[pos_y][pos_x] = protein.name_list[x+1]


def plot_best_protein(best_protein):
    """
    This function makes a visualization of a protein using matplotlib. 
    """
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


def find_optimum(times, protein_structure):
    """
    This function finds the optimal protein route, number of bonds and score, given a number of times to run the algorithm 
    and a protein structure as specified by the user. The function creates random protein routes and their scores x times 
    (indicated by the user), and remembers the highest score. It prints the corresponding route, score and number of bonds.
    The route that corresponds to the highest score is visualised in a scatterplot.  
    """

    # Set the highest score to zero and the best protein and grid to none
    highest_score = 0
    best_protein = None
    best_grid = None
    
    # Repeat the following process as many times as indicated by the user
    for i in range(times):

        # Call the protein structure
        protein = Protein(protein_structure)

        # Call the function to create a random protein 
        create_random(protein = protein)

        # If the latter function makes a valid route, calculate its score  
        if protein.wrong_protein == False:
            calculate_score(protein = protein)

            # If the found score is higher than the so far highest score, redefine the highest score, route, protein and grid 
            if protein.score > highest_score:
                highest_score = protein.score
                best_route = protein.route
                best_protein = protein
                best_grid = protein.grid
    
    # Print best route, score and bonds
    if best_protein != None:
        print("beste protein:", best_protein.route)
        print("beste score:", highest_score)
        print("bonds:", best_protein.bonds)
        plot_best_protein(best_protein)
    else:
        return 0

    
# if __name__ == "__main__":

#     # Asks user for the protein structure
#     protein_string = input("Please give the protein structure: ")
    
#     # Asks user for the number of times to run algorithm 
#     times = input("Please indicate the number of times to repeat the algoritm: ")

#     # Converts number of times to an integer
#     times_number = int(times)

#     # Calls function to find the optimal protein, giving along the times and protein structure as provided by the user
#     while True:
#         if find_optimum(times = times_number, protein_structure = protein_string) != 0:
#             break