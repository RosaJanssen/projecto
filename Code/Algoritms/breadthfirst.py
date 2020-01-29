from queue import Queue
import copy
from ..Classes.classes import Protein
from math import floor
import matplotlib.pyplot as plt
import numpy as np


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
    

def calculate_score(protein):
    """
    This function calculates the score of a protein route given a protein structure. It loops over the protein length to first find HH-bonds.
    The corresponding points are added to the score. The addition to the score is halved, in order to prevent duplicate bonds.
    The function does not look for CC- or CH-bonds, since the corresponding protein structures are too long to run. 
    """

    # Define x and y coordinates as starting positions (0,0; defined in Protein class)
    pos_x = protein.start_x
    pos_y = protein.start_y

    # Loop over the protein length
    for x in range(protein.length):

        # Check whether the aminoacid is hydrophobic 
        if protein.name_list[x] == "H":

            # If the aminoacid is hydrophobic, check the neighbouring hydrophobic aminoacids
            neighbours = check_neighbours_h(protein = protein, pos_x = pos_x, pos_y = pos_y)

            # Loop over the list of (directions to) hydrophobic neighbours 
            for i in neighbours:
                
                # If the first aminoacid is not connected to the hydrophobic neighbour, add 0.5 to the score 
                if x == 0:
                    if protein.route[x] != i:
                        protein.score += 1

                # If the final aminoacid is not connected to the hydrophobic neighbour, add 0.5 to the score
                elif x >= len(protein.route) - 1:
                    if protein.route[x-1] != - i:
                        protein.score += 1    
      
                # If any aminoacid (but the first and final) is not connected to the hydrophobic neighbour, add 0.5 to the score 
                elif protein.route[x] != i and protein.route[x-1] != -i:
                    protein.score += 1
                    protein.bonds.append(x)
        
        # As long as x is within the range of the route length, 
        if x < len(protein.route):

            # If the move is horizontal, update the x-coordinate 
            if abs(protein.route[x]) <= 1:
                pos_x += protein.route[x]
            # If the move is vertical, update the y-coordinate
            else:
                pos_y -= int(protein.route[x]/2)


def create_queue(protein_string):
    """
    This function creates a queue that lists all possible unique combinations of directions, i. e. routes. 
    Subsequently, it returns a list that only includes routes of the valid length.
    """

def create_queue(protein_string):
    depth = len(protein_string) - 1 
    # global queue
    queue = Queue()
    queue2 = []
    queue.put("")

    while not queue.empty():
        state = queue.get()
        
        
        if len(state) < depth:
            for i in ['A', 'B', 'C', 'D']:
                child = copy.deepcopy(state)
                child += i
                queue.put(child)
                if len(state) == depth - 1:
                    queue2.append(child)
    
    
    

    return queue2


def make_proteins(queue2, protein_string):
    """
    This function analyzes each individual route from the passed on list of routes (queue2). For each route,
    directions are transformed into positions in a coordinate list, starting in the middle of the grid (0,0). 
    Only "possible" positions can be added to the coordinate list, i.e. positions that are not yet in there. 
    Subsequently, the function calculates the score for only those routes that are "possible", 
    and from these returns the best protein with its attributes. 
    """

    # Set the highest score to zero and the best protein and grid to none
    highest_score = 0
    best_protein = None
    best_grid = None

    # Loop through the list of possible routes
    for i in queue2:

        # Initalize empty list for coordinates of the aminoacid
        coordinates = []    

        # Define protein object
        protein = Protein(protein_string)
        
        # Define x and y coordinates as starting positions (0,0; defined in Protein class) and append to coordinate list
        pos_x = protein.start_x
        pos_y = protein.start_y
        coordinates.append((pos_x, pos_y))

        # Make a list of the concerned possible route (i'th element in queue2) 
        route_list = list(i)

        # Loop through the number directions of the route
        for j in range(protein.length - 1):

            # Define direction concerned
            move = route_list[j]

            # If the move is upwards, adapt y-coordinate and append corresponding direction to protein route
            if move == "A":
                pos_y -= 1
                protein.route.append(2)

            # If the move is downward, adapt y-coordinate and append corresponding direction to protein route
            if move == "C":
                pos_y += 1
                protein.route.append(-2)

            # If the move is to the right, adapt x-coordinate and append corresponding direction to protein route
            if move == "D":
                pos_x -= 1
                protein.route.append(-1)

            # If the move is to the left, adapt x-coordinate and append corresponding direction to protein route
            if move == "B":
                pos_x += 1
                protein.route.append(1)

            # Place succeeding aminoacid at the updated x- and y-coordinate
            protein.grid[pos_y][pos_x] = protein.name_list[j+1]
            # Only add valid positions to the coordinate list (coordinates that do not yet exist)
            if (pos_x, pos_y) not in coordinates:
                coordinates.append((pos_x, pos_y))
        
        # If the coordinate list has the valid length, calculate the score of this route given the protein structure
        if len(coordinates) == protein.length:
            calculate_score(protein = protein)
            
            # If the found score is higher than the so far highest score, redefine the highest score, route, protein and grid 
            if protein.score >= highest_score:
                highest_score = protein.score
                best_route = protein.route
                best_protein = protein
                best_grid = protein.grid

    print("best protein", best_protein.route)
    print("best score", best_protein.score/2)
    print(best_protein.bonds)

    return(best_protein)
    

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


# if __name__ == "__main__":
#     """
#     This function first asks the user for the protein structure.
#     Subsequently, it creates the queue of all possible unique routes given the protein structure. From those, it selects only 
#     valid routes and from these finds and plots the best protein. 
#     """

#     # Asks user for the protein structure
#     protein_string = input("Please give the protein structure: ")

#     # Call function to create queue of possible routes given protein structure 
#     queue3 = create_queue(protein_string)

#     # Call function to find best protein 
#     best_protein = make_proteins(queue3)

#     # Call function to plot best protein
#     plot_best_protein(best_protein)