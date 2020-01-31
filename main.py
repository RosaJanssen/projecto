import csv
from Code.Algoritms import random
from Code.Algoritms import greedy_plotprotein
from Code.Algoritms import greedy_collectdata
from Code.Algoritms import breadthfirst
from Code.Algoritms import lookahead_plotprotein
from Code.Algoritms import lookahead_collectdata


if __name__ == "__main__":

    # ------------------------------------------------------------------ RANDOM ------------------------------------------------------------------
    """
    This function first asks the user to give a number of times to repeat the algoritm and the protein structure.
    Subsequently, it calls the find_optimum function giving along the times and protein structure as provided by the user. 
    """
  
    # Asks user for the number of times to run algorithm 
    #times = input("Please indicate the number of times to repeat the algoritm: ")

    # Converts number of times to an integer
    #times_number = int(times)

    # Asks user for the protein structure
    #given_protein = input("Please give the protein structure: ")

    # Calls the find_optimum function giving along the times and protein structure as provided by the user
    #random.find_optimum(times=times_number, protein_structure=given_protein)

    # ------------------------------------------------------------------ GREEDY ------------------------------------------------------------------

    # ------------------------------------------------------------------ PLOT PROTEIN ------------------------------------------------------------------
    """
    # Asks user for the protein structure
    protein_string = input("Please give the protein structure: ")
    
    # Asks user for the number of times to run algorithm 
    times = input("Please indicate the number of times to repeat the algoritm: ")

    # Converts number of times to an integer
    times_number = int(times)

    # Calls function to find the optimal protein, giving along the times and protein structure as provided by the user
    while True:
        if greedy_plotprotein.find_optimum(times = times_number, protein_structure = protein_string) != 0:
            break
    """
    # ------------------------------------------------------------------ COLLECT DATA ------------------------------------------------------------------

    # Open result file and initialize rows for the algoritm, route and stability 
    #f = open("resultfile_lookahead", "w")
    #with f:
    #    writer = csv.writer(f)
    #    writer.writerow(["Algoritme", "Route", "Stability"])

    # Ask user for the protein structure
    #protein_string = input("Please give the protein structure: ")
    
    # Ask user for the number of times to run algorithm 
    #times = input("Please indicate the number of times to repeat the algoritm: ")

    # Convert number of times to an integer
    #times_number = int(times)

    # Run the algoritm as many times as indicated by the user
    #for a in range(times_number):
    #    while True:
    #        if greedy_collectdata.find_optimum(times = times_number, protein_structure = protein_string) != 0:
    #            best_protein = greedy_collectdata.find_optimum(times = times_number, protein_structure = protein_string)
                # Write results of best proteins in the result file

                
    """

    # ------------------------------------------------------------------ BREADTHFIRST ------------------------------------------------------------------

    #"""
    #This function first asks the user for the protein structure.
    #Subsequently, it creates the queue of all possible unique routes given the protein structure. From those, it selects only 
    #valid routes and from these finds and plots the best protein. 
    #"""

    # Asks user for the protein structure
    # protein_string = input("Please give the protein structure: ")

    # # Call function to create queue of possible routes given protein structure 
    # queue3 = breadthfirst.create_queue(protein_string)

    # # Call function to find best protein 
    # best_protein = breadthfirst.make_proteins(queue3, protein_string=protein_string)

    # # Call function to plot best protein
    # breadthfirst.plot_best_protein(best_protein)

    # ------------------------------------------------------------------ LOOKAHEAD ------------------------------------------------------------------
    # ------------------------------------------------------------------ PLOT PROTEIN ------------------------------------------------------------------
    
    # Ask user for protein structure 
    protein_string = input("Please give the protein structure: ")
    
    # Set look ahead equal to four aminoacids
    lookahead = 4

    # Initializes empty strings to save route
    path = ""
    constant_string = ""

    # Iterate over the amino acids of the protein
    for i in range(len(protein_string)-1):

        # Make substring for the considered amino acids
        substring = protein_string[i:i+lookahead]

        # Create all possible routes
        queue = lookahead_plotprotein.create_queue(substring)

        # Add the previous steps to the route
        for j in range(len(queue)):
            queue[j] = constant_string + queue[j]

        # Find protein with highest stability
        best_protein = lookahead_plotprotein.make_proteins(queue, path, substring)

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
    
    # plot protein
    if best_protein != 0:
        lookahead_plotprotein.plot_best_protein(best_protein)
    

    # ------------------------------------------------------------------ COLLECT DATA ------------------------------------------------------------------

    # Open result file and initialize rows for the algoritm, route and stability 
    #f = open("resultfile_lookahead", "w")
    #with f:
    #    writer = csv.writer(f)
    #    writer.writerow(["Algoritme", "Route", "Stability"])
    
    # Ask user for the protein structure
    #protein_string = input("Please give the protein structure: ")
    #protein_string = input("Please give the protein structure: ")

    # Ask user for the number of times to run algorithm 
    #times = input("Please indicate the number of times to repeat the algoritm: ")

    # Convert number of times to an integer
    #times_number = int(times)

    # Call the generate data function and pass along the number of times to run the algoritm and the protein structure
    #lookahead_collectdata.generate_data(times = times_number, protein_string = protein_string)