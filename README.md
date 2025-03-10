# Stage-1
from Bio.Seq import Seq

def translate_dna_to_protein(dna_sequence):
    #"""Translates a DNA sequence into a protein sequence."""
    dna_seq = Seq("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
    protein_seq = dna_seq.translate()
    return str(protein_seq)

print(translate_dna_to_protein(dna_seq))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def logistic_growth(time, K=1000, r=0.1, N0=10, lag_random=5, exp_random=20):
    """
    Simulates a logistic growth curve.
    
    Parameters:
        time (array): Time points.
        K (int): Carrying capacity.
        r (float): Growth rate.
        N0 (int): Initial population size.
        lag_random (int): Randomized lag phase duration.
        exp_random (int): Randomized exponential phase duration.
        
          population (array): Population values at each time point.
    """
    t_lag = np.random.randint(0, lag_random)  # Randomized lag phase
    t_exp = np.random.randint(lag_random, exp_random)  # Randomized exponential phase

    population = []
    for t in time:
        if t < t_lag:  # Lag phase (slow growth)
            N = N0
        elif t < t_exp:  # Exponential phase
            N = N0 * np.exp(r * (t - t_lag))
        else:  # Logistic phase
            N = (K * N0 * np.exp(r * (t - t_lag))) / (K + N0 * (np.exp(r * (t - t_lag)) - 1))
        population.append(N)
    
    return np.array(population)

def generate_growth_dataframe(num_curves=100, time_range=100):
    """
    Generates a DataFrame with multiple logistic growth curves.
    
    Parameters:
        num_curves (int): Number of growth curves to simulate.
        time_range (int): Time range for simulation.
        
    Returns:
        DataFrame: A DataFrame where each column represents a different growth curve.
    """
    time_points = np.linspace(0, time_range, time_range)
    data = {}

    for i in range(num_curves):
        pop_growth = logistic_growth(time_points)
        data[f"Curve_{i+1}"] = pop_growth

    df = pd.DataFrame(data, index=time_points)
    df.index.name = "Time"
    return df

# Example usage
growth_df = generate_growth_dataframe()
print(growth_df.head())

def time_to_80_percent_max(population, time_points, K=1000):
    """
    Determines the time when population reaches 80% of carrying capacity.
    
    Parameters:
        population (array): Population growth data.
        time_points (array): Time points corresponding to population data.
        K (int): Carrying capacity.
    
    Returns:
        float: Time to reach 80% of K.
    """
    threshold = 0.8 * K
    for t, pop in zip(time_points, population):
        if pop >= threshold:
            return t
    return None  # If 80% of K is never reached

# Example usage
t_80 = time_to_80_percent_max(pop_growth, time_points)
print(f"Time to reach 80% of carrying capacity: {t_80}")

def hamming_distance(str1, str2):
    """
    Computes the Hamming distance between two strings.
    
    If the strings are of different lengths, the shorter one is padded.
    
    Parameters:
        str1 (str): First string.
        str2 (str): Second string.
    
    Returns:
        int: Hamming distance.
    """
    max_len = max(len(str1), len(str2))
    str1 = str1.ljust(max_len)  # Pad with spaces if needed
    str2 = str2.ljust(max_len)
    
    return sum(c1 != c2 for c1, c2 in zip(str1, str2))

# Example usage
slack_username = "YourSlackUsername"
twitter_handle = "YourTwitterHandle"
print(f"Hamming Distance: {hamming_distance(slack_username, twitter_handle)}")
