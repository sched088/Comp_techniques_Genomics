"""
Problems:
1. Design the HMM representation of protein secondary structure prediction, e.g. three states representing hexlices,
sheets and loops; emissions are the on the 20 amino acid residues.

LEARNING PROBLEM to get probabilities then DECODING PROBLEM!

2. Use maximum likelihood learning to train a HMM with the training data, called HMM_ml.
Note that both amino acid sequence and secondary structure sequence should be used to train HMM_ml.

3. Implement Baum-Welch algorithm for HMM training without using the given secondary structures in the training data.
Learn a HMM with the test data, called HMM_bw.

4. Apply HMM_ml and HMM_bw to predict the secondary structure of the test data.
Report your results by accuracy by the two HMMs. Explain your results.
"""

"""
Key: 
_ = loop
h = helix
e = sheet
"""
"""
1) Load in training data
2) Calculate theta (emission and transition probabilities)
    - akl
        - Need length of pi (state list per aa)
        - Need transition counts
    - ekb
        - Total number of AAs per each state
        - Total number of each AA per state
"""