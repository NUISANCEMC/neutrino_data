import sys
import numpy as np

#Function to calculate the chi-square from numpy arrays
def calc_chisq(h_data, h_mc, h_invcov):

    #Calculate the difference between data and MC and transform into a column vector
    diff = (h_data - h_mc).reshape(1, -1)

    #Use vector/matrix multiplication to calculate chi-square
    chisq = diff @ h_invcov @ diff.T
    return chisq.item()

#Load file paths from command line arguments
xsec_file = sys.argv[1]
inv_file  = sys.argv[2]

#Load cross-section values for data and MC
print("Opening {}".format(xsec_file))
with open(xsec_file, 'r') as file:
    csv = np.loadtxt(file, skiprows=1, delimiter=',')
    h_data = csv[:,1]
    h_mc   = csv[:,2]

#Load inverted covariance matrix
print("Opening {}".format(inv_file))
with open(inv_file, 'r') as file:
    h_invcov = np.loadtxt(file, delimiter=',')

print("Calculating...")
chisq = calc_chisq(h_data, h_mc, h_invcov)
print("Chi-sq: {:.4f}".format(chisq))
