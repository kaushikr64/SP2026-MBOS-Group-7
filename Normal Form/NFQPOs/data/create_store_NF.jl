"""
    Just a simple code used to create and store specific normal forms
"""

using NFCR3BP, JLD2

# System parameters
normalization_order = 16
mu = 0.012150585609624
system_name = "Earth-Moon"
lagrange_point = 1

# Enter Birkhoff for birkhoff (resonant yet unsupported), case sensitive
normalization_type = "RCM"

# Where the data is
storage_folder = "./data/"
storage_file =
    storage_folder *
    system_name *
    "/" *
    normalization_type *
    "/L$lagrange_point/Order$normalization_order" *
    ".jld2"


# Do I want to compute it anew and store it? If so uncomment the below line
normalform_parameters = (
    mu = mu,
    lagrange_point = lagrange_point,
    normalization_order = normalization_order,
    normalization_type = normalization_type,
)
compute_store_NF(storage_file, normalform_parameters)
