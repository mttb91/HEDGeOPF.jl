
using HEDGeOPF
# Import chosen LP and NLP solvers
import HiGHS, Ipopt

# Specify path to YAML configuration file
path = @__DIR__
# Generate dataset of AC-OPF instances
generate_dataset(path; filename = "settings.yaml")
