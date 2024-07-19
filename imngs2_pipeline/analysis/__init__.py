import sys
import os

# Get the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the parent directory
parent_dir = os.path.dirname(current_dir)

# Insert the parent directory at the beginning of sys.path
sys.path.insert(0, parent_dir)
