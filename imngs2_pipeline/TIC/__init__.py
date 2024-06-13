import sys
from pathlib import Path, PurePath

sys.path.insert(0, str(Path(PurePath(__file__)).parent.parent))  # Add the parent directory to the import path