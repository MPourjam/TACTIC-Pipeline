# End-to-end tests of the containerized pipeline

Install the test dependencies in end-to-end/requirements.txt and run the tests with pytest
**from the root of the repository**.
```bash
pip install -r test/end-to-end/requirements.txt
pytest
```

The test have a higher runtime, as the setup builds the latest container, and runs the 
pipeline on some test data. (21 minutes on my machine for the current test cases)
The test data is truncated to reduce the runtime of the tests and keep the file size small.

Currently, the results are just checked for the presence of the output files.
The tests cover multiple folder structures (flat, all files in one subfolder). This might 
be extended to more complex folder structures (multiple subfolders, nested subfolders, etc.).
