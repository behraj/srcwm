import subprocess, sys
subprocess.run([sys.executable, "-m", "pip", "install", "-q",
                "numpy", "pandas", "scipy", "scikit-learn", "h5py",
                "matplotlib", "cellxgene-census"], check=True)
import numpy, pandas, scipy, sklearn, h5py, matplotlib, cellxgene_census, tiledbsoma
print("numpy", numpy.__version__)
print("pandas", pandas.__version__)
print("scipy", scipy.__version__)
print("scikit-learn", sklearn.__version__)
print("h5py", h5py.__version__)
print("cellxgene-census", cellxgene_census.__version__)
print("tiledbsoma", tiledbsoma.__version__)
