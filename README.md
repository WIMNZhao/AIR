# AI-augmented R-group exploration in medicinal chemistry

This repo accompanies the preprint ["AI-augmented R-group exploration in medicinal chemistry"](https://www.biorxiv.org/content/10.1101/2024.09.23.614417v1).

### Setting up the environment 

Create a new conda environment and install rdkit:
`conda create -c conda-forge -n <your-env-name> rdkit`

Activate your environment and install the rest of the requirements:
`conda activate <your-env-name>`
`pip install -r requirements.txt`

### How to run

`cd data/drd2`

Train the model

`../../r_fit.py`

Predict the activity of test compounds 

`../../r_predict.py`

Peform the Free-Wilson analysis, preferably by training a new model on the entire dataset ref.smi+test.smi

`../../r_FW.py`

Get help 

`../../r_fit.py -h`

### Notes
The core is defined in SMARTS strings. The kekular structure of benzene C1=C-C=C-C=C1 is not valid.


