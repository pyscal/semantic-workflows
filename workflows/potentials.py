import pandas as pd
import os

filename_Fe_Ack = os.path.join(os.path.dirname(__file__), "potentials/Fe_Ack.eam")

EAM_Ackland_Fe = pd.DataFrame(
    {
        "Config": [["pair_style eam/alloy\n", "pair_coeff * * Fe_Ack.eam Fe\n"]],
        "Filename": [[filename_Fe_Ack]],
        "Model": ["EAM"],
        "Name": ["EAMAckNin"],
        "Species": [["Fe"]],
    }
)
