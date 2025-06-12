import pandas as pd

potential = pd.DataFrame(
    {
        "Config": [["pair_style eam/alloy\n", "pair_coeff * * Fe_Ack.eam Fe\n"]],
        "Filename": [["/u/smenon/Fe-workflows/Fe_Ack.eam"]],
        "Model": ["EAM"],
        "Name": ["EAMAckNin"],
        "Species": [["Fe"]],
    }
)
