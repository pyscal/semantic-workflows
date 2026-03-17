import json

nb_path = "wip_notebooks/06_defect_formation.ipynb"
nb = json.load(open(nb_path))
cells = nb["cells"]

for cell in cells:
    src = "".join(cell["source"])

    if "from workflows.pyiron.kg import KnowledgeDict" in src:
        new_src = src.replace(
            "from workflows.pyiron.kg import KnowledgeDict \n",
            "from conceptual_dictionary import ConceptualDict\n",
        )
        cell["source"] = [new_src]

    if src.strip() == "kg = KnowledgeDict()":
        cell["source"] = ["cd = ConceptualDict()\n"]

    if "kg=kg" in src:
        new_src = src.replace("kg=kg", "cdict=cd")
        cell["source"] = [new_src]

last_cell = cells[-1]
if "".join(last_cell["source"]).strip() == "":
    last_cell["source"] = ["cd.to_yaml('defect_formation.yaml')\n"]

with open(nb_path, "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print("Done. Verifying source cells:")
nb2 = json.load(open(nb_path))
for i, cell in enumerate(nb2["cells"]):
    src = "".join(cell["source"])
    if src.strip():
        print("cell", i + 1, ":", repr(src[:130]))
