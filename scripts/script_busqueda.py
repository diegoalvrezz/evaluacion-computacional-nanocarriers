import csv
import time
import urllib.parse
import requests

# ==========================
# LISTA DE MOLÉCULAS (60)
# ==========================

MOLECULAS = [

    # -------- MONÓMEROS / DEGRADACIÓN --------
    ("monomeros_degradacion", "L-Lactic acid"),
    ("monomeros_degradacion", "D-Lactic acid"),
    ("monomeros_degradacion", "Glycolic acid"),
    ("monomeros_degradacion", "Succinic acid"),
    ("monomeros_degradacion", "Glutaric acid"),
    ("monomeros_degradacion", "Adipic acid"),
    ("monomeros_degradacion", "6-Hydroxyhexanoic acid"),
    ("monomeros_degradacion", "6-Hydroxycaproic acid"),
    ("monomeros_degradacion", "epsilon-Caprolactone"),
    ("monomeros_degradacion", "Glucosamine"),
    ("monomeros_degradacion", "N-Acetylglucosamine"),
    ("monomeros_degradacion", "Ethylenediamine"),
    ("monomeros_degradacion", "Methyl acrylate"),
    ("monomeros_degradacion", "L-Lactide"),
    ("monomeros_degradacion", "Glycolide"),

    # -------- OLIGÓMEROS / UNIDADES --------
    ("oligomeros_unidades", "Diethylenetriamine"),
    ("oligomeros_unidades", "Triethylenetetramine"),
    ("oligomeros_unidades", "Tetraethylenepentamine"),
    ("oligomeros_unidades", "Pentaethylenehexamine"),
    ("oligomeros_unidades", "Tris(2-aminoethyl)amine"),
    ("oligomeros_unidades", "N,N'-Bis(2-aminoethyl)ethylenediamine"),
    ("oligomeros_unidades", "Methylamine"),
    ("oligomeros_unidades", "Ammonia"),
    ("oligomeros_unidades", "2-Aminoethanol"),
    ("oligomeros_unidades", "2-(2-Aminoethoxy)ethanol"),
    ("oligomeros_unidades", "Bis(2-aminoethyl) ether"),
    ("oligomeros_unidades", "Methyl 3-aminopropionate"),
    ("oligomeros_unidades", "Ethyl acrylate"),
    ("oligomeros_unidades", "Methyl methacrylate"),
    ("oligomeros_unidades", "Acrylic acid"),
    ("oligomeros_unidades", "Lactic acid dimer"),
    ("oligomeros_unidades", "Lactic acid trimer"),
    ("oligomeros_unidades", "Diglycolic acid"),
    ("oligomeros_unidades", "Chitobiose"),
    ("oligomeros_unidades", "Chitotriose"),

    # -------- FÁRMACOS MODELO --------
    ("farmacos_modelo", "Doxorubicin"),
    ("farmacos_modelo", "Paclitaxel"),
    ("farmacos_modelo", "5-Fluorouracil"),
    ("farmacos_modelo", "Methotrexate"),
    ("farmacos_modelo", "Curcumin"),
    ("farmacos_modelo", "Cisplatin"),
    ("farmacos_modelo", "Gemcitabine"),
    ("farmacos_modelo", "Sirolimus"),
    ("farmacos_modelo", "Dexamethasone"),
    ("farmacos_modelo", "Ibuprofen"),

    # -------- LIGANDOS TARGETING --------
    ("ligandos_targeting", "Folic acid"),
    ("ligandos_targeting", "Biotin"),
    ("ligandos_targeting", "D-Mannose"),
    ("ligandos_targeting", "D-Galactose"),
    ("ligandos_targeting", "N-Acetylneuraminic acid"),
    ("ligandos_targeting", "Glucuronic acid"),
    ("ligandos_targeting", "Glucosamine-6-sulfate"),
    ("ligandos_targeting", "Chondroitin sulfate disaccharide"),
    ("ligandos_targeting", "Arginylglycylaspartic acid"),
    ("ligandos_targeting", "Asparagine"),
    ("ligandos_targeting", "Arginine"),
    ("ligandos_targeting", "Aspartic acid"),
    ("ligandos_targeting", "4-Aminobenzoic acid"),
    ("ligandos_targeting", "2,2'-Iminodiethanol"),
    ("ligandos_targeting", "Triethylene glycol"),
]

# ==========================
# FUNCIONES PUBCHEM
# ==========================

def get_cid_by_name(name):
    q = urllib.parse.quote(name)
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{q}/cids/JSON"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return None
    data = r.json()
    return data["IdentifierList"]["CID"][0]


def get_property_txt(cid, prop):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{prop}/TXT"
    r = requests.get(url, timeout=30)
    if r.status_code != 200:
        return None
    value = r.text.strip()
    if value.lower().startswith("status"):
        return None
    return value


# ==========================
# GENERACIÓN CSV
# ==========================

output_file = "base_molecular_pubchem.csv"

with open(output_file, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["grupo", "nombre_entrada", "iupac_name", "smiles_pubchem", "cid", "status"]
    )
    writer.writeheader()

    for grupo, nombre in MOLECULAS:

        row = {
            "grupo": grupo,
            "nombre_entrada": nombre,
            "iupac_name": "",
            "smiles_pubchem": "",
            "cid": "",
            "status": "ok",
        }

        try:
            cid = get_cid_by_name(nombre)

            if cid is None:
                row["status"] = "not_found_name"
                writer.writerow(row)
                continue

            row["cid"] = cid

            iupac = get_property_txt(cid, "IUPACName")
            smiles = get_property_txt(cid, "CanonicalSMILES")

            row["iupac_name"] = iupac if iupac else ""
            row["smiles_pubchem"] = smiles if smiles else ""

            if row["smiles_pubchem"] == "":
                row["status"] = "smiles_missing"

        except Exception as e:
            row["status"] = f"error_{type(e).__name__}"

        writer.writerow(row)
        time.sleep(0.2)

print("CSV generado correctamente:", output_file)