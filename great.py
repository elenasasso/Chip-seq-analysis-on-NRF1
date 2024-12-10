import re

# Percorso del file
file_path = "20241125-public-4.0.4-4GyVoJ-hg38-all-gene.txt"

# Inizializza i contatori
promoter_count = 0
enhancer_count = 0

# Funzione per verificare la posizione della regione
def classify_region(region):
    try:
        distance = int(region)
        if -1000 <= distance <= 1000:
            return "promoter"
        elif -30000 <= distance <= -1001 or 1001 <= distance <= 30000:
            return "enhancer"
    except ValueError:
        pass
    return None

# Leggi il file e processa le righe
with open(file_path, "r") as file:
    for line in file:
        # Salta le righe di intestazione
        if line.startswith("#") or not line.strip():
            continue

        # Estrai le regioni
        fields = line.split("\t")
        gene_name = fields[0]
        regions = fields[1]

        # Analizza ogni regione associata
        for match in re.findall(r"\(([-+]?\d+)\)", regions):
            classification = classify_region(match)
            if classification == "promoter":
                promoter_count += 1
            elif classification == "enhancer":
                enhancer_count += 1

# Output dei risultati
print(f"Numero di geni nei promotori (-1 kb a +1 kb): {promoter_count}")
print(f"Numero di geni negli enhancer (oltre ±1 kb e entro ±30 kb): {enhancer_count}")
