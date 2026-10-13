"""Recorta amostras do GSE286313 (EPIC v2, sangue total) para data/amostras/.

Uso:
    python scripts/preparar_amostra_geo.py CAMINHO/GSE286313_MatrixProcessed_GPL33022.csv.gz \
                                           CAMINHO/GSE286313-GPL33022_series_matrix.txt.gz

Arquivos do GEO (baixar antes):
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286313/suppl/GSE286313_MatrixProcessed_GPL33022.csv.gz  (622 MB)
    https://ftp.ncbi.nlm.nih.gov/geo/series/GSE286nnn/GSE286313/matrix/GSE286313-GPL33022_series_matrix.txt.gz

O que faz:
    - mantém TODOS os probes e só as amostras de AMOSTRAS (escolhidas pela idade,
      antes de rodar o relógio, para cobrir 21 a 86 anos);
    - controle de qualidade: beta com detection p-value > 0,01 vira NaN;
    - arredonda os betas para 4 casas (reduz o arquivo; precisão muito acima do ruído do array);
    - salva betas (probes x amostras) e metadados (idade, sexo) com o GSM como nome da coluna.
"""
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SAIDA = ROOT / "data" / "amostras"
LIMIAR_DETP = 0.01

# Títulos no GEO, ordenados por idade (21 a 86 anos)
AMOSTRAS = [
    "CALERIE_5_EPICv2", "CALERIE_12_EPICv2", "CALERIE_7_EPICv2", "CALERIE_13_EPICv2",
    "CLHNS_3_EPICv2", "CALERIE_10_EPICv2", "CLHNS_5_EPICv2", "CLHNS_2_EPICv2",
    "CLHNS_13_EPICv2", "VHAS_4_EPICv1_EPICv2", "VHAS_19_EPICv1_EPICv2", "VHAS_7_EPICv1_EPICv2",
]


def ler_metadados(series_matrix):
    # Extrai título, GSM, sexo e idade das linhas !Sample_* da series matrix
    campos = {}
    with gzip.open(series_matrix, "rt", encoding="utf-8") as f:
        for linha in f:
            if linha.startswith("!Sample_"):
                partes = next(csv.reader([linha.rstrip("\n")], delimiter="\t"))
                campos.setdefault(partes[0], []).append(partes[1:])
    meta = pd.DataFrame({"gsm": campos["!Sample_geo_accession"][0], "titulo": campos["!Sample_title"][0]})
    for carac in campos["!Sample_characteristics_ch1"]:
        chave = carac[0].split(":")[0].strip().lower()
        meta[chave] = [c.split(":", 1)[1].strip() for c in carac]
    meta = meta.rename(columns={"sex": "sexo", "age": "idade", "tissue": "tecido"})
    meta["idade"] = meta["idade"].astype(float).round(1)
    return meta.set_index("titulo")


def main(matriz, series_matrix):
    meta = ler_metadados(series_matrix).loc[AMOSTRAS]
    colunas = ["ID_REF"] + [c for a in AMOSTRAS for c in (a, f"{a} Detection Pval")]
    df = pd.read_csv(matriz, usecols=colunas, index_col="ID_REF")

    betas = df[AMOSTRAS].copy()
    detp = df[[f"{a} Detection Pval" for a in AMOSTRAS]].to_numpy()
    mascarados = detp > LIMIAR_DETP
    betas = betas.mask(mascarados).round(4)
    betas.columns = meta["gsm"].to_numpy()
    betas.index.name = "probe"

    SAIDA.mkdir(parents=True, exist_ok=True)
    betas.to_csv(SAIDA / "GSE286313_EPICv2_betas.csv.gz", compression="gzip")
    meta.reset_index().set_index("gsm")[["titulo", "tecido", "sexo", "idade"]].to_csv(
        SAIDA / "GSE286313_EPICv2_metadados.csv")

    print(f"{betas.shape[0]:,} probes x {betas.shape[1]} amostras")
    print(f"Mascarados por detection p > {LIMIAR_DETP}: {mascarados.sum():,} valores "
          f"({mascarados.mean():.3%}); NaN total: {betas.isna().sum().sum():,}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
