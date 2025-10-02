from typing import List, Dict
import re
import csv
import sys
import requests

# Tabela para fazer complemento: A vira T, C vira G, etc.
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    """
    Faz o reverse-complement de uma sequência.
    Exemplo: revcomp("TGGTAA") = "TTACCA"
    """
    return seq.translate(COMP)[::-1]


def clean(seq: str) -> str:
    """
    Remove caracteres que não são DNA e deixa tudo maiúsculo.
    Exemplo: clean("atcg123XYZ") = "ATCG"
    """
    return re.sub(r"[^ACGTNacgtn]", "", seq).upper()


def find_hairpins(seq: str, K: int, min_total: int = 12, max_total: int = 20) -> List[Dict]:
    """
    Procura grampos na sequência.
    Grampo = PREFIXO + LOOP + SUFIXO, onde SUFIXO é o reverse-complement do PREFIXO
    """
    S = clean(seq)
    n = len(S)
    hits: List[Dict] = []

    # Para cada posição na sequência
    for i in range(n):
        # Testa diferentes tamanhos de loop
        for loop in range(3, K):
            L = 2 * K + loop  # Tamanho total
            j = i + L
            
            # Verifica se cabe na sequência e está no tamanho certo
            if L < min_total or L > max_total or j > n:
                continue

            # Pega o prefixo e sufixo
            prefix = S[i:i + K]
            suffix = S[i + K + loop:j]
            
            # Ignora se tem N
            if "N" in prefix or "N" in suffix:
                continue

            # Verifica se é um grampo
            if suffix == revcomp(prefix):
                hits.append({
                    "start": i + 1,
                    "end": j,
                    "loop": loop,
                    "length": L,
                    "substring": S[i:j],
                    "prefix": prefix,
                    "suffix": suffix
                })

    # Remove sobreposições (pega os maiores primeiro)
    hits.sort(key=lambda h: (-h["length"], h["start"]))
    chosen: List[Dict] = []
    last_end = 0
    
    for h in hits:
        if not chosen or h["start"] > last_end:
            chosen.append(h)
            last_end = h["end"]

    chosen.sort(key=lambda h: h["start"])
    return chosen


def fetch_fasta_region(accession: str, start: int, end: int) -> str:
    """
    Baixa uma parte da sequência do NCBI.
    """
    url = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           f"?db=nuccore&id={accession}&rettype=fasta&retmode=text"
           f"&strand=1&seq_start={start}&seq_stop={end}")
    
    resp = requests.get(url, timeout=30)
    txt = resp.text
    lines = txt.splitlines()
    
    # Pula o cabeçalho e pega só a sequência
    seq = clean("\n".join(lines[1:]))
    return seq


def print_hits(hits: List[Dict], label: str) -> None:
    """
    Mostra os grampos encontrados.
    """
    print(f"\n{label}  total={len(hits)}")
    for h in hits:
        print(
            f"pos {h['start']}-{h['end']:>5}  len={h['length']:<2}  loop={h['loop']:<2}  "
            f"hairpin='{h['substring']}'  prefix={h['prefix']}  suffix={h['suffix']}"
        )


def save_hits_csv(hits: List[Dict], path: str) -> None:
    """
    Salva os resultados em um arquivo CSV.
    """
    with open(path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["start", "end", "length", "loop", "prefix", "suffix", "substring"])
        for h in hits:
            w.writerow([h["start"], h["end"], h["length"], h["loop"], h["prefix"], h["suffix"], h["substring"]])


def main() -> int:
    """
    Programa principal.
    """
    
    # Parte 1: sequência do enunciado
    s = "ATCTTAAAAACTGGTAACGAACTTACCAATACGTACTCGTTTTTCACACACACGTCACGTGATTTGATCACTTTTT"
    
    for K in (6, 5):
        hits = find_hairpins(s, K)
        print_hits(hits, f"(1) Enunciado  K={K}")

    # Parte 2: Maribacter
    acc = "CP002157.1"
    a, b = 88450, 98458

    region = fetch_fasta_region(acc, a, b)
    K2 = 6
    hits2 = find_hairpins(region, K2)
    print_hits(hits2, f"(2) Maribacter {acc}:{a}-{b}  K={K2}")

    # Salva em CSV
    out_csv = f"maribacter_{acc}_{a}_{b}_K{K2}.csv"
    save_hits_csv(hits2, out_csv)
    print(f"\nCSV salvo: {out_csv}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
