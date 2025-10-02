# README — Detecção de Grampos (Hairpins) em DNA — explicação passo a passo

Este é o meu guia didático para explicar o script `grampos.py` que identifica grampos (hairpins) em sequências de DNA.
Escrevi em primeira pessoa para eu conseguir apresentar com segurança.

---

## 1) O que eu preciso encontrar

Um grampo é uma subcadeia com este formato:

prefixo (K bases) + arco (3..K-1 bases) + sufixo (K bases)

A regra principal é: o sufixo deve ser o reverse-complement do prefixo.

- complemento: A ↔ T, C ↔ G (N fica N)
- reverse: inverter a ordem
Exemplo: revcomp("TGGTAA") = "TTACCA".

Além disso:
- tamanho total do grampo entre 12 e 20;
- K > 4;
- arco entre 3 e K − 1;
- não reportar grampos que se sobrepõem.

---

## 2) O que o programa faz (Resumo)

1. Parte (1): roda na sequência do enunciado com K=6 e K=5 e imprime os grampos.
2. Parte (2): baixa do NCBI o intervalo 88450..98458 do acesso CP002157.1 (Maribacter) e procura grampos (padrão K=6). Imprime e salva um CSV.

As posições são 1-based (como é comum em Biologia).

---

## 3) Como executar

Requisitos:
- Python 3.8+
- requests

Instalação e execução:
```bash
python -m venv .venv
source .venv/bin/activate         # Linux/macOS
# .venv\Scripts\activate        # Windows

pip install requests
python grampos.py
```

---

## 4) Exemplo concreto da Parte (1) — passo a passo

Sequência do enunciado (abreviada visualmente):
S = ATCTTAAAAACTGGTAACGAACTTACCAATACGTACTCGTTTTTCACACACACGTCACGTGATTTGATCACTTTTT

### 4.1) Caso K = 6

Saída do programa:
pos 12-28  len=17  loop=5  hairpin='TGGTAACGAACTTACCA'
prefix = TGGTAA
suffix = TTACCA

Como chego nisso:
- início 1-based = 12 → índice 0-based i = 11
- prefix = S[11:17] = "TGGTAA"
- tentei loop = 3, 4, 5; funcionou com loop = 5
- L = 2*K + loop = 12 + 5 = 17 → j = i + L = 28 (exclusivo)
- suffix = S[11 + 6 + 5 : 28] = S[22:28] = "TTACCA"
- revcomp("TGGTAA") = "TTACCA" → válido e L está entre 12 e 20

### 4.2) Caso K = 5

Saída:
pos 59-71  len=13  loop=3  hairpin='GTGATTTGATCAC'
prefix = GTGAT
suffix = ATCAC

Cálculo:
- i = 58 (0-based) → start = 59 (1-based)
- prefix = S[58:63] = "GTGAT"
- loop = 3
- L = 2*5 + 3 = 13 → j = 71
- suffix = S[66:71] = "ATCAC"
- revcomp("GTGAT") = "ATCAC" → ok

Observação: o sufixo é mostrado na orientação original da sequência (eu comparo com o reverse-complement do prefixo por trás).

---

## 5) Como o algoritmo decide os grampos (função principal)

def find_hairpins(seq, K, min_total=12, max_total=20):
    S = clean(seq)  # A/C/G/T/N maiúsculo
    n = len(S)
    hits = []

    for i in range(n):               # posição inicial
        for loop in range(3, K):     # arco 3..K-1
            L = 2*K + loop
            j = i + L
            if L < min_total or L > max_total or j > n:
                continue

            prefix = S[i:i+K]
            suffix = S[i+K+loop:j]

            if "N" in prefix or "N" in suffix:
                continue

            if suffix == revcomp(prefix):
                hits.append({
                    "start": i+1, "end": j, "loop": loop, "length": L,
                    "substring": S[i:j], "prefix": prefix, "suffix": suffix
                })

Depois de coletar todos os candidatos, eu removo sobreposição:

    hits.sort(key=lambda h: (-h["length"], h["start"]))
    chosen, last_end = [], 0
    for h in hits:
        if not chosen or h["start"] > last_end:
            chosen.append(h)
            last_end = h["end"]
    chosen.sort(key=lambda h: h["start"])
    return chosen

- Prioriza maiores (critério razoável quando há conflito).
- Seleção gulosa garante que cada novo achado começa depois do fim do último.
- Ordeno por posição no final para imprimir na ordem natural.

Complexidade: O(n*K).

---

## 6) Parte (2): baixando o intervalo do NCBI

def fetch_fasta_region(accession, start, end, max_retries=3, pause=0.8):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    url = (f"{base}?db=nuccore&id={accession}&rettype=fasta&retmode=text"
           f"&strand=1&seq_start={start}&seq_stop={end}&tool=hairpin_script")
    # faz até 3 tentativas com pequena espera, checa header '>' do FASTA,
    # remove o cabeçalho e limpa a sequência com a função clean()

Detalhes:
- uso `strand=1` (fio direto);
- verifico se veio mesmo FASTA (`>` no começo);
- o conteúdo é limpado para deixar só A/C/G/T/N em maiúsculo;
- tem backoff simples para evitar erro de limite da API.

---

## 7) O que sai na tela e no CSV

Impressão típica:
(1) Enunciado  K=6  total=1
pos 12-   28  len=17  loop=5  hairpin='TGGTAACGAACTTACCA'  prefix=TGGTAA  suffix=TTACCA

CSV (Parte 2):
start,end,length,loop,prefix,suffix,substring
123,139,17,5,TGGTAA,TTACCA,TGGTAACGAACTTACCA

Arquivo salvo (K=6):
maribacter_CP002157.1_88450_98458_K6.csv

---

## 8) Por que K=6 por padrão na Parte (2)?

Porque com o limite 12–20, os K possíveis e seus comprimentos L são:
- K=5 → L=13..14
- K=6 → L=15..17  (equilíbrio bom)
- K=7 → L=17..20
- K=8 → L=19..20
- K≥9 → L≥21 (fora do limite)

Se eu quiser testar todos, posso iterar sobre (5,6,7,8) e gerar um CSV para cada.

---

## 9) Relação com a aula de Programação Dinâmica

Na aula usamos DP para alinhamento geral (com gaps, mismatches).  
Aqui eu resolvi com **varredura simples** porque o enunciado exige **pareamento perfeito** no “pescoço” e tamanhos pequenos.  
Se eu quisesse permitir erros, eu trocaria a comparação direta por um alinhamento DP entre `prefix` e `revcomp(suffix)`.

---

## 10) Checklist para eu apresentar

- Definir hairpin e reverse-complement.
- Mostrar o esquema K | loop | K e o limite 12–20.
- Fazer a conta do exemplo K=6 (pos 12–28) e do K=5 (59–71).
- Explicar 1-based vs 0-based.
- Dizer como removo sobreposição priorizando maiores.
- Mostrar o CSV gerado.

---

## 11) Perguntas comuns

- **Por que ignorar N?** Para evitar falsos positivos no pescoço.
- **Daria para fazer com DP?** Sim, se permitir mismatches/gaps. Aqui não precisa.
- **Por que as posições são 1-based?** Convenção biológica; fica mais natural para leitura.

FIM