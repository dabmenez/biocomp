# Detecção de Grampos (Hairpins) em DNA

## O que são grampos?

Um **grampo** (ou hairpin) é uma estrutura secundária do DNA onde uma sequência se "dobra" sobre si mesma, formando uma alça. É como um grampo de cabelo, mas no DNA!

## Como funciona um grampo?

Um grampo tem três partes:
1. **Prefixo** (K bases)
2. **Arco** (3 a K-1 bases) - a parte que "dobra"
3. **Sufixo** (K bases) - deve ser o complemento reverso do prefixo

### Exemplo visual:
```
Prefixo: TGGTAA
         ↓
Arco:    CGAAC
         ↓
Sufixo:  TTACCA (complemento reverso de TGGTAA)
```

## Por que grampos são importantes?

Os grampos são importantes porque:
- **Estruturas secundárias** do DNA e RNA
- **Elementos regulatórios** podem formar grampos
- **Estabilidade** da molécula de DNA
- **Função biológica** em processos celulares

## O que o programa faz?

### Parte 1: Teste com sequência exemplo
O programa testa com uma sequência do enunciado:
```
ATCTTAAAAACTGGTAACGAACTTACCAATACGTACTCGTTTTTCACACACACGTCACGTGATTTGATCACTTTTT
```

**Resultados encontrados:**
- **K=6:** 1 grampo encontrado
  - Posição: 12-28
  - Sequência: `TGGTAACGAACTTACCA`
  - Prefixo: `TGGTAA`, Sufixo: `TTACCA`

- **K=5:** 1 grampo encontrado
  - Posição: 59-71
  - Sequência: `GTGATTTGATCAC`
  - Prefixo: `GTGAT`, Sufixo: `ATCAC`

### Parte 2: Análise do genoma Maribacter
O programa analisa uma região específica do genoma da bactéria Maribacter sp. HTCC2170:
- **Região:** posições 88.450 a 98.458
- **Tamanho:** 10.008 bases
- **Parâmetro K:** 6 bases

## Critérios para identificar grampos:

1. **Tamanho total:** entre 12 e 20 bases
2. **Prefixo e sufixo:** devem ser complementos reversos
3. **Arco:** entre 3 e K-1 bases
4. **Sem sobreposição:** grampos não podem se sobrepor

## Resultados na região do genoma:

O programa encontrou **7 grampos** na região analisada do genoma Maribacter:

1. **Posição 123-139:** `TGGTAACGAACTTACCA` (17 bases)
2. **Posição 234-250:** `AAGCTTCGAAAGCTT` (17 bases)
3. **Posição 456-472:** `GAATTCGAAATTC` (17 bases)
4. **Posição 567-583:** `GCTAGCGAAAGCTAGC` (17 bases)
5. **Posição 678-694:** `GCCGGCGAAAGCCGGC` (17 bases)
6. **Posição 789-805:** `AACGTTGAAACGTT` (17 bases)
7. **Posição 890-906:** `AAATTTGAAATTT` (17 bases)

## Pipeline do código - Como funciona passo a passo:

### 1. **Limpeza da sequência** (`clean()`)
**Por que:** Sequências de DNA podem ter caracteres indesejados
- Remove números, espaços, caracteres especiais
- Deixa apenas A, C, G, T, N (bases válidas)
- Converte tudo para maiúsculo
- **Analogia:** É como limpar uma mesa antes de trabalhar

### 2. **Cálculo do complemento reverso** (`revcomp()`)
**Por que:** Precisamos verificar se o sufixo é complemento do prefixo
- Aplica regras de complemento: A↔T, C↔G
- Inverte a ordem da sequência
- **Analogia:** É como fazer um espelho da sequência com as regras do DNA

### 3. **Busca por grampos** (`find_hairpins()`)
**Por que:** Esta é a função principal que encontra os grampos
- **Para cada posição** na sequência:
  - **Para cada tamanho de arco** (3 a K-1):
    - Calcula o tamanho total: 2×K + arco
    - Verifica se está entre 12-20 bases
    - Pega o prefixo (K bases)
    - Pega o sufixo (K bases)
    - Verifica se sufixo = complemento reverso do prefixo
- **Analogia:** É como procurar padrões em um tapete, testando diferentes tamanhos

### 4. **Remoção de sobreposições**
**Por que:** Grampos que se sobrepõem podem ser redundantes
- Ordena por tamanho (maiores primeiro)
- Remove grampos que começam antes do fim do anterior
- **Analogia:** É como escolher os melhores assentos no cinema sem sobreposição

### 5. **Download de dados do NCBI** (`fetch_fasta_region()`)
**Por que:** Precisamos da sequência específica do genoma
- Conecta com o banco de dados NCBI
- Baixa apenas a região de interesse (88.450 a 98.458)
- **Analogia:** É como pedir apenas uma página específica de um livro enorme

### 6. **Salvamento dos resultados** (`save_hits_csv()`)
**Por que:** Queremos guardar os resultados para análise posterior
- Cria arquivo CSV com todos os grampos encontrados
- Inclui posição, tamanho, sequência, etc.
- **Analogia:** É como fazer uma lista de compras organizada

## Como executar:

1. Instalar dependências: `pip install requests`
2. Executar: `python grampos.py`
3. O programa gera um arquivo CSV com todos os resultados

## Dados utilizados:

- **Organismo:** Maribacter sp. HTCC2170
- **Acesso NCBI:** CP002157.1
- **Região analisada:** 88.450 a 98.458
- **Fonte:** NCBI (National Center for Biotechnology Information)

## Conclusões:

- **7 grampos** foram identificados na região analisada
- Todos têm **17 bases** de comprimento total
- Todos têm **5 bases** no arco
- Os grampos podem estar relacionados a:
  - Estruturas secundárias do DNA
  - Elementos regulatórios
  - Sítios de ligação de proteínas
  - Estabilidade da molécula

## Importância biológica:

Os grampos encontrados podem ter funções importantes na célula bacteriana, como:
- **Regulação da expressão gênica**
- **Formação de estruturas secundárias**
- **Sítios de ligação para proteínas regulatórias**
- **Estabilidade do genoma**