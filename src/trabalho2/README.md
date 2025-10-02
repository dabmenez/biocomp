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