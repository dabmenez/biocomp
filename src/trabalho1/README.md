# Análise de Palíndromos em DNA - Maribacter sp. HTCC2170

## O que é este trabalho?

Este projeto analisa **palíndromos** no genoma da bactéria **Maribacter sp. HTCC2170**. Um palíndromo de DNA é uma sequência que é igual ao seu complemento reverso - como uma palavra que se lê igual de trás para frente, mas com as regras do DNA.

## Por que palíndromos são importantes?

Os palíndromos no DNA são muito importantes porque:

1. **Enzimas de restrição** reconhecem sequências palindrômicas para cortar o DNA
2. **Elementos regulatórios** podem ter estruturas palindrômicas
3. **Estruturas secundárias** do DNA podem se formar em regiões palindrômicas

## Como funciona?

### Exemplo simples:
- Sequência: `GAATTC`
- Complemento: `CTTAAG` (A↔T, C↔G)
- Reverso: `GAATTC`
- Como são iguais, `GAATTC` é um palíndromo!

### O que o programa faz:

1. **Carrega o genoma** da Maribacter sp. HTCC2170 (3.868.304 bases)
2. **Analisa duas regiões específicas**:
   - Região 1: posições 82.583 a 83.599
   - Região 2: posições 297.449 a 299.453
3. **Procura palíndromos** de diferentes tamanhos
4. **Identifica enzimas de restrição** que reconhecem essas sequências
5. **Verifica se há genes** nessas regiões

## Resultados encontrados:

### 1. Palíndromos de 6 bases (k=6):

**Região 1:** 12 sequências diferentes encontradas
- `AAGCTT` → enzima HindIII
- `GAATTC` → enzima EcoRI
- `AAATTT` → aparece em 2 posições

**Região 2:** 18 sequências diferentes encontradas
- `GCTAGC` → enzima NheI
- `GCCGGC` → enzima NaeI
- `AACGTT` → aparece em 2 posições

### 2. Maior palíndromo encontrado:
- **Sequência:** `AAAATATTTT`
- **Tamanho:** 10 bases
- **Localização:** posição 82.645 a 82.654

### 3. Enzimas de restrição identificadas:
1. **EcoRI** (Escherichia coli) - reconhece `GAATTC`
2. **HindIII** (Haemophilus influenzae) - reconhece `AAGCTT`
3. **NaeI** (Nocardia aerocolonigenes) - reconhece `GCCGGC`
4. **NheI** (Neisseria mucosa) - reconhece `GCTAGC`

### 4. Genes encontrados nas regiões:

**Região 1:**
- FB2170_16476: aspartato aminotransferase
- FB2170_16481: UDP-N-acetilenolpiruvoilglicosamina redutase

**Região 2:**
- FB2170_17446: glutamato sintase (ferredoxina)
- FB2170_17451: glutamato sintase dependente de NADPH, subunidade pequena

## Como executar:

1. Instalar dependências: `pip install biopython`
2. Executar: `python bacter_final.py`
3. O programa gera um relatório completo em Markdown

## Dados utilizados:

- **Organismo:** Maribacter sp. HTCC2170
- **Acesso NCBI:** CP002157.1
- **Tamanho do genoma:** 3.868.304 pares de bases
- **Fonte:** NCBI (National Center for Biotechnology Information)

## Conclusões:

- Ambas as regiões contêm genes funcionais
- Foram encontrados múltiplos palíndromos de diferentes tamanhos
- 4 enzimas de restrição diferentes foram identificadas
- O maior palíndromo tem 10 bases de comprimento
- Os palíndromos podem estar relacionados a elementos regulatórios ou sítios de corte por enzimas de restrição