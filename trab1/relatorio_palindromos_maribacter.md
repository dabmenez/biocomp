# Relatório de Análise de Palíndromos - Maribacter sp. HTCC2170

**Genoma:** CP002157.1 (3,868,304 bp)

Este relatório apresenta os resultados da análise de palíndromos maximais em duas regiões específicas do genoma da bactéria Maribacter sp. HTCC2170.

---

## Região 1: 82583..83599

Esta região possui 1017 pares de bases.

### Genes Encontrados
A região contém os seguintes genes:

- **FB2170_16476**: aspartate aminotransferase (posições 81393..82586)
- **FB2170_16481**: UDP-N-acetylenolpyruvoylglucosamine reductase (posições 82583..83599)

### Palíndromos Encontrados

**Palíndromos de 4 bases:** 15 sequências diferentes
- TTAA (posições: [82607, 82671, 82745, 83129, 83298, 83312, 83427])
- GCGC (posições: [82676, 83166])
- AATT (posições: [82695, 82739, 82985, 83054, 83181, 83274, 83319, 83533, 83563])
- ... e mais 12 sequências

**Palíndromos de 6 bases:** 12 sequências diferentes
- AAGCTT (posições: [82635])
- AAATTT (posições: [82715, 82851])
- TAATTA (posições: [82746])
- ... e mais 9 sequências

**Palíndromos de 8 bases:** 2 sequências diferentes
- ATTGCAAT (posições: [82690, 83049])
- TTTATAAA (posições: [83556])

**Palíndromos de 10 bases:** 2 sequências diferentes
- AAAATATTTT (posições: [82645])
- CTTTTAAAAG (posições: [82829])

### Resumo Geral

Distribuição dos palíndromos por tamanho:
- 10 bases: 2 sequências
- 8 bases: 3 sequências
- 6 bases: 13 sequências
- 4 bases: 43 sequências
- 2 bases: 184 sequências

**Maior palíndromo encontrado:**
- Sequência: AAAATATTTT
- Tamanho: 10 bases
- Localização: 82645..82654

### Sítios de Enzimas de Restrição

Alguns palíndromos correspondem a sítios de enzimas de restrição conhecidas:

- GAATTC → EcoRI (de Escherichia coli R)
- AAGCTT → HindIII (de Haemophilus influenzae Rd)

---

## Região 2: 297449..299453

Esta região possui 2005 pares de bases.

### Genes Encontrados
A região contém os seguintes genes:

- **FB2170_17446**: glutamate synthase (ferredoxin) (posições 294775..299283)
- **FB2170_17451**: NADPH-dependent glutamate synthase, small subunit (posições 299292..300758)

### Palíndromos Encontrados

**Palíndromos de 4 bases:** 16 sequências diferentes
- CATG (posições: [297465, 297961, 298290, 298462, 298704, 298936, 299291])
- GCGC (posições: [297557, 299000, 299172])
- AATT (posições: [297561, 297839, 298277, 298413, 298484, 298491, 298643, 298648, 298760, 298897, 298971, 299117, 299126, 299146, 299185, 299267])
- ... e mais 13 sequências

**Palíndromos de 6 bases:** 18 sequências diferentes
- AACGTT (posições: [297518, 298880])
- TTCGAA (posições: [297597])
- GTTAAC (posições: [297945, 298249])
- ... e mais 15 sequências

**Palíndromos de 8 bases:** 7 sequências diferentes
- TTGCGCAA (posições: [297701])
- GGAATTCC (posições: [297713])
- TGATATCA (posições: [297737])
- ... e mais 4 sequências

### Resumo Geral

Distribuição dos palíndromos por tamanho:
- 8 bases: 7 sequências
- 6 bases: 21 sequências
- 4 bases: 81 sequências
- 2 bases: 332 sequências

**Maior palíndromo encontrado:**
- Sequência: TTGCGCAA
- Tamanho: 8 bases
- Localização: 297701..297708

### Sítios de Enzimas de Restrição

Alguns palíndromos correspondem a sítios de enzimas de restrição conhecidas:

- GAATTC → EcoRI (de Escherichia coli R)
- GCTAGC → NheI (de Neisseria mucosa)
- GCCGGC → NaeI (de Nocardia aerocolonigenes)

---

## Resultados Principais

### Maior Palíndromo Identificado

O maior palíndromo encontrado em ambas as regiões é:
- Sequência: AAAATATTTT
- Tamanho: 10 bases

Este palíndromo não corresponde a nenhuma enzima de restrição conhecida.

### Enzimas de Restrição Encontradas

Total de 4 enzimas diferentes foram identificadas:

- EcoRI (sítio: GAATTC) - Escherichia coli R
- HindIII (sítio: AAGCTT) - Haemophilus influenzae Rd
- NaeI (sítio: GCCGGC) - Nocardia aerocolonigenes
- NheI (sítio: GCTAGC) - Neisseria mucosa

## Respostas às Perguntas Solicitadas

### 1. Palíndromos maximais de tamanho k

O programa foi testado com diferentes valores de k (4, 6, 8, 10, 12, 14, 16, 18, 20).

**Resultados para k=6:**
- Região 1 (82583..83599): 12 sequências diferentes
  - AAGCTT (posições: [82635])
  - AAATTT (posições: [82715, 82851])
  - TAATTA (posições: [82746])
  - ATTAAT (posições: [82781])
  - GATATC (posições: [82799])
  - ACATGT (posições: [82908])
  - AATATT (posições: [82952])
  - CAGCTG (posições: [82993])
  - ATATAT (posições: [83101])
  - GAATTC (posições: [83324])
  - AACGTT (posições: [83436])
  - CCATGG (posições: [83488])
- Região 2 (297449..299453): 18 sequências diferentes
  - AACGTT (posições: [297518, 298880])
  - TTCGAA (posições: [297597])
  - GTTAAC (posições: [297945, 298249])
  - TTGCAA (posições: [298118])
  - TTTAAA (posições: [298170])
  - ACATGT (posições: [298177])
  - ATTAAT (posições: [298315])
  - AAATTT (posições: [298532])
  - TAATTA (posições: [298583])
  - ACCGGT (posições: [298654])
  - GCCGGC (posições: [298663])
  - CGCGCG (posições: [298691])
  - TGATCA (posições: [298932])
  - TTATAA (posições: [299031, 299369])
  - GCTAGC (posições: [299163])
  - TATATA (posições: [299277])
  - GAATTC (posições: [299316])
  - ACGCGT (posições: [299357])

### 2. Maior palíndromo maximal encontrado

O maior palíndromo maximal encontrado é: AAAATATTTT
- Tamanho: 10 bases

Teste de maximalidade (verificando até que tamanho existem palíndromos):
- 2 bases: Encontrados
- 4 bases: Encontrados
- 6 bases: Encontrados
- 8 bases: Encontrados
- 10 bases: Encontrados
- 12 bases: Nenhum encontrado

### 3. Enzimas de restrição encontradas

Foram encontradas 4 enzimas de restrição diferentes:

1. EcoRI (sítio: GAATTC) - Escherichia coli R
2. HindIII (sítio: AAGCTT) - Haemophilus influenzae Rd
3. NaeI (sítio: GCCGGC) - Nocardia aerocolonigenes
4. NheI (sítio: GCTAGC) - Neisseria mucosa

### 4. Análise de CDS/ORF

Ambos os trechos analisados contêm genes anotados:

**Região 1 (82583..83599):**
- FB2170_16476: aspartate aminotransferase
- FB2170_16481: UDP-N-acetylenolpyruvoylglucosamine reductase
**Região 2 (297449..299453):**
- FB2170_17446: glutamate synthase (ferredoxin)
- FB2170_17451: NADPH-dependent glutamate synthase, small subunit

---
