# Análise de Palíndromos Maximais - Maribacter sp. HTCC2170

Este programa realiza análise computacional de palíndromos maximais no genoma da bactéria Maribacter sp. HTCC2170, respondendo às perguntas específicas da tarefa de biologia computacional.

## Como Funciona

### Conceito de Palíndromo em DNA
Um palíndromo de DNA é uma sequência que é igual ao seu complemento reverso. Por exemplo:
- `TATA` é palíndromo porque seu complemento reverso é `ATAT` (que é igual a `TATA` quando lido de trás para frente)
- `GAATTC` é palíndromo porque seu complemento reverso é `GAATTC`

### Algoritmo Implementado
O programa usa uma abordagem de expansão centrada para encontrar palíndromos maximais:

1. **Palíndromos de comprimento ímpar**: Expande a partir de cada base individual
2. **Palíndromos de comprimento par**: Expande a partir de cada par de bases adjacentes
3. **Verificação de maximalidade**: Testa se o palíndromo pode ser estendido mantendo a propriedade palindrômica

### Estrutura do Código

```python
# Funções principais:
- load_genome_data()           # Carrega dados dos arquivos locais
- is_palindrome()              # Verifica se sequência é palíndromo
- find_maximal_palindromes_of_length_k()  # Encontra palíndromos de tamanho k
- find_all_maximal_palindromes()          # Encontra todos os palíndromos maximais
- check_cds_overlap()          # Verifica sobreposição com genes
- map_to_restriction_enzymes() # Mapeia para enzimas de restrição
- generate_report()            # Gera relatório em Markdown
```

## Como Usar

### Execução Completa (Recomendado)
```bash
python bacter_final.py
```
Executa análise completa das duas regiões especificadas e gera relatório.

### Análise Específica
```bash
# Buscar palíndromos de tamanho específico
python bacter_final.py --k 6 --intervals 82583-83599 297449-299453

# Encontrar apenas o maior palíndromo
python bacter_final.py --find-largest --intervals 82583-83599 297449-299453
```

## Dependências

```bash
pip install biopython
```

## Arquivos Necessários

- `maribacter_HTCC2170.fasta` - Sequência genômica (3.7 MB)
- `maribacter_HTCC2170.gb` - Anotações GenBank (7.6 MB)

## Resultados Obtidos

### 1. Palíndromos Maximais de Tamanho k=6

**Região 1 (82583..83599):** 12 palíndromos únicos
- `AAGCTT` → posição 82635 (HindIII)
- `AAATTT` → posições 82715, 82851
- `GAATTC` → posição 83324 (EcoRI)
- E mais 9 palíndromos...

**Região 2 (297449..299453):** 18 palíndromos únicos
- `AACGTT` → posições 297518, 298880
- `GCTAGC` → posição 299163 (NheI)
- `GCCGGC` → posição 298663 (NaeI)
- E mais 15 palíndromos...

### 2. Maior Palíndromo Maximal Encontrado

- **Sequência:** `AAAATATTTT`
- **Tamanho:** 10 bp
- **Posição:** 82645..82654 (Região 1)
- **Teste de maximalidade:** k=2 a k=10 encontrados, k=12+ não encontrados

### 3. Enzimas de Restrição Identificadas

4 enzimas diferentes encontradas:
1. **EcoRI** - sítio: `GAATTC` (Escherichia coli R)
2. **HindIII** - sítio: `AAGCTT` (Haemophilus influenzae Rd)
3. **NaeI** - sítio: `GCCGGC` (Nocardia aerocolonigenes)
4. **NheI** - sítio: `GCTAGC` (Neisseria mucosa)

### 4. Análise de CDS/ORF

**Ambos os trechos contêm genes anotados:**

**Região 1 (82583..83599):**
- FB2170_16476: aspartate aminotransferase
- FB2170_16481: UDP-N-acetylenolpyruvoylglucosamine reductase

**Região 2 (297449..299453):**
- FB2170_17446: glutamate synthase (ferredoxin)
- FB2170_17451: NADPH-dependent glutamate synthase, small subunit

## Dados Genômicos

- **Organismo:** Maribacter sp. HTCC2170
- **Acesso NCBI:** CP002157.1
- **Tamanho:** 3,868,304 bp
- **Fonte:** NCBI (National Center for Biotechnology Information)

## Saída do Programa

```
Executando análise completa da tarefa...
Carregando dados do genoma Maribacter sp. HTCC2170 dos arquivos locais...
Genoma carregado: CP002157.1, comprimento 3,868,304 bp

Análise de 2 região(ões) do genoma Maribacter sp. HTCC2170
Tamanho do genoma: 3,868,304 bp

Gerando relatório completo...
Relatório salvo em: relatorio_palindromos_maribacter.md

Resumo do relatório:
- Análise completa de palíndromos maximais
- Verificação de CDS/ORF nos trechos
- Identificação de enzimas de restrição
- Respostas a todas as perguntas da tarefa
```

## Relatório Gerado

O programa gera automaticamente `relatorio_palindromos_maribacter.md` com:
- Análise detalhada de cada região
- Palíndromos encontrados por tamanho
- Maior palíndromo maximal identificado
- Enzimas de restrição mapeadas
- Respostas completas às perguntas da tarefa
