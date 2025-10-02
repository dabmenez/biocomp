# Biologia Computacional - Estudos

Repositório com projetos de **Biologia Computacional** - análise de sequências de DNA, palíndromos e grampos.

## Estrutura

```
biologia/
├── data/                    # Dados genômicos (FASTA, GenBank)
├── src/                     # Código fonte
│   ├── trabalho1/           # Análise de Palíndromos
│   └── trabalho2/           # Detecção de Grampos
├── results/                 # Resultados e relatórios
└── scripts/                 # Scripts utilitários
```

## Projetos

### 1. Análise de Palíndromos Maximais
**Localização:** `src/trabalho1/`

- Identifica palíndromos maximais em sequências de DNA
- Organismo: Maribacter sp. HTCC2170 (CP002157.1)
- Mapeia para enzimas de restrição
- Analisa sobreposição com genes (CDS/ORF)

**Executar:**
```bash
cd src/trabalho1
python bacter_final.py
```

### 2. Detecção de Grampos (Hairpins)
**Localização:** `src/trabalho2/`

- Identifica estruturas de grampos em DNA
- Tamanho configurável (K)
- Exporta resultados para CSV

**Executar:**
```bash
cd src/trabalho2
python grampos.py
```

## Configuração

1. **Instalar dependências:**
```bash
pip install -r requirements.txt
```

2. **Configurar ambiente:**
```bash
chmod +x scripts/setup_environment.sh
./scripts/setup_environment.sh
```

## Dados

- **Organismo:** Maribacter sp. HTCC2170
- **Acesso NCBI:** CP002157.1
- **Tamanho:** 3,868,304 bp
- **Arquivos:** `data/maribacter_HTCC2170.fasta` e `data/maribacter_HTCC2170.gb`

## Resultados

- **Palíndromos:** Maior encontrado: `AAAATATTTT` (10 bp)
- **Enzimas:** EcoRI, HindIII, NaeI, NheI
- **Grampos:** Detectados com K=6, tamanho 12-20 bases

---
*Projetos acadêmicos de Biologia Computacional*# biocomp
