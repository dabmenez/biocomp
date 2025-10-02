#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
bacter_final.py

Análise Computacional de Palíndromos Maximais e Características Genômicas
em Maribacter sp. HTCC2170 (CP002157).

Este programa realiza as seguintes tarefas:
1. Obtém o registo genômico completo de Maribacter sp. HTCC2170 do NCBI.
2. Identifica todos os palíndromos maximais de um tamanho k específico em dois loci genômicos.
3. Determina o maior palíndromo maximal encontrado nessas regiões.
4. Verifica se os loci contêm Sequências Codificantes (CDS).
5. Mapeia os palíndromos encontrados para sítios de reconhecimento de enzimas de restrição.

Uso:
    python bacter_final.py --k 6 --intervals 82583-83599 297449-299453
    python bacter_final.py --find-largest  # Para encontrar o maior palíndromo
"""

import sys
import argparse
import os
from collections import defaultdict
from Bio import SeqIO
from io import StringIO

# URLs removidas - programa agora usa arquivos locais

# Mapa de complemento para DNA
COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def rev_comp(s: str) -> str:
    """Retorna o complemento reverso de uma sequência de DNA."""
    return s.translate(COMP)[::-1]

def is_palindrome(s: str) -> bool:
    """Verifica se uma sequência é um palíndromo (igual ao seu complemento reverso)."""
    return s.upper() == rev_comp(s).upper()

def load_genome_data():
    """Carrega os dados do genoma dos arquivos locais."""
    print("Carregando dados do genoma Maribacter sp. HTCC2170 dos arquivos locais...")
    
    try:
        # Carregar arquivo FASTA
        if not os.path.exists("maribacter_HTCC2170.fasta"):
            print("Erro: Arquivo maribacter_HTCC2170.fasta não encontrado!")
            print("Execute o programa uma vez para baixar os arquivos.")
            sys.exit(1)
        
        with open("maribacter_HTCC2170.fasta", "r", encoding="utf-8") as f:
            fasta_record = next(SeqIO.parse(f, "fasta"))
        
        # Carregar arquivo GenBank
        if not os.path.exists("maribacter_HTCC2170.gb"):
            print("Erro: Arquivo maribacter_HTCC2170.gb não encontrado!")
            print("Execute o programa uma vez para baixar os arquivos.")
            sys.exit(1)
        
        with open("maribacter_HTCC2170.gb", "r", encoding="utf-8") as f:
            gb_record = next(SeqIO.parse(f, "genbank"))
        
        sequence = str(fasta_record.seq)
        print(f"Genoma carregado: {fasta_record.id}, comprimento {len(sequence):,} bp")
        return sequence, gb_record
        
    except Exception as e:
        print(f"Erro ao carregar dados do genoma: {e}")
        print("Certifique-se de que os arquivos maribacter_HTCC2170.fasta e maribacter_HTCC2170.gb estão presentes.")
        sys.exit(1)

def find_maximal_palindromes_of_length_k(seq, k):
    """
    Encontra todos os palíndromos maximais de tamanho exato k.
    
    Args:
        seq (str): Sequência de DNA
        k (int): Tamanho desejado dos palíndromos
        
    Returns:
        dict: {sequência_palíndromo: [posições_início_1based]}
    """
    n = len(seq)
    found_positions = defaultdict(list)
    
    for i in range(0, n - k + 1):
        sub = seq[i:i+k]
        if is_palindrome(sub):
            # Verificar se é maximal (não pode ser estendido)
            can_extend = False
            if i > 0 and i + k < n:
                # Tentar estender para a esquerda e direita
                extended = seq[i-1:i+k+1]
                if is_palindrome(extended):
                    can_extend = True
            
            if not can_extend:
                found_positions[sub.upper()].append(i + 1)  # 1-based
    
    return found_positions

def find_all_maximal_palindromes(seq):
    """
    Encontra todos os palíndromos maximais de qualquer tamanho.
    
    Args:
        seq (str): Sequência de DNA
        
    Returns:
        list: Lista de tuplas (início, fim, sequência) em coordenadas 0-based
    """
    palindromes = []
    n = len(seq)
    
    # Palíndromos de comprimento ímpar (centro único)
    for i in range(n):
        l, r = i, i
        while l >= 0 and r < n and seq[l].upper() == rev_comp(seq[r]).upper():
            l -= 1
            r += 1
        # O palíndromo maximal é de l+1 a r-1
        if r - l > 1:  # Pelo menos 2 bases
            palindromes.append((l + 1, r, seq[l+1:r]))
    
    # Palíndromos de comprimento par (centro entre duas bases)
    for i in range(n - 1):
        l, r = i, i + 1
        while l >= 0 and r < n and seq[l].upper() == rev_comp(seq[r]).upper():
            l -= 1
            r += 1
        if r - l > 2:  # Pelo menos 4 bases
            palindromes.append((l + 1, r, seq[l+1:r]))
    
    # Remover duplicatas e ordenar
    return sorted(list(set(palindromes)))

def check_cds_overlap(gb_record, start, end):
    """
    Verifica se uma região genômica se sobrepõe a alguma CDS anotada.
    
    Args:
        gb_record: Registro GenBank
        start (int): Posição inicial (1-based)
        end (int): Posição final (1-based)
        
    Returns:
        list: Lista de informações sobre CDS encontradas
    """
    cds_found = []
    
    for feature in gb_record.features:
        if feature.type == "CDS":
            # Coordenadas da feature (0-based no Biopython)
            f_start = int(feature.location.start) + 1  # Converter para 1-based
            f_end = int(feature.location.end)
            
            # Verificar sobreposição
            if not (f_end < start or f_start > end):
                locus_tag = feature.qualifiers.get("locus_tag", ["N/A"])[0]
                product = feature.qualifiers.get("product", ["N/A"])[0]
                cds_found.append({
                    'locus_tag': locus_tag,
                    'product': product,
                    'start': f_start,
                    'end': f_end
                })
    
    return cds_found

def map_to_restriction_enzymes(palindromes):
    """
    Mapeia sequências palindrômicas para enzimas de restrição conhecidas.
    
    Args:
        palindromes (list): Lista de sequências palindrômicas
        
    Returns:
        dict: Mapeamento de palíndromos para informações da enzima
    """
    # Base de dados de enzimas de restrição (REBASE simplificada)
    restriction_enzymes = {
        "AAGCTT": ("HindIII", "Haemophilus influenzae Rd"),
        "CTGCAG": ("PstI", "Providencia stuartii"),
        "GAATTC": ("EcoRI", "Escherichia coli R"),
        "GGATCC": ("BamHI", "Bacillus amyloliquefaciens H"),
        "GTCGAC": ("SalI", "Streptomyces albus G"),
        "GTATAC": ("BstEII", "Bacillus stearothermophilus EII"),
        "GGCGCC": ("NarI", "Nocardia argentinensis"),
        "CATATG": ("NdeI", "Neisseria denitrificans"),
        "CCCGGG": ("SmaI", "Serratia marcescens"),
        "GCGGCCGC": ("NotI", "Nocardia otitidiscaviarum"),
        "TCTAGA": ("XbaI", "Xanthomonas badrii"),
        "GCTAGC": ("NheI", "Neisseria mucosa"),
        "GGTACC": ("KpnI", "Klebsiella pneumoniae"),
        "GAGCTC": ("SacI", "Streptomyces achromogenes"),
        "AGATCT": ("BglII", "Bacillus globigii"),
        "AAGCTT": ("HindIII", "Haemophilus influenzae Rd"),
        "TTAATTAA": ("AseI", "Aquifex aeolicus"),
        "GCCGGC": ("NaeI", "Nocardia aerocolonigenes"),
    }
    
    enzyme_matches = {}
    unique_palindromes = set(palindromes)
    
    for pal in unique_palindromes:
        if pal in restriction_enzymes:
            enzyme_matches[pal] = restriction_enzymes[pal]
    
    return enzyme_matches

def analyze_region(sequence, gb_record, start, end, k=None):
    """
    Analisa uma região específica do genoma.
    
    Args:
        sequence (str): Sequência completa do genoma
        gb_record: Registro GenBank
        start (int): Posição inicial (1-based)
        end (int): Posição final (1-based)
        k (int, optional): Tamanho específico de palíndromos a buscar
    """
    print(f"\n{'='*60}")
    print(f"ANÁLISE DA REGIÃO {start}..{end}")
    print(f"{'='*60}")
    
    # Extrair subsequência
    subseq = sequence[start-1:end]
    print(f"Comprimento da região: {len(subseq)} bp")
    
    # Verificar CDS
    print(f"\n--- ANÁLISE DE CDS/ORF ---")
    cds_results = check_cds_overlap(gb_record, start, end)
    if cds_results:
        print("✓ REGIÃO CONTÉM CDS/ORF ANOTADOS:")
        for cds in cds_results:
            print(f"  • Locus: {cds['locus_tag']}")
            print(f"    Produto: {cds['product']}")
            print(f"    Localização: {cds['start']}..{cds['end']}")
    else:
        print("✗ Nenhuma CDS/ORF anotada encontrada nesta região")
    
    # Buscar palíndromos
    if k is not None:
        print(f"\n--- PALÍNDROMOS MAXIMAIS DE TAMANHO {k} ---")
        palindromes = find_maximal_palindromes_of_length_k(subseq, k)
        
        if palindromes:
            print(f"Encontrados {len(palindromes)} palíndromos únicos:")
            for pal, positions in palindromes.items():
                # Converter posições locais para globais
                global_positions = [start - 1 + pos for pos in positions]
                print(f"  • {pal} => posições no genoma: {global_positions}")
        else:
            print(f"Nenhum palíndromo maximal de tamanho {k} encontrado")
    
    # Buscar todos os palíndromos maximais
    print(f"\n--- TODOS OS PALÍNDROMOS MAXIMAIS ---")
    all_palindromes = find_all_maximal_palindromes(subseq)
    
    if all_palindromes:
        # Agrupar por tamanho
        by_size = defaultdict(list)
        for start_pos, end_pos, seq in all_palindromes:
            by_size[len(seq)].append((start_pos, end_pos, seq))
        
        print(f"Encontrados {len(all_palindromes)} palíndromos maximais:")
        for size in sorted(by_size.keys(), reverse=True):
            palindromes_of_size = by_size[size]
            print(f"  Tamanho {size}: {len(palindromes_of_size)} palíndromos")
            
            # Mostrar apenas os primeiros 5 de cada tamanho
            for i, (start_pos, end_pos, seq) in enumerate(palindromes_of_size[:5]):
                global_start = start + start_pos
                global_end = start + end_pos - 1
                print(f"    • {seq} (posição {global_start}..{global_end})")
            
            if len(palindromes_of_size) > 5:
                print(f"    ... e mais {len(palindromes_of_size) - 5} palíndromos")
        
        # Encontrar o maior
        largest = max(all_palindromes, key=lambda x: len(x[2]))
        largest_global_start = start + largest[0]
        largest_global_end = start + largest[1] - 1
        print(f"\nMAIOR PALÍNDROMO MAXIMAL:")
        print(f"  Sequência: {largest[2]}")
        print(f"  Tamanho: {len(largest[2])} bp")
        print(f"  Posição: {largest_global_start}..{largest_global_end}")
    
    # Mapear para enzimas de restrição
    print(f"\n--- ENZIMAS DE RESTRIÇÃO ---")
    if all_palindromes:
        all_sequences = [seq for _, _, seq in all_palindromes]
        enzyme_matches = map_to_restriction_enzymes(all_sequences)
        
        if enzyme_matches:
            print(f"Encontradas {len(enzyme_matches)} sequências que correspondem a enzimas de restrição:")
            count = 0
            for pal, (enzyme, organism) in enzyme_matches.items():
                if count < 4:  # Mostrar apenas as primeiras 4
                    print(f"  • {pal} => {enzyme} (origem: {organism})")
                    count += 1
            if len(enzyme_matches) > 4:
                print(f"  ... e mais {len(enzyme_matches) - 4} enzimas")
        else:
            print("Nenhuma sequência corresponde a enzimas de restrição conhecidas")

def find_largest_palindrome(sequence, gb_record, regions):
    """
    Encontra o maior palíndromo maximal em todas as regiões especificadas.
    """
    print(f"\n{'='*60}")
    print("BUSCA PELO MAIOR PALÍNDROMO MAXIMAL")
    print(f"{'='*60}")
    
    all_palindromes = []
    
    for start, end in regions:
        subseq = sequence[start-1:end]
        palindromes = find_all_maximal_palindromes(subseq)
        
        # Converter para coordenadas globais
        for pal_start, pal_end, seq in palindromes:
            global_start = start + pal_start
            global_end = start + pal_end - 1
            all_palindromes.append((global_start, global_end, seq))
    
    if all_palindromes:
        # Encontrar o maior
        largest = max(all_palindromes, key=lambda x: len(x[2]))
        
        print(f"MAIOR PALÍNDROMO MAXIMAL ENCONTRADO:")
        print(f"  Sequência: {largest[2]}")
        print(f"  Tamanho: {len(largest[2])} bp")
        print(f"  Posição: {largest[0]}..{largest[1]}")
        
        # Verificar em qual região está
        for i, (start, end) in enumerate(regions):
            if start <= largest[0] <= end:
                print(f"  Região: {i+1} ({start}..{end})")
                break
        
        # Mapear para enzimas de restrição
        enzyme_matches = map_to_restriction_enzymes([largest[2]])
        if enzyme_matches:
            pal, (enzyme, organism) = next(iter(enzyme_matches.items()))
            print(f"  Enzima de restrição: {enzyme} (origem: {organism})")
        else:
            print(f"  Não corresponde a nenhuma enzima de restrição conhecida")
    else:
        print("Nenhum palíndromo maximal encontrado nas regiões especificadas")

def generate_report(sequence, gb_record, regions):
    """
    Gera um relatório completo em Markdown com todas as análises.
    """
    report = []
    
    # Cabeçalho do relatório
    report.append("# Relatório de Análise de Palíndromos - Maribacter sp. HTCC2170")
    report.append("")
    report.append(f"**Genoma:** CP002157.1 ({len(sequence):,} bp)")
    report.append("")
    report.append("Este relatório apresenta os resultados da análise de palíndromos maximais em duas regiões específicas do genoma da bactéria Maribacter sp. HTCC2170.")
    report.append("")
    report.append("---")
    report.append("")
    
    # Análise de cada região
    all_palindromes = []
    all_restriction_enzymes = set()
    
    for i, (start, end) in enumerate(regions):
        report.append(f"## Região {i+1}: {start}..{end}")
        report.append("")
        
        # Informações básicas
        subseq = sequence[start-1:end]
        report.append(f"Esta região possui {len(subseq)} pares de bases.")
        report.append("")
        
        # Análise de CDS/ORF
        report.append("### Genes Encontrados")
        cds_results = check_cds_overlap(gb_record, start, end)
        if cds_results:
            report.append("A região contém os seguintes genes:")
            report.append("")
            for cds in cds_results:
                report.append(f"- **{cds['locus_tag']}**: {cds['product']} (posições {cds['start']}..{cds['end']})")
            report.append("")
        else:
            report.append("Nenhum gene foi encontrado nesta região.")
            report.append("")
        
        # Palíndromos de diferentes tamanhos
        report.append("### Palíndromos Encontrados")
        report.append("")
        
        # Testar diferentes tamanhos de k
        for k in [4, 6, 8, 10, 12, 14, 16, 18, 20]:
            palindromes_k = find_maximal_palindromes_of_length_k(subseq, k)
            if palindromes_k:
                report.append(f"**Palíndromos de {k} bases:** {len(palindromes_k)} sequências diferentes")
                for pal, positions in list(palindromes_k.items())[:3]:  # Mostrar apenas os primeiros 3
                    global_positions = [start - 1 + pos for pos in positions]
                    report.append(f"- {pal} (posições: {global_positions})")
                if len(palindromes_k) > 3:
                    report.append(f"- ... e mais {len(palindromes_k) - 3} sequências")
                report.append("")
        
        # Todos os palíndromos maximais
        all_pals = find_all_maximal_palindromes(subseq)
        if all_pals:
            # Agrupar por tamanho
            by_size = defaultdict(list)
            for start_pos, end_pos, seq in all_pals:
                by_size[len(seq)].append((start_pos, end_pos, seq))
            
            report.append("### Resumo Geral")
            report.append("")
            report.append("Distribuição dos palíndromos por tamanho:")
            for size in sorted(by_size.keys(), reverse=True):
                palindromes_of_size = by_size[size]
                report.append(f"- {size} bases: {len(palindromes_of_size)} sequências")
            
            # Encontrar o maior
            largest = max(all_pals, key=lambda x: len(x[2]))
            largest_global_start = start + largest[0]
            largest_global_end = start + largest[1] - 1
            report.append("")
            report.append(f"**Maior palíndromo encontrado:**")
            report.append(f"- Sequência: {largest[2]}")
            report.append(f"- Tamanho: {len(largest[2])} bases")
            report.append(f"- Localização: {largest_global_start}..{largest_global_end}")
            report.append("")
            
            # Coletar para análise geral
            for _, _, seq in all_pals:
                all_palindromes.append(seq)
        
        # Enzimas de restrição
        report.append("### Sítios de Enzimas de Restrição")
        report.append("")
        if all_pals:
            all_sequences = [seq for _, _, seq in all_pals]
            enzyme_matches = map_to_restriction_enzymes(all_sequences)
            
            if enzyme_matches:
                report.append("Alguns palíndromos correspondem a sítios de enzimas de restrição conhecidas:")
                report.append("")
                for pal, (enzyme, organism) in enzyme_matches.items():
                    report.append(f"- {pal} → {enzyme} (de {organism})")
                    all_restriction_enzymes.add((enzyme, pal, organism))
                report.append("")
            else:
                report.append("Nenhum palíndromo corresponde a enzimas de restrição conhecidas.")
                report.append("")
        
        report.append("---")
        report.append("")
    
    # Análise geral
    report.append("## Resultados Principais")
    report.append("")
    
    # Maior palíndromo de todas as regiões
    if all_palindromes:
        largest_overall = max(all_palindromes, key=len)
        report.append("### Maior Palíndromo Identificado")
        report.append("")
        report.append(f"O maior palíndromo encontrado em ambas as regiões é:")
        report.append(f"- Sequência: {largest_overall}")
        report.append(f"- Tamanho: {len(largest_overall)} bases")
        report.append("")
        
        # Verificar se corresponde a enzima de restrição
        enzyme_matches = map_to_restriction_enzymes([largest_overall])
        if enzyme_matches:
            pal, (enzyme, organism) = next(iter(enzyme_matches.items()))
            report.append(f"Este palíndromo corresponde ao sítio da enzima {enzyme}.")
        else:
            report.append("Este palíndromo não corresponde a nenhuma enzima de restrição conhecida.")
        report.append("")
    
    # Resumo das enzimas encontradas
    if all_restriction_enzymes:
        report.append("### Enzimas de Restrição Encontradas")
        report.append("")
        report.append(f"Total de {len(all_restriction_enzymes)} enzimas diferentes foram identificadas:")
        report.append("")
        for enzyme, site, organism in sorted(all_restriction_enzymes):
            report.append(f"- {enzyme} (sítio: {site}) - {organism}")
        report.append("")
    
    # Respostas às perguntas específicas
    report.append("## Respostas às Perguntas Solicitadas")
    report.append("")
    
    report.append("### 1. Palíndromos maximais de tamanho k")
    report.append("")
    report.append("O programa foi testado com diferentes valores de k (4, 6, 8, 10, 12, 14, 16, 18, 20).")
    report.append("")
    report.append("**Resultados para k=6:**")
    for i, (start, end) in enumerate(regions):
        subseq = sequence[start-1:end]
        palindromes_k6 = find_maximal_palindromes_of_length_k(subseq, 6)
        if palindromes_k6:
            report.append(f"- Região {i+1} ({start}..{end}): {len(palindromes_k6)} sequências diferentes")
            for pal, positions in palindromes_k6.items():
                global_positions = [start - 1 + pos for pos in positions]
                report.append(f"  - {pal} (posições: {global_positions})")
        else:
            report.append(f"- Região {i+1} ({start}..{end}): Nenhum palíndromo de 6 bases")
    report.append("")
    
    report.append("### 2. Maior palíndromo maximal encontrado")
    report.append("")
    if all_palindromes:
        largest = max(all_palindromes, key=len)
        report.append(f"O maior palíndromo maximal encontrado é: {largest}")
        report.append(f"- Tamanho: {len(largest)} bases")
        report.append("")
        report.append("Teste de maximalidade (verificando até que tamanho existem palíndromos):")
        for k in range(2, 22, 2):
            found_any = False
            for start, end in regions:
                subseq = sequence[start-1:end]
                pals = find_maximal_palindromes_of_length_k(subseq, k)
                if pals:
                    found_any = True
                    break
            if found_any:
                report.append(f"- {k} bases: Encontrados")
            else:
                report.append(f"- {k} bases: Nenhum encontrado")
                break
    report.append("")
    
    report.append("### 3. Enzimas de restrição encontradas")
    report.append("")
    if all_restriction_enzymes:
        report.append(f"Foram encontradas {len(all_restriction_enzymes)} enzimas de restrição diferentes:")
        report.append("")
        for i, (enzyme, site, organism) in enumerate(sorted(all_restriction_enzymes), 1):
            report.append(f"{i}. {enzyme} (sítio: {site}) - {organism}")
    else:
        report.append("Nenhuma enzima de restrição foi encontrada.")
    report.append("")
    
    report.append("### 4. Análise de CDS/ORF")
    report.append("")
    report.append("Ambos os trechos analisados contêm genes anotados:")
    report.append("")
    for i, (start, end) in enumerate(regions):
        cds_results = check_cds_overlap(gb_record, start, end)
        if cds_results:
            report.append(f"**Região {i+1} ({start}..{end}):**")
            for cds in cds_results:
                report.append(f"- {cds['locus_tag']}: {cds['product']}")
        else:
            report.append(f"**Região {i+1} ({start}..{end}):** Nenhum gene encontrado")
    report.append("")
    
    report.append("---")
    report.append("")
    
    return "\n".join(report)

def main():
    parser = argparse.ArgumentParser(
        description="Análise de palíndromos maximais no genoma Maribacter sp. HTCC2170",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  python bacter_final.py                    # Executa análise completa e gera relatório
  python bacter_final.py --k 6 --intervals 82583-83599 297449-299453
  python bacter_final.py --find-largest --intervals 82583-83599 297449-299453
        """
    )
    
    parser.add_argument("--k", type=int, help="Tamanho específico dos palíndromos a buscar")
    parser.add_argument("--intervals", nargs="+",
                        help="Intervalos a analisar no formato start-end (ex: 82583-83599)")
    parser.add_argument("--find-largest", action="store_true",
                        help="Encontrar o maior palíndromo maximal em todas as regiões")
    parser.add_argument("--generate-report", action="store_true", default=True,
                        help="Gerar relatório completo em Markdown (padrão)")
    
    args = parser.parse_args()
    
    # Se não especificou argumentos, executar análise completa
    if not args.k and not args.find_largest and not args.intervals:
        print("Executando análise completa da tarefa...")
        regions = [(82583, 83599), (297449, 299453)]
    else:
        # Validar argumentos
        if not args.intervals:
            print("Erro: Especifique --intervals")
            sys.exit(1)
        
        if args.k and args.k % 2 != 0:
            print("Erro: k deve ser um número par")
            sys.exit(1)
        
        # Parse dos intervalos
        regions = []
        for interval in args.intervals:
            if "-" not in interval:
                print(f"Erro: Formato de intervalo inválido: {interval}")
                sys.exit(1)
            
            try:
                start, end = map(int, interval.split("-"))
                regions.append((start, end))
            except ValueError:
                print(f"Erro: Intervalo inválido: {interval}")
                sys.exit(1)
    
    # Carregar dados do genoma
    sequence, gb_record = load_genome_data()
    
    # Verificar se os intervalos estão dentro do genoma
    genome_length = len(sequence)
    for start, end in regions:
        if start < 1 or end > genome_length or start > end:
            print(f"Erro: Intervalo {start}-{end} está fora dos limites do genoma (1-{genome_length})")
            sys.exit(1)
    
    print(f"\nAnálise de {len(regions)} região(ões) do genoma Maribacter sp. HTCC2170")
    print(f"Tamanho do genoma: {genome_length:,} bp")
    
    # Se solicitou análise específica, executar
    if args.k or args.find_largest:
        # Analisar cada região
        for i, (start, end) in enumerate(regions):
            analyze_region(sequence, gb_record, start, end, args.k)
        
        # Se solicitado, encontrar o maior palíndromo
        if args.find_largest:
            find_largest_palindrome(sequence, gb_record, regions)
        
        print(f"\n{'='*60}")
        print("ANÁLISE CONCLUÍDA")
        print(f"{'='*60}")
    
    # Gerar relatório completo
    if args.generate_report:
        print("\nGerando relatório completo...")
        report = generate_report(sequence, gb_record, regions)
        
        # Salvar relatório
        with open("relatorio_palindromos_maribacter.md", "w", encoding="utf-8") as f:
            f.write(report)
        
        print("Relatório salvo em: relatorio_palindromos_maribacter.md")
        print("\nResumo do relatório:")
        print("- Análise completa de palíndromos maximais")
        print("- Verificação de CDS/ORF nos trechos")
        print("- Identificação de enzimas de restrição")
        print("- Respostas a todas as perguntas da tarefa")

if __name__ == "__main__":
    main()
